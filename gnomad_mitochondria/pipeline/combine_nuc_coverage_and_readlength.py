import argparse
import logging
import math
import os
import re
import sys

import hail as hl

from os.path import dirname
from gnomad.utils.slack import slack_notifications
from hail.utils.java import info

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Annotate coverage")
logger.setLevel(logging.INFO)


def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int, min_partitions: int) -> hl.MatrixTable:
    """
    Hierarchically join together MatrixTables in the provided list.

    :param mts: List of MatrixTables to join together
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable
    """
    # Convert the MatrixTables to tables where entries are an array of structs
    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        # Calculate the number of jobs to run based on the chunk size
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []

        for i in range(n_jobs):
            # Grab just the tables for the given job
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )

            # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            merged = merged.checkpoint(f"stage_{stage}_job_{i}_pre.ht", overwrite=True)
            # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )

            # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            next_stage.append(
                merged.checkpoint(
                    os.path.join(temp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"Completed stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


def main(args):  # noqa: D103
    input_tsv = args.input_tsv
    idxstats = args.idxstats
    read_len = args.read_length
    output_ht = args.output_ht
    temp_dir = args.temp_dir
    chunk_size = args.chunk_size
    overwrite = args.overwrite

    if args.overwrite == False and hl.hadoop_exists(output_ht):
        logger.warning(
            "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
            output_ht,
        )
    # Ensure that user supplied ht extension for output_ht
    if not output_ht.endswith(".ht"):
        sys.exit("Path supplied as output_ht must end with .ht extension")

    logger.info(
        "Reading in individual coverage files as matrix tables and adding to a list of matrix tables..."
    )

    ht = hl.import_table(input_tsv)
    if args.head is not None:
        ht = ht.head(args.head)
    ht = ht.filter(ht[idxstats] != "")
    ht = ht.filter(ht[read_len] != "")
    idxstats_load = {}
    lens_load = {}
    df = ht.to_pandas()
    for _, row in df.iterrows():
        idxstats_load[row["s"]] = row[idxstats]
        lens_load[row["s"]] = row[read_len]

    mt_list = []
    ht_list = []
    for sample in idxstats_load.keys():
        path_idx = idxstats_load[sample]
        path_lens = lens_load[sample]

        ht_idx = hl.import_table(path_idx, no_header=True, impute=True, min_partitions=args.n_read_partitions
                                ).rename({'f0': 'chrom', 'f1': 'chrom_len', 'f2':'mapped_reads', 'f3': 'unmampped_reads'}
                                ).key_by('chrom')
        ht_idx = ht_idx.annotate(s=sample)
        mt_idx = ht_idx.to_matrix_table(row_key=['chrom'], col_key=['s'], row_fields=['chrom_len'])#.persist()
        mt_list.append(mt_idx)

        ht_len = hl.import_table(path_lens, impute=True).annotate(s=sample).key_by('s')
        ht_list.append(ht_len)

    if len(ht_list) > 1:
        ht_read_lens = ht_list[0].union(*ht_list[1:len(ht_list)])
    else:
        ht_read_lens = ht_list[0]
    ht_read_lens = ht_read_lens.checkpoint(os.path.join(temp_dir, f"read_lengths_temp.ht"), overwrite=True)
    logger.info("Joining individual coverage mts...")

    cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=args.n_read_partitions)
    n_samples = cov_mt.count_cols()
    n_samples_len = ht_read_lens.count()

    if n_samples != n_samples_len:
        raise ValueError("ERROR: number of cols in coverage MT must equal number of rows in read length mt")

    logger.info("Adding coverage annotations...")
    # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
    cov_mt = cov_mt.annotate_cols(**ht_read_lens[cov_mt.col_key])
    cov_mt = cov_mt.annotate_rows(chrom_len = mt_list[0].rows()[cov_mt.row_key].chrom_len)
    cov_mt.show()
    cov_mt.describe()

    output_mt = re.sub(r"\.ht$", ".mt", output_ht)
    output_tsv_cov = re.sub(r"\.ht$", "_mapped_read_counts.tsv", output_ht)
    output_filtered = re.sub(r"\.ht$", "_nuclearcoverages.tsv", output_ht)

    if args.filter_to_main_chrom:
        filter_to = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
        if args.include_mtdna_cov:
            filter_to = filter_to + ['chrM']
    else:
        filter_to = list(set(cov_mt_filt.chrom.collect()) - set(['chrM']))
        if args.include_mtdna_cov:
            filter_to = filter_to + ['chrM']

    if not args.hail_only:
        logger.info("Writing per_sample information...")
        cov_mt_filt = cov_mt.filter_rows(hl.literal(filter_to).contains(cov_mt.chrom))
        cov_mt_filt = cov_mt_filt.annotate_cols(total_mapped_reads = hl.agg.sum(cov_mt_filt.mapped_reads))
        cov_ht_filt = cov_mt_filt.annotate_cols(nuclear_coverage = cov_mt_filt.total_mapped_reads * cov_mt_filt.READ_LENGTH / hl.agg.sum(cov_mt_filt.chrom_len)).cols()
        cov_ht_filt.export(output_filtered)

    logger.info("Writing coverage mt and ht...")
    cov_mt.write(output_mt, overwrite=overwrite)
    cov_ht = cov_mt.cols()
    cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)
    if not args.hail_only:
        cov_mt.mapped_reads.export(output_tsv_cov)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script combines individual mitochondria coverage files and outputs a hail table with coverage annotations"
    )
    parser.add_argument(
        "-i",
        "--input-tsv",
        help="Input file with files to combine. Expects sample name to be 's'.",
        required=True,
    )
    parser.add_argument(
        "--idxstats",
        help="Column containing idxstats data.",
        required=True,
    )
    parser.add_argument(
        "--read-length",
        help="Column containing read length data."
    )
    parser.add_argument(
        "--head", default=None, type=int,
        help="Max number of records to process."
    )
    parser.add_argument(
        '--filter-to-main-chrom', action='store_true', help="If enabled, will filter to main chromosomes only."
    )
    parser.add_argument(
        "--include-mtdna-cov", action='store_true', help='If enabled will include mtDNA in genome coverage estimates.'
    )
    parser.add_argument(
        "-o", "--output-ht", help="Name of ht to write output", required=True
    )
    parser.add_argument(
        "-t",
        "--temp-dir",
        help="Temporary directory to use for intermediate outputs",
        required=True,
    )
    parser.add_argument(
        "--slack-token", help="Slack token that allows integration with slack",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to",
    )
    parser.add_argument(
        "--chunk-size",
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    parser.add_argument(
        "--n-read-partitions", type=int, help="The number of partitions to use when reading tsvs. This should be 1 if the files are small.", default=1
    )
    parser.add_argument(
        "--hail-only", action='store_true', help='Skip generating flat files.'
    )

    args = parser.parse_args()

    # Both a slack token and slack channel must be supplied to receive notifications on slack
    if args.slack_channel and args.slack_token:
        with slack_notifications(args.slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
