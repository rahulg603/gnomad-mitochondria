import argparse
import logging
import math
import os
import re
import sys
import dxpy
import pyspark

import hail as hl

from os.path import dirname
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
            if min_partitions > 10:
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
    # start SQL session and initialize constants
    my_database = dxpy.find_one_data_object(name=args.dx_init.lower())["id"]
    input_ht = f'dnax://{my_database}/{args.input_ht}/'
    output_ht = f'dnax://{my_database}/{args.output_ht}'
    temp_dir = f'dnax://{my_database}/{args.temp_dir}/'
    chunk_size = args.chunk_size
    overwrite = args.overwrite
    keep_targets = args.keep_targets
    sc = pyspark.SparkContext()
    spark = pyspark.sql.SparkSession(sc)
    hl.init(sc=sc, tmp_dir=temp_dir)
    hl._set_flags(no_whole_stage_codegen='1')

    if args.overwrite == False and hl.hadoop_exists(output_ht):
        logger.warning(
            "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
            output_ht,
        )
    # Ensure that user supplied ht extension for output_ht
    if not output_ht.endswith(".ht"):
        sys.exit("Path supplied as output_ht must end with .ht extension")

    mt_list = []
    logger.info(
        "Reading in individual coverage files as matrix tables and adding to a list of matrix tables..."
    )
    idx = 0
    paths = hl.read_table(input_ht)
    pairs_for_coverage = paths.annotate(pairs = (paths.batch, paths.coverage)).pairs.collect()
    for batch, base_level_coverage_metrics in pairs_for_coverage:
        idx+=1
        mt = hl.import_matrix_table(
            'file://' + base_level_coverage_metrics,
            delimiter="\t",
            row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
            row_key=["chrom", "pos"],
            min_partitions=args.n_read_partitions,
        )
        if not keep_targets:
            mt = mt.drop("target")
        else:
            mt = mt.key_rows_by(*["chrom", "pos", "target"])
        mt = mt.key_cols_by().rename({"x": "coverage", 'col_id':'s'}).key_cols_by('s')
        mt = mt.annotate_cols(batch = batch)

        mt_list.append(mt)
        if idx % 10 == 0:
            logger.info(f"Imported batch {str(idx)}...")

    logger.info("Joining individual coverage mts...")
    cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size, min_partitions=args.n_read_partitions)
    n_samples = cov_mt.count_cols()

    logger.info("Adding coverage annotations...")
    # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
    cov_mt = cov_mt.annotate_rows(
        locus=hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"),
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
        over_1000=hl.float((hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)),
    )
    cov_mt.show()

    cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")

    output_mt = re.sub(r"\.ht$", ".mt", output_ht)
    output_tsv = re.sub(r"\.ht$", ".tsv", output_ht)
    output_samples = re.sub(r"\.ht$", "_sample_level.txt", output_ht)

    if not args.hail_only:
        logger.info("Writing sample level coverage...")
        sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
        sample_mt.coverage.export(output_samples)

    logger.info("Writing coverage mt and ht...")
    cov_mt.write(output_mt, overwrite=overwrite)
    cov_ht = cov_mt.rows()
    cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)
    if not args.hail_only:
        cov_ht.export(output_tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script combines individual mitochondria coverage files and outputs a hail table with coverage annotations"
    )
    parser.add_argument(
        "-i",
        "--input-ht",
        help="Input ht with paths to coverage file.",
        required=True,
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
        "--chunk-size",
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    parser.add_argument(
        "--keep-targets", help="Will add an annotation for target from the coverage file", action="store_true"
    )
    parser.add_argument(
        "--n-read-partitions", type=int, help="The number of partitions to use when reading tsvs. This should be 1 if the files are small.", default=1
    )
    parser.add_argument(
        "--hail-only", action='store_true', help='Skip generating flat files.'
    )
    parser.add_argument(
        "--dx-init", type=str, required=True, help='SQL database path for use in DNAnexus.'
    )

    args = parser.parse_args()
    main(args)
