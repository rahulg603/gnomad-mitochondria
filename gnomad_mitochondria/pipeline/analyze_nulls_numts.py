import argparse
import hail as hl

def prune_by_intervals(mt, bases):
    ht = mt.rows()
    ht_locs = ht.group_by(ht.target).aggregate(min_pos = hl.agg.min(ht.locus.position) + bases, 
                                               max_pos = hl.agg.max(ht.locus.position) - bases)
    
    if ht_locs.filter(ht_locs.min_pos > ht_locs.max_pos).count() > 0:
        raise ValueError('ERROR: there are intervals where the number of bases to drop is greater than the total interval size.')
    
    ht_annot = ht.annotate(**ht_locs[ht.target])
    ht_annot = ht_annot.filter((ht_annot.locus.position >= ht_annot.min_pos) & (ht_annot.locus.position <= ht_annot.max_pos))
    mt = mt.semi_join_rows(ht_annot)
    return mt


def annotate_by_pois(mt):
    # Determines if a sample has coverage that is too high or low using a Poisson distribution
    mt = mt.annotate_entries(p_hi = hl.ppois(mt.coverage, mt.mean_coverage, lower_tail=False),
                             p_lo = hl.ppois(mt.coverage, mt.mean_coverage, lower_tail=True))
    mt = mt.annotate_entries(is_hi = mt.p_hi <= 0.025, is_lo = mt.p_lo <= 0.025)
    return mt


def get_per_target_stats(mt):
    mt = mt.group_rows_by(mt.target
          ).aggregate_rows(N = hl.agg.count()
          ).aggregate_entries(N_hi = hl.agg.count_where(mt.is_hi), 
                              N_lo = hl.agg.count_where(mt.is_lo)).result()
    mt = mt.annotate_entries(prop_hi = mt.N_hi/mt.N, prop_lo = mt.N_lo/mt.N)
    return mt


def get_nulls(mt_target, mt_null, tmp, n_pick=1000):
    """ Gets per-target null distributions.
    Outputs a MatrixTable with a vector of null values for each NUMT target and individual.
    """
    mt_null_agg = mt_null.group_rows_by(mt_null.target
                        ).aggregate_entries(null_vec_hi = hl.cumulative_sum(hl.agg.collect(hl.int32(mt_null.is_hi))),
                                            null_vec_lo = hl.cumulative_sum(hl.agg.collect(hl.int32(mt_null.is_lo)))).result()
    mt_null_agg = mt_null_agg.repartition(1000).checkpoint(tmp + 'tmp_null_1.mt', overwrite=True)
    mt_null_agg = mt_null_agg.annotate_rows(rand_idx = hl.int(mt_null_agg.target.split('_')[1]))
    mt_null_agg = mt_null_agg.filter_rows(mt_null_agg.rand_idx <= n_pick)
    ht_null_agg = mt_null_agg.annotate_cols(null_hi_s=hl.agg.collect(mt_null_agg.null_vec_hi),
                                            null_lo_s=hl.agg.collect(mt_null_agg.null_vec_lo)).cols()
    ht_null_agg = ht_null_agg.repartition(500).checkpoint(tmp + 'tmp_null_2.ht', overwrite=True)
                                            
    mt_target = mt_target.annotate_cols(null_hi_s = ht_null_agg[mt_target.s].null_hi_s,
                                        null_lo_s = ht_null_agg[mt_target.s].null_lo_s)
    mt_target = mt_target.annotate_cols(null_hi_s = hl.map(lambda x: hl.enumerate(x), mt_target.null_hi_s),
                                        null_lo_s = hl.map(lambda x: hl.enumerate(x), mt_target.null_lo_s))
    mt_target = mt_target.repartition(1000).checkpoint(tmp + 'tmp_null_3.mt', overwrite=True)

    def get_first(vec, n):
        return hl.filter(lambda z: z[0] < n, vec)
    
    mt_target = mt_target.annotate_entries(null_vec_hi = hl.map(lambda x: hl.filter(lambda q: q[1], get_first(x, mt_target.N)).length(),
                                                                mt_target.null_hi_s),
                                           null_vec_lo = hl.map(lambda x: hl.filter(lambda q: q[1], get_first(x, mt_target.N)).length(),
                                                                mt_target.null_lo_s)
                        ).select_cols()
    mt_target = mt_target.checkpoint(tmp + 'tmp_null_4.mt', overwrite=True)
    return mt_target


def main(args):  # noqa: D103
    tmp = args.tmp
    qc = args.qc
    bases = args.prune_interval
    nulls = args.nulls
    numts = args.numts
    numt_only = args.numt_only

    tmp = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/tmp/"
    qc = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/qc_result_sample.tsv"
    bases = 500
    nulls = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/test_null/annotated_file_coverage_null.mt"
    numts = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/test_null/annotated_file_coverage_numts.mt"

    tmp = "gs://ccdg-4day-temp/tmp/"
    numts = 'gs://ccdg/rgupta/testing_nulls/annotated_file_coverage_numts.mt'
    nulls = 'gs://ccdg/rgupta/testing_nulls/annotated_file_coverage_null.mt'
    qc = "gs://ccdg/rgupta/testing_nulls/qc_result_sample.tsv"

    ht_qc = hl.import_table(qc, impute=True)
    ht_qc = ht_qc.select(s = ht_qc['entity:qc_result_sample_id'], mean_coverage=ht_qc.mean_coverage).key_by('s')
    
    mt_numt = prune_by_intervals(hl.read_matrix_table(numts), bases)
    mt_numt = mt_numt.annotate_cols(mean_coverage = ht_qc[mt_numt.col_key].mean_coverage)
    mt_numt = annotate_by_pois(mt_numt)
    mt_numt_ct = get_per_target_stats(mt_numt)

    if not numt_only:
        mt_null = prune_by_intervals(hl.read_matrix_table(nulls), bases)
        mt_null = mt_null.annotate_cols(mean_coverage = ht_qc[mt_null.col_key].mean_coverage)
        mt_null = annotate_by_pois(mt_null)
        mt_null_ct = get_nulls(mt_numt_ct, mt_null, tmp)
        mt_numt_ct = mt_numt_ct.annotate_entries(null_hi = mt_null_ct[mt_numt.row_key, mt_numt.col_key].null_vec_hi,
                                                null_lo = mt_null_ct[mt_numt.row_key, mt_numt.col_key].null_vec_lo)
        mt_numt_ct = mt_numt_ct.select_entries(p_val_hi = mt_numt_ct.null_hi.filter(lambda x: x > mt_numt_ct.prop_hi).length()/mt_numt_ct.null_hi.length(),
                                            p_val_lo = mt_numt_ct.null_lo.filter(lambda x: x > mt_numt_ct.prop_lo).length()/mt_numt_ct.null_lo.length())
        ht_out = mt_numt_ct.entries()
        ht_out.export(args.o)
    else:
        mt_numt_ct.entries().export(args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script combines individual mitochondria coverage files and outputs a hail table with coverage annotations"
    )
    parser.add_argument(
        "-i",
        "--input-tsv",
        help="Input file with coverage files to combine in tab-delimited format of participant_id, base_level_coverage_metrics, sample",
        required=True,
    )
    parser.add_argument(
        "-o", "--output-ht", help="Name of ht to write output", required=True
    )
    parser.add_argument(
        "--prune-interval", help="Removes the first and last x bases from each interval, given by this term.", type=int, default=500
    )
    parser.add_argument(
        "--numt-only", help="If enabled, outputs a NUMT only table of the statistic as a function of sample and target.", action='store_true'
    )

    args = parser.parse_args()
    main(args)
