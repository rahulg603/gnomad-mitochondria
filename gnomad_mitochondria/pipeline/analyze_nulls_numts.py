import argparse
import hail as hl

suffix_dict = {'':0.025, '_5e_neg5': 5e-5, '_1e_neg6':1e-6, '_1e_neg4':1e-4, '_bonf':3e-9}

def prune_by_intervals(mt, bases, dont_parse_name=False):
    ht_annot = mt.rows()

    if dont_parse_name:
        ht_locs = ht_annot.group_by(ht_annot.target
                         ).aggregate(min_pos = hl.agg.min(ht_annot.locus.position) + bases, 
                                     min_pos_500 = hl.agg.min(ht_annot.locus.position) + 500,
                                     max_pos = hl.agg.max(ht_annot.locus.position) - bases,
                                     max_pos_500 = hl.agg.max(ht_annot.locus.position) - 500)
        ht_annot = ht_annot.annotate(**ht_locs[ht_annot.target])

    else:
        ht_annot = ht_annot.annotate(split_target = ht_annot.target.split('\\|'))
        ht_annot = ht_annot.annotate(split_positions = hl.map(lambda x: x.replace('^.+_(?=[0-9]{1,20}_[0-9]{1,20}_500bp)','').replace('_500bp','').split('_')[0:2], ht_annot.split_target))
        ht_annot = ht_annot.annotate(expected_min_pos = hl.min(hl.map(lambda x: hl.int32(x[0]), ht_annot.split_positions)),
                                    expected_max_pos = hl.max(hl.map(lambda x: hl.int32(x[1]), ht_annot.split_positions)))
        ht_annot = ht_annot.annotate(max_pos = ht_annot.expected_max_pos + 500 - bases,
                                    min_pos = ht_annot.expected_min_pos - 500 + bases)

    if ht_annot.filter(ht_annot.min_pos > ht_annot.max_pos).count() > 0:
        raise ValueError('ERROR: there are intervals where the number of bases to drop is greater than the total interval size.')
    
    ht_annot = ht_annot.filter(((ht_annot.locus.position >= ht_annot.min_pos)) & \
                               (ht_annot.locus.position <= ht_annot.max_pos))
    mt = mt.semi_join_rows(ht_annot)
    return mt


def annotate_by_pois(mt):
    # Determines if a sample has coverage that is too high or low using a Poisson distribution
    mt = mt.annotate_entries(p_hi = hl.ppois(mt.coverage, mt.mean_coverage, lower_tail=False),
                             p_lo = hl.ppois(mt.coverage, mt.mean_coverage, lower_tail=True))
    mt = mt.annotate_entries(**{"is_hi" + k: mt.p_hi <= v for k,v in suffix_dict.items()},
                             **{"is_lo" + k: mt.p_lo <= v for k,v in suffix_dict.items()})

    return mt


def get_per_target_stats(mt):
    suffixes = suffix_dict.keys()
    mt = mt.group_rows_by(mt.target
          ).aggregate_rows(N = hl.agg.count()
          ).aggregate_entries(**{'N_hi' + suffix: hl.agg.count_where(mt['is_hi' + suffix]) for suffix in suffixes},
                              **{'N_lo' + suffix: hl.agg.count_where(mt['is_lo' + suffix]) for suffix in suffixes}).result()
    mt = mt.annotate_entries(**{'prop_hi' + suffix: mt['N_hi' + suffix]/mt.N for suffix in suffixes},
                             **{'prop_lo' + suffix: mt['N_lo' + suffix]/mt.N for suffix in suffixes})
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

    # tmp = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/tmp/"
    # qc = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/qc_result_sample.tsv"
    # bases = 500
    # nulls = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/test_null/annotated_file_coverage_null.mt"
    # numts = "gs://fc-secure-f0dd1b4e-d639-4a0c-8712-6d02e9d8981a/test_null/annotated_file_coverage_numts.mt"


    ht_qc = hl.import_table(qc, impute=True)
    ht_qc = ht_qc.select(s = ht_qc['entity:qc_result_sample_id'], mean_coverage=ht_qc.mean_coverage).key_by('s')
    
    mt_numt = prune_by_intervals(hl.read_matrix_table(numts), bases)
    mt_numt = mt_numt.annotate_cols(mean_coverage = ht_qc[mt_numt.col_key].mean_coverage)
    mt_numt = annotate_by_pois(mt_numt)
    mt_numt_ct = get_per_target_stats(mt_numt)

    if not numt_only:
        mt_null = prune_by_intervals(hl.read_matrix_table(nulls), bases, dont_parse_name=True)
        mt_null = mt_null.annotate_cols(mean_coverage = ht_qc[mt_null.col_key].mean_coverage)
        mt_null = annotate_by_pois(mt_null)
        mt_null_ct = get_nulls(mt_numt_ct, mt_null, tmp)
        mt_numt_ct = mt_numt_ct.annotate_entries(null_hi = mt_null_ct[mt_numt.row_key, mt_numt.col_key].null_vec_hi,
                                                null_lo = mt_null_ct[mt_numt.row_key, mt_numt.col_key].null_vec_lo)
        mt_numt_ct = mt_numt_ct.select_entries(p_val_hi = mt_numt_ct.null_hi.filter(lambda x: x > mt_numt_ct.prop_hi).length()/mt_numt_ct.null_hi.length(),
                                            p_val_lo = mt_numt_ct.null_lo.filter(lambda x: x > mt_numt_ct.prop_lo).length()/mt_numt_ct.null_lo.length())
        ht_out = mt_numt_ct.entries()
        ht_out.export(args.output_tsv)
    else:
        mt_numt_ct.entries().export(args.output_tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script combines individual mitochondria coverage files and outputs a hail table with coverage annotations"
    )
    parser.add_argument(
        "--nulls", help="MT with imported null distribution files.",
    )    
    parser.add_argument(
        "--numts", help="MT with imported NUMT coverage files.", required=True,
    )
    parser.add_argument(
        "-o", "--output-tsv", help="Name of tsv to write output", required=True
    )
    parser.add_argument(
        "--prune-interval", help="Removes the first and last x bases from each interval, given by this term.", type=int, default=500
    )
    parser.add_argument(
        "--qc", help="Path to QC data. Required to annotate mean coverage and comptue statistics.", required=True
    )
    parser.add_argument(
        "--tmp", help="Folder for temporary files.", required=True
    )
    parser.add_argument(
        "--numt-only", help="If enabled, outputs a NUMT only table of the statistic as a function of sample and target.", action='store_true'
    )

    args = parser.parse_args()
    main(args)
