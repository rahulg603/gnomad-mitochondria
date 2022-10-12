# Pipeline to filter sample_stats file and add haplogroup / mtcn annotations

import hail as hl
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--sample-stats', required=True, help='Sample stats flat file.')
parser.add_argument('--contamination-drop', required=True, help='Samples that are contaminated, IDENTIFIED BY RUNNING add_annotation.py WITHOUT SAMPLE FILTERING.')
parser.add_argument('--out', required=True, help='gs:// path for filtered sample stats table')


if __name__ == '__main__':
    args = parser.parse_args()
    stat_ht = hl.import_table(args.sample_stats, impute=True, key='s')
    contamination_ht = hl.import_table(args.contamination_drop, impute=True, key='s')

    # filter out contaminated samples
    stat_ht = stat_ht.filter(contamination_ht[stat_ht.s].keep)

    # filter out samples with overlapping homoplasmies
    stat_ht = stat_ht.filter(stat_ht.num_mt_overlaps == 0)
    
    # add some annotations
    stat_ht = stat_ht.annotate(
        stat_ht=hl.if_else(
            stat_ht.major_haplogroup.startswith("HV")
            | stat_ht.major_haplogroup.startswith("L"),
            stat_ht.major_haplogroup[0:2],
            stat_ht.major_haplogroup[0],
        ),
        mtcn = 2 * stat_ht.mean_coverage / stat_ht.nuc_mean_coverage,
        mtcn_median = 2 * stat_ht.median_coverage / stat_ht.nuc_mean_coverage
    )

    # export
    stat_ht.export(args.out)