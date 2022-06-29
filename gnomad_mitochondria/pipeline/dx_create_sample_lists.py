import pandas as pd
import argparse
import dxpy
import re
import os
import multiprocessing
import glob
from functools import reduce

FUSE_PREFIX = '/mnt/project/'

def fuse_find_data_objects(folder, suffix, recursive=True):
    path_search = FUSE_PREFIX + folder + '**/*' + suffix
    identified_objects = glob.glob(path_search, recursive=recursive)
    return identified_objects


def make_keyed_df(key, lst, key_to_suffix):
    df = pd.DataFrame({key: lst})
    df['s'] = df[key].apply(os.path.basename)
    df['s'] = df['s'].str.replace(key_to_suffix[key]+'$',"")
    return df


def produce_fuse_file_table(input_links, suffix_key, enforce_nonmissing=True):
    # output a pandas table with columns for all downloaded samples with sample IDs
    trimmed_inputs = [make_keyed_df(k, v, suffix_key) for k, v in input_links.items()]
    joint_table = reduce(lambda x, y: pd.merge(x, y, on = 's', how='outer'), trimmed_inputs)
    if enforce_nonmissing:
        # verify that all tables are non-missing
        tf_null = joint_table.isnull().any().any()
        if tf_null:
            raise ValueError('ERROR: all df values should be non-missing.')
    
    return joint_table


def reader1(args):
    idx, f = args
    df = pd.read_csv(f, index_col=0, header=None, sep='\t').transpose().assign(newcolthis=idx)
    return df


def reader2(f):
    return pd.read_csv(f, index_col=None, header=0, sep='\t')


def import_and_cat_tables(directory_df, batch_col, path_col, new_batch_col, append_ids_and_t=False, max_threads=8, filter_by=None, enforce_nonmiss=False):
    batches = [x for _, x in directory_df[batch_col].iteritems()]
    directories = [x for _, x in directory_df[path_col].iteritems()]
    p = multiprocessing.Pool(processes=max_threads)

    to_subset_to = filter_by if filter_by is not None else batches
    empty_df = pd.DataFrame({new_batch_col:[]})

    if append_ids_and_t:
        df_from_each_file = p.map(reader1, [(idx, f) for idx, f in zip(batches, directories) if idx in to_subset_to])
        if len(df_from_each_file) == 0:
            concatenated_df = empty_df
        else:
            concatenated_df = pd.concat(df_from_each_file, ignore_index=True, axis=0)
            concatenated_df = concatenated_df.rename({'newcolthis': new_batch_col}, axis=1)
    else:
        df_from_each_file = p.map(reader2, [f for idx, f in zip(batches, directories) if idx in to_subset_to])
        if len(df_from_each_file) == 0:
            concatenated_df = empty_df
        else:
            concatenated_df = pd.concat(df_from_each_file, ignore_index=True, axis=0)

    # ensure that all ids are found in the new df
    if enforce_nonmiss:
        new_ids = [x for _, x in concatenated_df[new_batch_col].iteritems()]
        tf_found_new = all([x_old in new_ids for x_old in batches])
        tf_found_old = all([x_new in batches for x_new in new_ids])
        if not (tf_found_new and tf_found_old):
            raise ValueError('ERROR: the list of individuals in directory_df must be the same as that obtained from the read tables.')
    
    return concatenated_df.reset_index(drop=True)



def obj_to_path(obj):
    this_path = f"file:///{FUSE_PREFIX}{obj['describe']['folder']}/{obj['describe']['name']}"
    return this_path


def get_batch_name(path):
    return os.path.basename(os.path.dirname(path))


def get_file_type(path):
    thisbase = os.path.basename(path)
    if thisbase == 'batch_analysis_statistics.tsv':
        thisclass = "stats"
    elif thisbase == 'batch_merged_mt_calls.vcf.bgz':
        thisclass = "variants"
    elif thisbase == 'batch_merged_mt_coverage.tsv.bgz':
        thisclass = "coverage"
    elif thisbase == 'batch_idxstats_metrics.tsv.gz':
        thisclass = "idxstats"
    else:
        raise ValueError('ERROR: files do not fit spec.')
    return thisclass


parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, help='DX path to batch folders.', required=True)
parser.add_argument('--subfolder-match', type=str, help='Only include files in subflolders matching this regex.')

if __name__ == '__main__':
    args = parser.parse_args()
    data_obj = dxpy.find_data_objects(classname='file', name='batch*', name_mode='glob', folder=args.path, describe=True)
    if args.subfolder_match is not None:
        data_obj = [x for x in data_obj if re.search(f'{args.path}/{args.subfolder_match}', x['describe']['folder'])]

    file_paths = [obj_to_path(x) for x in data_obj]
    batch_ids = [get_batch_name(x) for x in file_paths]
    types = [get_file_type(x) for x in file_paths]

    df = pd.DataFrame({'path': file_paths, 'type': types, 'batch_id': batch_ids})

    # combine stats tables here
    stats_table = import_and_cat_tables(df[df.type == 'stats',:], 'batch_id', 'path', 'batch_id', enforce_nonmiss=True)

    # read qc metrics
    data_dict_qc = {'qc': fuse_find_data_objects(qc_stats_folder, qc_suffix, False)}
    downloaded_qc_files = produce_fuse_file_table(data_dict_qc, {'qc':qc_suffix})
    qc_table = import_and_cat_tables(downloaded_qc_files, 's', 'qc', 's', append_ids_and_t=True, filter_by=list(stats_table['s']))
    final_info_table = stats_table.merge(qc_table, how='outer', on='s')

    # output files
    table_cov.to_csv(coverage_calling_output, sep='\t', index=False)
    table_vcf.to_csv(vcf_merging_output, sep='\t', index=False)

    # output to sql
    ht_cov = hl.Table.from_pandas(table_cov)
    ht_cov.export(f'dnax://{my_database}/{coverage_calling_output}')
    ht_vcf = hl.Table.from_pandas(table_vcf)
    ht_vcf.export(f'dnax://{my_database}/{vcf_merging_output}')