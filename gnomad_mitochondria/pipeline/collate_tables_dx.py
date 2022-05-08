import hail as hl
import pandas as pd
import argparse
import dxpy
import psutil, multiprocessing
import json
import sys, os, glob
import collections

from dxpy.utils import file_load_utils
from dxpy.bindings.download_all_inputs import _parallel_file_download, _get_num_parallel_threads, _create_dirs, \
                                              _sequential_file_download
from functools import reduce

FUSE_PREFIX = '/mnt/project/'


def make_input_json(d, fl_out="sample.json"):
    f = open(fl_out, "w")
    json.dump(d, f)
    f.close()
    return fl_out


def make_keyed_df(key, lst, key_to_suffix):
    df = pd.DataFrame({key: lst})
    df['s'] = df[key].apply(os.path.basename)
    df['s'] = df['s'].str.replace(key_to_suffix[key]+'$',"")
    return df


def reader1(args):
    idx, f = args
    df = pd.read_csv(f, index_col=0, header=None, sep='\t').transpose().assign(newcolthis=idx)
    return df


def reader2(f):
    return pd.read_csv(f, index_col=None, header=0, sep='\t')


def custom_get_job_input_filenames(input_dict):
    """Extract list of files, returns a set of directories to create, and
    a set of files, with sources and destinations. The paths created are
    relative to the input directory.

    Note: we go through file names inside arrays, and create a
    separate subdirectory for each. This avoids clobbering files when
    duplicate filenames appear in an array.
    """

    files = collections.defaultdict(list)  # dictionary, with empty lists as default elements
    dirs = []  # directories to create under <idir>

    # Local function for adding a file to the list of files to be created
    # for example:
    #    iname == "seq1"
    #    subdir == "015"
    #    value == { "$dnanexus_link": {
    #       "project": "project-BKJfY1j0b06Z4y8PX8bQ094f",
    #       "id": "file-BKQGkgQ0b06xG5560GGQ001B"
    #    }
    # will create a record describing that the file should
    # be downloaded into seq1/015/<filename>
    def add_file(iname, subdir, value, descr):
        if not dxpy.is_dxlink(value):
            return
        if descr is None:
            handler = dxpy.get_handler(value)
            if not isinstance(handler, dxpy.DXFile):
                return
            thisnm = handler.name
            thisid = handler.id
        else:
            handler = descr
            thisnm = handler['name']
            thisid = handler['id']
        filename = file_load_utils.make_unix_filename(thisnm)
        print('ding2')
        trg_dir = iname
        if subdir is not None:
            trg_dir = os.path.join(trg_dir, subdir)
        files[iname].append({'trg_fname': os.path.join(trg_dir, filename),
                             'handler': handler,
                             'src_file_id': thisid})
        dirs.append(trg_dir)

    # An array of inputs, for a single key. A directory
    # will be created per array entry. For example, if the input key is
    # FOO, and the inputs are {A, B, C}.vcf then, the directory structure
    # will be:
    #   <idir>/FOO/00/A.vcf
    #   <idir>/FOO/01/B.vcf
    #   <idir>/FOO/02/C.vcf
    def add_file_array(input_name, link_tuples):
        link_tuples = list(link_tuples)
        num_files = len(link_tuples)
        if num_files == 0:
            return
        num_digits = len(str(num_files - 1))
        dirs.append(input_name)
        for i, (link, descr) in enumerate(link_tuples):
            subdir = str(i).zfill(num_digits)
            add_file(input_name, subdir, link, descr)

    for input_name, value in list(input_dict.items()):
        if isinstance(value, zip):
            # This is a file array
            add_file_array(input_name, value)
        else:
            add_file(input_name, None, value, None)

    ## create a dictionary of the all non-file elements
    rest_hash = {key: val for key, val in list(input_dict.items()) if key not in files}
    return dirs, files, rest_hash


def custom_download_all_inputs(input_links, suffix_key, exclude=None, parallel=False, max_threads=8, enforce_nonmissing=True):
    '''
    :param exclude: List of input variables that should not be downloaded.
    :type exclude: Array of strings
    :param parallel: Should we download multiple files in parallel? (default: False)
    :type filename: boolean
    :param max_threads: If parallel is True, how many threads should be used
        to download files? (default: 8)
    :type append: int
    :returns: dict of lists of strings where each key is the input variable
                and each list element is the full path to the file that has
                been downloaded.

    A custom version of the download_all_inputs function that does not assume an input JSON file.
    By convention, if an input parameter "FOO" has value

        {"$dnanexus_link": "file-xxxx"}

    and filename INPUT.TXT, then the linked file will be downloaded into the
    path:

        $HOME/in/FOO/INPUT.TXT

    If an input is an array of files, then all files will be placed into
    numbered subdirectories under a parent directory named for the
    input. For example, if the input key is FOO, and the inputs are {A, B,
    C}.vcf then, the directory structure will be:

        $HOME/in/FOO/0/A.vcf
                     1/B.vcf
                     2/C.vcf

    Zero padding is used to ensure argument order. For example, if there are
    12 input files {A, B, C, D, E, F, G, H, I, J, K, L}.txt, the directory
    structure will be:

        $HOME/in/FOO/00/A.vcf
                     ...
                     11/L.vcf
    '''

    # Input directory, where all inputs are downloaded
    idir = file_load_utils.get_input_dir()
    try:
        dirs, inputs, rest = custom_get_job_input_filenames(input_links)
    except IOError:
        msg = 'Error: Could not find the input json file: {0}.\n'.format('<<DICT INPUT>>')
        msg += '       This function should only be called from within a running job.'
        print(msg)
        raise

    print('Filenames obtained.')
    # Exclude directories
    # dirs contain all folders (e.g. $HOME/in/FOO) and their sub folders (e.g. $HOME/in/FOO/1, $HOME/in/FOO/2, etc.)
    # If the main folder is excluded, its sub-folder would also be excluded from dirs_to_create
    dirs_to_create = []
    for d in dirs:
        keep = True
        if (exclude is not None) and (d is not None):
            if (d.split('/')[0] in exclude):
                keep = False
        if keep:
            dirs_to_create.append(d)

    # Create the directory structure, in preparation for download.
    # Allows performing the download in parallel.
    _create_dirs(idir, dirs_to_create)
    print('Directories created.')

    # Remove excluded inputs
    if exclude:
        inputs = file_load_utils.filter_dict(inputs, exclude)

    # Convert to a flat list of elements to download
    to_download = []
    for ival_list in inputs.values():
        to_download.extend(ival_list)

    # Download the files
    if parallel:
        total_mem = psutil.virtual_memory().total >> 20  # Total RAM in MB
        num_cores = multiprocessing.cpu_count()
        max_num_parallel_downloads = _get_num_parallel_threads(max_threads, num_cores, total_mem)
        sys.stderr.write("Downloading files using {} threads".format(max_num_parallel_downloads))
        _parallel_file_download(to_download, idir, max_num_parallel_downloads)
    else:
        _sequential_file_download(to_download, idir)

    # output a pandas table with columns for all downloaded samples with sample IDs
    trimmed_inputs = [make_keyed_df(k, [idir + '/' + subdict['trg_fname'] for subdict in v], suffix_key) for k, v in inputs.items()]
    joint_table = reduce(lambda x, y: pd.merge(x, y, on = 's', how='outer'), trimmed_inputs)
    if enforce_nonmissing:
        # verify that all tables are non-missing
        tf_null = joint_table.isnull().any().any()
        if tf_null:
            raise ValueError('ERROR: all df values should be non-missing.')
    
    return joint_table


def import_and_cat_tables(directory_df, id_col, path_col, new_id_col, append_ids_and_t=False, max_threads=8, filter_by=None, enforce_nonmiss=False):
    ids = [x for _, x in directory_df[id_col].iteritems()]
    directories = [x for _, x in directory_df[path_col].iteritems()]
    p = multiprocessing.Pool(processes=max_threads)

    to_subset_to = filter_by if filter_by is not None else ids
    empty_df = pd.DataFrame({new_id_col:[]})

    if append_ids_and_t:
        df_from_each_file = p.map(reader1, [(idx, f) for idx, f in zip(ids, directories) if idx in to_subset_to])
        if len(df_from_each_file) == 0:
            concatenated_df = empty_df
        else:
            concatenated_df = pd.concat(df_from_each_file, ignore_index=True, axis=0)
            concatenated_df = concatenated_df.rename({'newcolthis': new_id_col}, axis=1)
    else:
        df_from_each_file = p.map(reader2, [f for idx, f in zip(ids, directories) if idx in to_subset_to])
        if len(df_from_each_file) == 0:
            concatenated_df = empty_df
        else:
            concatenated_df = pd.concat(df_from_each_file, ignore_index=True, axis=0)

    # ensure that all ids are found in the new df
    if enforce_nonmiss:
        new_ids = [x for _, x in concatenated_df[new_id_col].iteritems()]
        tf_found_new = all([x_old in new_ids for x_old in ids])
        tf_found_old = all([x_new in ids for x_new in new_ids])
        if not (tf_found_new and tf_found_old):
            raise ValueError('ERROR: the list of individuals in directory_df must be the same as that obtained from the read tables.')
    
    return concatenated_df.reset_index(drop=True)


def fuse_find_data_objects(folder, suffix, recursive=True):
    path_search = FUSE_PREFIX + folder + '**/*' + suffix
    identified_objects = glob.glob(path_search, recursive=recursive)
    return identified_objects


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


def run_describe(lst):
    """
    dxpy.describe() only allows up to 1000 elements; here we loop across the size of list
    """
    n = len(lst)
    n_split = 999
    split_lists = [lst[i * n_split:(i + 1) * n_split] for i in range((n + n_split - 1) // n_split )] 
    return [y for x in split_lists for y in dxpy.describe(x)]


def main(pipeline_output_folder, vcf_suffix, coverage_suffix, mtstats_suffix, qc_stats_folder, qc_suffix,
         vcf_merging_output, coverage_calling_output):
    # download mito pipeline data
    data_dict = {'vcf': fuse_find_data_objects(pipeline_output_folder, vcf_suffix, False), 
                 'coverage': fuse_find_data_objects(pipeline_output_folder, coverage_suffix, False),
                 'stats': fuse_find_data_objects(pipeline_output_folder, mtstats_suffix, False)}
    downloaded_files = produce_fuse_file_table(data_dict, {'vcf':vcf_suffix, 'coverage':coverage_suffix, 'stats':mtstats_suffix})

    # read qc metrics
    data_dict_qc = {'qc': fuse_find_data_objects(qc_stats_folder, qc_suffix, False)}
    downloaded_qc_files = produce_fuse_file_table(data_dict_qc, {'qc':qc_suffix})

    # import stats and qc and merge all into table
    stats_table = import_and_cat_tables(downloaded_files, 's', 'stats', 's', enforce_nonmiss=True)
    qc_table = import_and_cat_tables(downloaded_qc_files, 's', 'qc', 's', append_ids_and_t=True, filter_by=list(stats_table['s']))
    final_analysis_table = stats_table.merge(qc_table, how='outer', on='s').merge(downloaded_files, how='inner', on='s')

    # produce tables for coverage and VCF merging
    final_analysis_table["s_2"] = final_analysis_table["s"]
    table_cov = final_analysis_table[["s", "coverage", "s_2"]].rename({'s':'collaborator_participant_id', 'coverage':'base_coverage'}).rename({'s_2':'s'})
    table_cov['coverage'] = 'file:///' + table_cov['coverage']
    table_vcf = final_analysis_table.rename({'s_2':'entity:participant_id'})
    table_vcf['vcf'] = 'file:///' + table_vcf['vcf']

    # output
    table_cov.to_csv(coverage_calling_output, sep='\t', index=False)
    table_vcf.to_csv(vcf_merging_output, sep='\t', index=False)


parser = argparse.ArgumentParser()
parser.add_argument('--pipeline-output-folder', type=str, required=True, 
                    help="Folder containing folders, each of which should contain pipeline outputs. Do not include project name. Should contain VCF, coverage, and diagnostic files.")
parser.add_argument('--vcf-merging-output', type=str, required=True, 
                    help="Local path to tsv to output which will contain all fields and will be usable for VCF merging.")
parser.add_argument('--coverage-calling-output', type=str, required=True, 
                    help="Local path to tsv to output which will contain 3 fields pointing to mtDNA per-base coverage for merging.")

parser.add_argument('--vcf-suffix', type=str, default='.self.ref.final.split.selfToRef.final.vcf',
                    help="Suffix of each final VCF to import.")
parser.add_argument('--coverage-suffix', type=str, default='_per_base_coverage.appended.liftedOver.tsv',
                    help="Suffix of each coverage tsv to import.")
parser.add_argument('--mtstats-suffix', type=str, default='_mtanalysis_diagnostic_statistics.tsv',
                    help="Suffix of each mtPipeline statistics file to import.")
parser.add_argument('--qc-stats-folder', type=str, default='/Bulk/Whole genome sequences/Concatenated QC Metrics/',
                    help="Folder containing folders, each of which should contain QC data from WGS. Also supports a single folder with relevant files in it. Do not include project name.")
parser.add_argument('--qc-suffix', type=str, default='.qaqc_metrics',
                    help="Suffix of each WGS QC file to import. Assumes that this file contains multiple rows for a single sample.")

pipeline_output_folder = '220403_mitopipeline_v2_2_ukb_trial/'
vcf_suffix = '.self.ref.final.split.selfToRef.final.vcf'
coverage_suffix = '_per_base_coverage.appended.liftedOver.tsv'
mtstats_suffix = '_mtanalysis_diagnostic_statistics.tsv'
qc_stats_folder = '/Bulk/Whole genome sequences/Concatenated QC Metrics/'
qc_suffix = '.qaqc_metrics'
vcf_merging_output = 'tab_vcf_merging.tsv'
coverage_calling_output = 'tab_coverage.tsv'

if __name__ == '__main__':
    args = parser.parse_args()
    main(**vars(args))