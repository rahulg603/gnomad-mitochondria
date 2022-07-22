import pandas as pd
import argparse
import os, re
from datetime import datetime
from google.cloud import storage


def check_table(df, subj_col):
    tf_dupe = any(df.duplicated(subset=[subj_col]))
    if tf_dupe:
        raise ValueError('ERROR: dataframe cannot have duplicates.')


parser = argparse.ArgumentParser()
parser.add_argument('--database-stats', type=str, default="mt_pipeline_single_2_5_stats.tsv",
                    help="Name of file containing all completed sample-level statistics. If doesn't exist, will create newly.")
parser.add_argument('--database-paths', type=str, default='mt_pipeline_single_2_5_paths.tsv',
                    help="Name of file containing paths to all completed VCF and coverage files. If doesn't exist, will create newly.")
parser.add_argument('--new-stats', type=str, required=True,
                    help="Name of file with sample-level stats to add to the database.")
parser.add_argument('--new-paths', type=str, required=True,
                    help="Name of file with paths to completed VCF and coverage file to add to the database.")

if __name__ == "__main__":
    args = parser.parse_args()
    bucket_name = re.sub('^gs://', '', os.getenv("WORKSPACE_BUCKET"))
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    new_suff = '_' + datetime.now(tz=None).strftime('%d_%b_%y_%H.%M.%S')
    
    df_stats = pd.read_csv(f"gs://{bucket_name}/{args.new_stats}", sep='\t')
    check_table(df_stats, 's')
    if storage.Blob(bucket=bucket, name=f"gs://{bucket_name}/{args.database_stats}").exists(storage_client):
        df_stats_exist = pd.read_csv(f"gs://{bucket_name}/{args.database_stats}", sep='\t')
        this_path_spl = os.path.splitext(f"gs://{bucket_name}/{args.database_stats}")
        df_stats_exist.to_csv(this_path_spl[0] + new_suff + this_path_spl[1], sep='\t', index=False)
        
        df_out = pd.concat([df_stats_exist, df_stats], axis=0)
        check_table(df_out, 's')
        df_out.to_csv(f"gs://{bucket_name}/{args.database_stats}", sep='\t', index=False)
    else:
        df_stats.to_csv(f"gs://{bucket_name}/{args.database_stats}", sep='\t', index=False)

    df_paths = pd.read_csv(f"gs://{bucket_name}/{args.new_paths}", sep='\t')
    check_table(df_paths, 'batch')
    if storage.Blob(bucket=bucket, name=args.database_paths).exists(storage_client):
        df_paths_exist = pd.read_csv(f"gs://{bucket_name}/{args.database_paths}", sep='\t')
        this_path_spl = os.path.splitext(f"gs://{bucket_name}/{args.database_paths}")
        df_paths_exist.to_csv(this_path_spl[0] + new_suff + this_path_spl[1], sep='\t', index=False)
        
        df_out_paths = pd.concat([df_paths_exist, df_paths], axis=0)
        check_table(df_out_paths, 'batch')
        df_out_paths.to_csv(f"gs://{bucket_name}/{args.database_paths}", sep='\t', index=False)
    else:
        df_paths.to_csv(f"gs://{bucket_name}/{args.database_paths}", sep='\t', index=False)