#!/bin/bash

#### PARAMETERS
export JOBLIMIT=5000 # how many jobs to allow simulataneously
export PORTID=8094
export USE_MEM=80
export SQL_DB_NAME="local_cromwell_run.db" # name of local SQL database
export outputFold=221105_trial
export n_max=2

#### INSTALL DEPENDENCIES
pip install pyhocon
git clone https://github.com/rahulg603/gnomad-mitochondria.git


#### DOWNLOAD DATA
curl https://personal.broadinstitute.org/rahul/validate_vcf/validate_vcf.wdl -o validate_vcf.wdl
curl https://github.com/broadinstitute/cromwell/releases/download/77/cromwell-77.jar -o cromwell-77.jar -L
curl https://github.com/broadinstitute/cromwell/releases/download/77/womtool-77.jar -o womtool-77.jar -L
curl -s "https://get.sdkman.io" -o install_sdkman.sh
bash install_sdkman.sh


#### CONFIGURE
gcloud auth list --format json | jq -r .[0].account


#### GET FILEs
gsutil -u $GOOGLE_PROJECT ls gs://fc-aou-datasets-controlled/v5/microarray/vcf/single_sample/*.vcf.gz | head -n${n_max} | sort --version-sort  > vcf_file_list.txt
gsutil -u $GOOGLE_PROJECT ls gs://fc-aou-datasets-controlled/v5/microarray/vcf/single_sample/*.vcf.gz.tbi | head -n${n_max} | sort --version-sort  > vcf_tbi_file_list.txt


#### MUNGE INPUTS AND PREPARE COMMANDS
python <<'CODE'
import os, re
import json
import pandas as pd
from copy import deepcopy
from pyhocon import ConfigFactory, ConfigTree, HOCONConverter
from google.cloud import storage

# get globals
mem = int(os.getenv('USE_MEM'))
joblimit = int(os.getenv('JOBLIMIT'))
port = int(os.getenv('PORTID'))
sql_db = os.getenv("SQL_DB_NAME")
bucket = os.getenv("WORKSPACE_BUCKET")
project = os.getenv("GOOGLE_PROJECT")
path_indiv_save = os.getenv("outputFold") + '/'
n_max = int(os.getenv("n_max"))


# set up cromwell directories
cromwell_test_workdir = bucket + "/"  # Later, "cromwell-executions" will be appended to this for cromwell-workflow storage.
output_bucket = bucket + "/" + os.getenv("outputFold")  # This is where the output of the WDL will be.
 
print(f'Workspace bucket: {bucket}')
print(f'Workspace project: {project}')
print(f'Workspace cromwell working bucket: {cromwell_test_workdir}')
print(f'Workspace output bucket: {output_bucket}')

# set up filenames
options_filename = "options.json"
wdl_filename = "validate_vcf.wdl"

# create options file
# options_content = f'{{\n  "jes_gcs_root": "{output_bucket}",\n  "workflow_failure_mode": "NoNewCalls"\n}}\n'
options_content = f'{{\n  "jes_gcs_root": "{output_bucket}"\n}}\n'

fp = open(options_filename, 'w')
fp.write(options_content)
fp.close()
print(options_content)

# create cromwell configuration, lifting the job limit and exporting to local SQL database
with open('/home/jupyter/cromwell.conf', 'r') as f:
    input_conf = f.read()
include_str_repl = 'include required(classpath("application"))\n\n'
input_conf_rm = input_conf.replace(include_str_repl, '')
cromwell_config_file = ConfigFactory.parse_string(input_conf_rm)
cromwell_config_file['system'] = ConfigTree({'new-workflow-poll-rate': 1,
                                             'max-concurrent-workflows': joblimit,
                                             'max-workflow-launch-count': 400,
                                             'job-rate-control': ConfigTree({'jobs': 50,
                                                                             'per': '3 seconds'})})
cromwell_config_file['backend']['providers']['PAPIv2-beta']['config']['concurrent-job-limit'] = joblimit
cromwell_config_file['backend']['providers']['PAPIv2-beta']['config']['genomics']['enable-fuse'] = True
cromwell_config_file['database'] = ConfigTree({'profile': "slick.jdbc.HsqldbProfile$",
                                               'insert-batch-size': 6000,
                                               'db': ConfigTree({'driver':"org.hsqldb.jdbcDriver", 
                                                                 'url':f'jdbc:hsqldb:file:{sql_db};shutdown=false;hsqldb.default_table_type=cached;hsqldb.tx=mvcc;hsqldb.large_data=true;hsqldb.lob_compressed=true;hsqldb.script_format=3;hsqldb.result_max_memory_rows=20000',
                                                                 'connectionTimeout': 300000})})
with open('/home/jupyter/cromwell.new.conf', 'w') as f:
    f.write(include_str_repl + HOCONConverter.to_hocon(cromwell_config_file))

cromwell_run_cmd = f'source "/home/jupyter/.sdkman/bin/sdkman-init.sh" &&  sdk install java 11.0.14-tem && echo "Validating WDL..." && java -jar womtool-77.jar validate _WDL_FILE_ && java -Xmx{str(mem)}g -classpath ".:sqlite-jdbc.jar" -Dconfig.file=/home/jupyter/cromwell.new.conf -Dwebservice.port={str(port)} -jar cromwell-77.jar server'
cromwell_run_cmd_final = cromwell_run_cmd.replace("_WDL_FILE_", wdl_filename)

with open(f"cromwell_startup_script.sh", "w") as text_file:
    text_file.write("#!/bin/bash\n")
    text_file.write(cromwell_run_cmd_final + '\n')

json_collection = []
this_json = {
  "validate_vcf.ordered_vcf_shards_idx_list": "vcf_tbi_file_list.txt",
  "validate_vcf.ref_fasta": "gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "validate_vcf.ref_fai": "gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "validate_vcf.ref_dict": "gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict",
  "validate_vcf.ordered_vcf_shards_list": "vcf_file_list.txt"
}
json_collection.append(this_json)

batch_json_filename = f"batch_input_allofus.json"
with open(batch_json_filename, 'w') as f:
    json.dump(json_collection, f)
    
with open('ct_submissions.txt', 'w') as f:
    f.write(str(len(json_collection)))

batch_cromwell_cmd = f'curl -X POST "http://localhost:{str(port)}/api/workflows/v1/batch" -H "accept: application/json" -F workflowSource=@_WDL_FILE_ -F workflowInputs=@_INPUTS_ -F workflowOptions=@_OPTIONS_FILE_'
batch_cromwell_cmd = batch_cromwell_cmd.replace("_WDL_FILE_", wdl_filename)
batch_cromwell_cmd = batch_cromwell_cmd.replace("_INPUTS_", batch_json_filename)
batch_cromwell_cmd = batch_cromwell_cmd.replace("_OPTIONS_FILE_", options_filename)

count_n_submit = "submission_count=$(grep -o 'Submitted' batch_submission_ids.txt | wc -l)\n"
test_ct = 'if [ "$submission_count" -ne "$(cat ct_submissions.txt)" ]; then echo "ERROR: submission count is incorrect."; exit 1; fi\n'
get_batch_ids = 'cat batch_submission_ids.txt | sed \'s/{"id"://g\' | sed \'s/","status":"Submitted"}//g\' | sed \'s/"//g\' | sed \'s/,/\\n/g\' | sed \'s/\\[//g\' | sed \'s/\\]//g\' > ordered_batch_ids.txt\n'

with open(f"cromwell_submission_script_batch.sh", "w") as text_file:
    text_file.write("#!/bin/bash\n")
    text_file.write(batch_cromwell_cmd + ' | tee batch_submission_ids.txt\n')
    text_file.write(count_n_submit)
    text_file.write(test_ct)
    text_file.write(get_batch_ids)
    text_file.write('echo "" >> ordered_batch_ids.txt\n')
CODE


#### LAUNCH SERVER AS A SUBPROCESS
chmod 777 cromwell_startup_script.sh
setsid ./cromwell_startup_script.sh > cromwell_server_stdout.log 2>cromwell_server_stderr.log &


#### CREATE MONITORING COMMAND
sleep 120
echo "Server started."
echo "Here is the tail of the current stdout.log. Examine this to make sure the server is running:"
tail -n10 cromwell_server_stdout.log
echo ""
echo "Run cromwell_submission_script_batch.sh to submit the desired jobs."
echo ""