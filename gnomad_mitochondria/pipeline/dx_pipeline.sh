#!/bin/bash

pip3 install gnomad
git clone https://github.com/rahulg603/gnomad-mitochondria.git
git clone https://github.com/broadinstitute/gnomad_qc.git
git clone https://github.com/broadinstitute/gnomad_methods.git
mv gnomad_qc gnomad_qc_hold
cd gnomad_qc_hold
mv gnomad_qc ../
cd ..
mv gnomad-mitochondria gnomad_mitochondria_hold
mv gnomad_mitochondria_hold/gnomad_mitochondria ./
cd gnomad_mitochondria/pipeline/

python collate_tables_dx.py \
--pipeline-output-folder '220403_mitopipeline_v2_2_ukb_trial' \
--vcf-merging-output 'tab_vcf_merging.tsv' \
--coverage-calling-output 'tab_coverage.tsv'

mkdir pipelineworkingdir
mkdir pipelineworkingdir/tmp
mkdir pipelineworkingdir/coverage
mkdir pipelineworkingdir/vcf

python annotate_coverage.py \
-i 'tab_coverage.tsv' \
-o pipelineworkingdir/coverage/testing_pipeline_trial_ukb_coverage.ht \
-t pipelineworkingdir/tmp/ \
--overwrite

python combine_vcfs.py \
-p 'tab_vcf_merging.tsv' \
-c pipelineworkingdir/coverage/testing_pipeline_trial_ukb_coverage.mt \
-v vcf \
-a '/mnt/project/reference/grch38_genome/blacklist_sites.hg38.chrM.bed' \
-a-ref GRCh38 \
-o pipelineworkingdir/vcf/ \
-t /pipelineworkingdir/tmp/ \
-f 220506_testing_ukb_pipeline_variants \
--overwrite --include-extra-v2-fields