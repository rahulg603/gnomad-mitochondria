# mitochondria
Scripts to analyze mitochondria data. These scripts combine and process the output from running the mitochondria mode of Mutect2 through [Terra](https://terra.bio/) on several samples. The most recent version of the WDL is available [here](https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl).

## Step 1: annotate_coverage.py

The WDL outputs one file per sample containing per base coverage across the entire mitochondria calculated from Picard's CollectWgsMetrics. Running annotate_coverage.py will combine the per base coverage files across many samples into a mt(MatrixTable), ht(HailTable), and tsv file, and will also calculate the following aggregate statistics per base:
* Mean coverage
* Median coverage
* Fraction of samples with > 100x coverage
* Fraction of samples with > 10000x coverage

Required inputs:
* input-tsv: Input file with coverage files to combine in tab-delimited format of participant_id, path to the per base coverage for that sample, sample name
* output-ht: Name of ht to write output
* temp-dir: Temporary directory to use for intermediate outputs


The mt that is output is used as input to combine_vcfs.py.

## Step 2: combine_vcfs.py

The combine_vcfs.py script takes individual sample vcfs and combines them into one vcf/mt. This script:

* Combines individual sample vcfs into one vcf
* Removes the "possible_numt", "mt_many_low_hets", and "FAIL" FT filters because these filters were found to have low performance but still may be present if running earlier versions of the Mutect2 mitochondria pipeline
* For sites without a call in a particular sample, sets the call to missing if the sample has coverage <= minimum_homref_coverage (default: 100x) at that position, otherwise sets the call to homoplasmic reference
* Applies the "artifact_prone_site" filter to any SNP or deletion that spans a known problematic site supplied in the "artifact_prone_sites_path" bed file

Required inputs:
* participant-data: Participant data (the downloaded data tab from Terra), should be a tab-delimited file with at minimum columns for 'entity:participant_id'(sample name with prohibited-characters replaced with underscores), 's'(sample name), and VCF output (path to the VCF output by mutect, name of this column is supplied to the vcf_col_name parameter)
* participants-to-subset: Path to txt file of participant_ids to which the data should be subset (file should contain header and one line for each participant_id matching the 'entity:participant_id's supplied in Terra
* coverage-mt-path: Path to MatrixTable of sample-level coverage (per-sample and per-base, can be generated by running annotate_coverage.py
* vcf-col-name: Name of column in participant data file that contains the path to the VCF output by Mutect2
* artifact-prone-sites-path: Path to BED file of artifact-prone sites to flag in the FILTER column
* output-bucket: Path to bucket to which results should be written
* temp-dir: Temporary directory to use for intermediate outputs
* file-name: File name to use for output files (will be used the the .vcf.bgz and .mt outputs)

## Step 3: add_annotations.py

Annotations to the mt/vcf are added by running the add_annotations.py script. Some of the main functions in this script include:
* Sets the GT to heteroplasmic (0/1) or homoplasmic (1/1) based on the supplied min-het-threshold (default: 0.10) and min-hom-threshold (0.95)
* Adds sample annotations from data download in Terra ("participant_id", "contamination", "freemix_percentage", "major_haplogroup", "wgs_median_coverage", "mt_mean_coverage"). "contamination", "major_haplogroup", and "mt_mean_coverage" are output by the Mutect2 mitochondria pipeline, whereas "freemix_percentage" (from VerifyBamID) and "wgs_median_coverage" (from CollectWgsMetrics) should be uploaded by the user (as well as "age" and "pop" if the "subset-to-gnomad-release" argument cannot be used). Column "s" (Sample ID) is used as the key for joining.
* Annotates a variant as "hap_defining_variant" (based on Phylotree), and adds in-silico predictors of variant severity (from PON-mt-tRNA and MitoTIP)
* Adds sample metadata, at minimum inferred nuclear ancestry ("pop") should be supplied
* Filters out samples with low (< 50) or high (> 500) estimated mitochondria copy number
* Filters out samples with contamination > 2% (based on nuclear contamination, mitochondria contamination, and internal algorithm)
* Adds VEP and dbSNP annotations
* Applies two variant-level filters ("indel_stack" and "ngp"). Definitions of the variant-level filters:
	* "artifact_prone_site": Variant overlaps a site that is known to be problematic. Applied in combine_vcfs.py.
	* "indel_stack": All samples with the given indel were multiallelic for heteroplasmic indels at the site.
	* "npg": No pass genotype. No sample had a pass genotype for the variant.
* Flags "common_low_heteroplasmy" variants, which are variants found at an overall frequency of > 0.001 across all samples with a PASS genotype and heteroplasmy level > 0 and < 50% (includes variants < vaf-filter-threshold  which are subsequently filtered)
* Sets genotypes < vaf-filter-threshold (default: 0.01) to homoplasmic reference
* Sets non-PASS genotypes and genotypes with heteroplasmy < min-het-threshold (default: 0.10) to missing
* Outputs summary of variant statistics and applied filters ("stats.txt" and "stats_pass.txt")

Required inputs:
* mt-path: Path to combined mt
* output-dir: Path to directory to which output should be written
* participant-data: Output file that results from Terra data download
* vep-results: MatrixTable path to output vep results (either the existing results or where to ouput new vep results if also setting run_vep)

Output will include an annotated mt.vcf, sites-only ht/vcf, and simplified txt file containing only the variant-identifying information as well as "filters", "AC_hom", "AC_het", "AF_hom", "AF_het", "AN", "max_hl".

NOTE: vaf-filter-threshold should match what was supplied as the vaf_filter_threshold when running the Mutect2 mitochondria pipeline.



An example of the VEP init file can be found in the resources directory. Example for spinning up a cluster to run VEP:
```
hailctl dataproc start CLUSTER_NAME \
--worker-machine-type=n1-highmem-4  \
--worker-boot-disk-size=200 \
--preemptible-worker-boot-disk-size=200 \
--metadata=VEP_CONFIG_PATH=/vep_data/vep-gcloud.json,VEP_CONFIG_URI=file:///vep_data/vep-gcloud.json,VEP_REPLICATE=us \
--init=PATH_TO_VEP_INIT \
--requester-pays-allow-all \
--project PROJECT_NAME \
--num-preemptible-workers 15 \
--max-idle 2h \
--properties=spark:spark.speculation=true
```

Scripts assume samples have been processed on GRCh38. More information on details of the pipeline and detailed descriptions of annotations can be found on our [blog](https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/).
