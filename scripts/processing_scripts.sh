#!/bin/bash

# This script serves for processing the raw sequencing data of the PSC study (Czech cohort). 
# The Norwegian cohort was processed with the same script but with different library names. 
# The merging of these two cohorts was done afterward with a custom R script.
# This script includes the 1. preprocessing, 2. denoising and taxonomic assignment, 3. merging and filtering the data, 4. exporting
# The data are processed in per-library basis and merged afterwards.

# Requirements:
# - bbmap
# - cutadapt
# - vsearch
# - qiime2
# - pre-built taxonomic classifier (see building_rescript_classifier.sh)

# conda environment (Optional)
conda activate /path/to/conda_environment

# set TMPDIR (Optional)
export TMPDIR=/path/to/tmp

# Set path for inputs and outputs
path=/path/to/raw/data
path_results=/path/to/results/directory
path_base=/path/to/base/directory
path_bbmap=/path/to/bbmap
path_db=/path/to/silva_dbs

##############################################################
# STAGE 1: PREPROCESSING
##############################################################

# iterate through each library
for lib_name in cz_lib2 cz_lib3 cz_lib4 cz_lib5 cz_lib6 cz_lib7 cz_lib8 
do

	# move to the library directory
	cd ${path}/${lib_name}/

	## Removing adapters and phiX reads using bbamp
	# make directory for the results
	mkdir ${path_results}/${lib_name}
	mkdir ${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk

	# iterate through every sammple
	for read1 in *R1_001.fastq*
	do	
		read2=$(echo $read1| sed 's/R1_/R2_/')

		# call bbduk to remove adapters and phix 
		${path_bbmap}bbduk.sh \
			ref=${path_bbmap}resources/adapters.fa \
			in1=${read1} in2=${read2} out=stdout.fq \
			k=23 hdist=1 tbo cf=TRUE ftm=5 \
			2> ${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk/remove_adaptors.log | \
		
		${path_bbmap}bbduk.sh \
			ref=${path_bbmap}resources/phix174_ill.ref.fa.gz \
			in=stdin.fq int=t \
			out1=${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk/noadaptors_nophix_${read1} \
			out2=${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk/noadaptors_nophix_${read2} \
			k=31 hdist=1 2> ${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk/remove_phiX.log
	done

	## Removing primers with cutadapt, they are saved in the CSV file primers.csv

	# make directory for the results
	mkdir ${path_results}/${lib_name}/trimmed
	
	# move to the results from the last step
	cd ${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk

	# iterate through every sample
	for read1 in *R1_001.fastq*
	do
		read2=$(echo $read1| sed 's/R1_/R2_/')
		sample1=$(echo $read1| sed 's/noadaptors_nophix_//')
		sample2=$(echo $read2| sed 's/noadaptors_nophix_//')
		command1="cutadapt "
		cd $path_base
		for item in $(cat primers.csv) 
		do
			FWD=$(echo $item | cut -f2 -d ",")
			REV=$(echo $item | cut -f5 -d ",")

			command1+="-g ^${FWD} -G ^${REV} "
		done
		command1+="--discard-untrimmed -j 5 -o ${path_results}/${lib_name}/trimmed/trimmed_${sample1} \
		-p ${path_results}/${lib_name}/trimmed/trimmed_${sample2} \
		${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk/${read1} \
		${path_results}/${lib_name}/remove_adapters_and_phiX_bbduk/${read2}"

		eval $command1
		
	done


	## Merging paired-end reads using bbmerge

	# make directory for the results
	mkdir ${path_results}/${lib_name}/merge_pairs_bbmerge
	
	# move to the results from the last step
	cd ${path_results}/${lib_name}/trimmed

	# iterate through every sample
	for read1 in *R1_001.fastq*
	do

		read2=$(echo $read1| sed 's/R1/R2/')
		sample1=$(echo $read1| sed 's/trimmed_//')
		sample2=$(echo $read2| sed 's/trimmed_//')
		sample=$(echo $read1| sed 's/_L001_R1_001.fastq.gz//')
		
		${path_bbmap}bbmerge.sh \
			in1=${path_results}/${lib_name}/trimmed/${read1} \
			in2=${path_results}/${lib_name}/trimmed/${read2} \
			out=stdout.fq qtrim=r trimq=15 maxlength=500 mininsert=350 ihist=../merge_pairs_bbmerge/${sample}_insert_histogram.ihist | \
		${path_bbmap}reformat.sh -Xmx20g -int=f -maxns=0 \
			in=stdin.fq \
			out=${path_results}/${lib_name}/merge_pairs_bbmerge/${sample}_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz \
		lhist=${path_results}/${lib_name}/merge_pairs_bbmerge/${sample}_histogram_after_N_filter.txt;
	done

	## Orienting the merged reads using VSEARCH

	# make directory for the results
	mkdir ${path_results}/${lib_name}/oriented

	# move to the results from the last step
	cd ${path_results}/${lib_name}/merge_pairs_bbmerge

	# iterate through every sample
	for read in *mergedpairs.fastq.gz
	do
		sample=$(echo $read| sed 's/trimmed_//')
		sample=$(echo $sample| sed 's/_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz//')
		vsearch --orient $read \
			--db ${path_db}/silva_nrr99_v138.1_orienting.udb \
			--fastqout ${path_results}/${lib_name}/oriented/${sample}_oriented.fq \
			--tabbedout ${path_results}/${lib_name}/oriented/orient.tx
	done



##############################################################
# STAGE 2: DENOISING and TAXONOMIC ASSIGNMENT
##############################################################

## conda environment (Optional)
conda activate /path/to/qiime_environment

for lib_name in cz_lib2 cz_lib3 cz_lib4 cz_lib5 cz_lib6 cz_lib7 cz_lib8 
do
	## IMPORT TO QIIME

	# move to the results from the last step
	cd ${path_results}/${lib_name}/oriented

	# make directory for the results
	mkdir ${path_results}/${lib_name}/qiime2
	mkdir ${path_results}/${lib_name}/output

	# create manifest
	echo -e "sample-id\tabsolute-filepath\tdirection" > ../qiime2/manifest_${lib_name}.tsv
	for read in *-oriented.fq; do
		sample=$(echo $read| sed 's/-oriented.fq//')
		echo -e "${sample}\t${path_results}${lib_name}/oriented/${read}\tforward" >>  ../qiime2/manifest_${lib_name}.tsv
	done

	# import to qiime
	cd ${path_results}/${lib_name}/qiime2

	qiime tools import \
    	--input-path manifest_${lib_name}.tsv \
    	--type 'SampleData[SequencesWithQuality]' \
    	--input-format SingleEndFastqManifestPhred33V2 \
    	--output-path ${lib_name}.qza
	
	# create summarized visualization
    	qiime demux summarize --i-data ${lib_name}.qza --o-visualization output/${lib_name}_imported_seqs.qzv --p-n 10000

	# use DEBLUR to denoise the reads
    	qiime deblur denoise-16S \
	--i-demultiplexed-seqs ${lib_name}.qza \
	--p-trim-length 400 \
	--p-sample-stats \
	--o-representative-sequences output/${lib_name}_ASV_seqs.qza \
	--o-table output/${lib_name}_ASV_abundance.qza \
	--o-stats output/${lib_name}_stats_deblur.qza
	
	# use feature-classifier classify-sklearn to assign the taxonomy
	qiime feature-classifier classify-sklearn \
	--i-reads ${path_results}/${lib_name}/qiime2/output/${lib_name}_ASV_seqs.qza \
	--i-classifier ${path_db}/silva-138.1-ssu-nr99-319f-806r-classifier.qza \
	--p-n-jobs -1 \
	--o-classification ${path_results}/${lib_name}/qiime2/output/${lib_name}_ASV_taxonomy.qza

done

##############################################################
# STAGE 3: MERGING LIBRARIES AND FILTERING
#############################################################

cd ${path_results}

qiime feature-table merge \
    --i-tables cz_lib2/qiime2/output/cz_lib2_ASV_abundance.qza \
    --i-tables cz_lib3/qiime2/output/cz_lib3_ASV_abundance.qza\
    --i-tables cz_lib4/qiime2/output/cz_lib4_ASV_abundance.qza\
    --i-tables cz_lib5/qiime2/output/cz_lib5_ASV_abundance.qza\
    --i-tables cz_lib6/qiime2/output/cz_lib6_ASV_abundance.qza\
    --i-tables cz_lib7/qiime2/output/cz_lib7_ASV_abundance.qza\
    --i-tables cz_lib8/qiime2/output/cz_lib8_ASV_abundance.qza\
    --o-merged-table ${path_results}/ASV_abundance_merged.qza

qiime feature-table merge-seqs \
    --i-data cz_lib2/qiime2/output/cz_lib2_ASV_seqs.qza \
    --i-data cz_lib3/qiime2/output/cz_lib3_ASV_seqs.qza \
    --i-data cz_lib4/qiime2/output/cz_lib4_ASV_seqs.qza \
    --i-data cz_lib5/qiime2/output/cz_lib5_ASV_seqs.qza \
    --i-data cz_lib6/qiime2/output/cz_lib6_ASV_seqs.qza \
    --i-data cz_lib7/qiime2/output/cz_lib7_ASV_seqs.qza \
    --i-data cz_lib8/qiime2/output/cz_lib8_ASV_seqs.qza \
    --o-merged-data ${path_results}/ASV_seqs_merged.qza

qiime feature-table merge-taxa \
    --i-data cz_lib2/qiime2/output/cz_lib2_ASV_taxonomy.qza \
    --i-data cz_lib3/qiime2/output/cz_lib3_ASV_taxonomy.qza \
    --i-data cz_lib4/qiime2/output/cz_lib4_ASV_taxonomy.qza \
    --i-data cz_lib5/qiime2/output/cz_lib5_ASV_taxonomy.qza \
    --i-data cz_lib6/qiime2/output/cz_lib6_ASV_taxonomy.qza \
    --i-data cz_lib7/qiime2/output/cz_lib7_ASV_taxonomy.qza \
    --i-data cz_lib8/qiime2/output/cz_lib8_ASV_taxonomy.qza \
    --o-merged-data ${path_results}/ASV_taxonomy_merged.qza


qiime taxa filter-table \
	--i-table ASV_abundance_merged.qza \
	--i-taxonomy ASV_taxonomy_merged.qza \
	--p-include 'o__' \
	--p-mode 'contains' \
	--p-exclude 'mitochondria,chloroplast,o__;' \
	--o-filtered-table ASV_abundance_filtered.qza

##############################################################
# STAGE 4: EXPORT
#############################################################

qiime tools export \
--input-path ASV_abundance_filtered.qza \
--output-path results

qiime tools export \
--input-path ASV_seqs_merged.qza \
--output-path results

qiime tools export \
--input-path ASV_taxonomy_merged.qza \
--output-path results