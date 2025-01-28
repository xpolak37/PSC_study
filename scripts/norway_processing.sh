conda activate /home/povp/ppola/conda_envs/qiime2-amplicon-2024.2
path_bbmap=/home/povp/ppola/bbmap/
path=/home/povp/ppola/norsko/
path_base=/home/povp/ppola/
cd $path
for lib_name in lib20 lib21 lib27 lib29 lib30 lib32 lib48
do	
	# Removing Adapters and PhiX

	mkdir $lib_name
	mkdir ${lib_name}/remove_adapters_and_phiX_bbduk/
	${path_bbmap}/bbduk.sh ref=${path_bbmap}/resources/adapters.fa in1=data/${lib_name}_R1.fastq.gz in2=data/${lib_name}_R2.fastq.gz \
	out=stdout.fq k=23 hdist=1 tbo cf=TRUE ftm=5 \
	2> ${lib_name}/remove_adapters_and_phiX_bbduk/remove_adaptors.log | ${path_bbmap}/bbduk.sh \
	ref=${path_bbmap}/resources/phix174_ill.ref.fa.gz in=stdin.fq int=t \
	out1=${lib_name}/remove_adapters_and_phiX_bbduk/${lib_name}_noadaptors_nophix_R1.fastq.gz \
	out2=${lib_name}/remove_adapters_and_phiX_bbduk/${lib_name}_noadaptors_nophix_R2.fastq.gz \
	k=31 hdist=1 2> ${lib_name}/remove_adapters_and_phiX_bbduk/remove_phiX.log

	# Demultiplexing

	ulimit -n 10000
	mkdir ${lib_name}/demultiplex_cutadapt
	cutadapt -j 0 -e 1 --no-indels --discard-untrimmed --action none -g ^file:data/barcodes_V3.fa -G ^file:data/barcodes_V4.fa \
	-o ${lib_name}/demultiplex_cutadapt/{name1}-{name2}.R1.fastq.gz -p ${lib_name}/demultiplex_cutadapt/{name1}-{name2}.R2.fastq.gz \
	${lib_name}/remove_adapters_and_phiX_bbduk/${lib_name}_noadaptors_nophix_R2.fastq.gz ${lib_name}/remove_adapters_and_phiX_bbduk/${lib_name}_noadaptors_nophix_R1.fastq.gz > ${lib_name}/demultiplex_cutadapt/demultiplex_cutadapt.log
	
	# create pattern file for mmv to rename demultiplexed files based on barcodes
	# I dont have mmv and do not have admin rights

	perl -ane 'print "$F[2]-$F[4].R1.fastq.gz\t$F[0]_noadaptors_nophix_R1.fastq.gz\n"' \
	data/mapping_${lib_name}.txt | tail -n +2 > ${lib_name}/pattern_for_renaming_demuxed_files_cutadapt.txt
	perl -ane 'print "$F[2]-$F[4].R2.fastq.gz\t$F[0]_noadaptors_nophix_R2.fastq.gz\n"' \
	data/mapping_${lib_name}.txt | tail -n +2 >> ${lib_name}/pattern_for_renaming_demuxed_files_cutadapt.txt
	cd ${lib_name}/demultiplex_cutadapt
	
	# Read the pattern file line by line 
	while IFS=$'\t' read -r old_pattern new_pattern; do
    	# Use mv command with wildcard to rename the files
    	mv "$old_pattern"* "$new_pattern"
	done < ../pattern_for_renaming_demuxed_files_cutadapt.txt
	
	#remove barcode combinations not in mapping file
	rm V3.*-V4*
	cd ../../

	# Read Counts
	echo -e "sampleid\tbases_in_R1\treads_in_R1\tavg_read_length_in_R1" > ${lib_name}/readcount_after_demultiplex_${lib_name}.tsv
	for i in $(tail -n +2 data/mapping_${lib_name}.txt | cut -f 1); do
	echo -en "$i\t" >> ${lib_name}/readcount_after_demultiplex_${lib_name}.tsv;
	zcat ${lib_name}/demultiplex_cutadapt/${i}_noadaptors_nophix_R1.fastq.gz | paste - - - - | cut -f 2 | perl -ne 'END {if ($. == 0) {print "0\t0\t0\n"} else {print "$c\t$.\t", sprintf("%.0f", $c/$.), "\n"}} $c+=(length($_)-1)' - >> \
	${lib_name}/readcount_after_demultiplex_${lib_name}.tsv;
	done

	# trimming primers - missing for loop

	mkdir ${lib_name}/trim_barcodes_and_primers_cutadapt
	mkdir ${lib_name}/trim_barcodes_and_primers_cutadapt/logs
	for i in $(tail -n +2 data/mapping_${lib_name}.txt | cut -f 1); do
	cutadapt -g ACTCCTACGGGAGGCAGCAG -G GGACTACHVGGGTWTCTAAT -e 0.1 --overlap 20 --discard-untrimmed -m 250 -j 0 \
	-o ${lib_name}/trim_barcodes_and_primers_cutadapt/${i}_noadaptors_nophix_trimmedprimers_R1.fastq.gz \
	-p ${lib_name}/trim_barcodes_and_primers_cutadapt/${i}_noadaptors_nophix_trimmedprimers_R2.fastq.gz \
	${lib_name}/demultiplex_cutadapt/${i}_noadaptors_nophix_R1.fastq.gz ${lib_name}/demultiplex_cutadapt/${i}_noadaptors_nophix_R2.fastq.gz \
	> ${lib_name}/trim_barcodes_and_primers_cutadapt/logs/${i}_cutadapt.log;
	done
	
	# quality trimming, merging
	mkdir ${lib_name}/merge_pairs_bbmerge
	mkdir ${lib_name}/merge_pairs_bbmerge/hist
	mkdir ${lib_name}/merge_pairs_bbmerge/logs
	for i in $(tail -n +2 data/mapping_${lib_name}.txt | cut -f 1); do
	${path_bbmap}/bbmerge.sh in1=${lib_name}/trim_barcodes_and_primers_cutadapt/${i}_noadaptors_nophix_trimmedprimers_R1.fastq.gz \
	in2=${lib_name}/trim_barcodes_and_primers_cutadapt/${i}_noadaptors_nophix_trimmedprimers_R2.fastq.gz \
	out=stdout.fq qtrim=r trimq=15 maxlength=500 mininsert=350 ihist=${lib_name}/merge_pairs_bbmerge/hist/${i}_insert_histogram.ihist 2> ${lib_name}/merge_pairs_bbmerge/logs/${i}_log.txt \
	| ${path_bbmap}/reformat.sh -Xmx20g -int=f -maxns=0 in=stdin.fq out=${lib_name}/merge_pairs_bbmerge/${i}_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz \
	lhist=${lib_name}/merge_pairs_bbmerge/hist/${i}_histogram_after_N_filter.txt 2>> ${lib_name}/merge_pairs_bbmerge/logs/${i}_log.txt
	done
rm	
	# Read Counts

	echo -e "sampleid\tbases\tmerged_contigs\tavg_contig_length" > ${lib_name}/readcount_after_demultiplex_trimmedprimers_mergepairs_${lib_name}.tsv
	for i in $(tail -n +2 data/mapping_${lib_name}.txt | cut -f 1); do
	echo -en "$i\t" >> ${lib_name}/readcount_after_demultiplex_trimmedprimers_mergepairs_${lib_name}.tsv;
	zcat ${lib_name}/merge_pairs_bbmerge/${i}_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz | paste - - - - | cut -f 2 | perl -ne 'END {if ($. == 0) {print "0\t0\t0\n"} else {print "$c\t$.\t", sprintf("%.0f", $c/$.), "\n"}} $c+=(length($_)-1)' - >> ${lib_name}/readcount_after_demultiplex_trimmedprimers_mergepairs_${lib_name}.tsv;
	done
	
	##################### qiime2 ###########################
	
	mkdir ${lib_name}/qiime
	mkdir ${lib_name}/qiime/output

	# importing
	echo -e "sample-id\tabsolute-filepath\tdirection" > ${lib_name}/qiime/manifest_${lib_name}.tsv
	for i in $(tail -n +2 data/mapping_${lib_name}.txt | cut -f 1); do
	echo -e "${i}\t${path}/${lib_name}/merge_pairs_bbmerge/${i}_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz\tforward" >>  ${lib_name}/qiime/manifest_${lib_name}.tsv
	done
	
	qiime tools import \
	--input-path ${lib_name}/qiime/manifest_${lib_name}.tsv \
    	--type 'SampleData[SequencesWithQuality]' \
    	--input-format SingleEndFastqManifestPhred33V2 \
    	--output-path ${lib_name}/qiime/${lib_name}.qza

	qiime demux summarize --i-data ${lib_name}/qiime/${lib_name}.qza --o-visualization ${lib_name}/qiime/${lib_name}_imported_seqs.qzv --p-n 10000
	
	qiime deblur denoise-16S \
	--i-demultiplexed-seqs ${lib_name}/qiime/${lib_name}.qza \
	--p-trim-length 400 \
	--p-sample-stats \
	--min_reads 1 \
	--o-representative-sequences ${lib_name}/qiime_minreads/output/${lib_name}_ASV_seqs.qza \
	--o-table ${lib_name}/qiime_minreads/output/${lib_name}_ASV_abundance.qza \
	--o-stats ${lib_name}/qiime_minreads/output/${lib_name}_stats_deblur.qza \
	--p-jobs-to-start 20

	qiime feature-classifier classify-sklearn \
	--i-reads ${lib_name}/qiime/output/${lib_name}_ASV_seqs.qza \
	--i-classifier ${path_base}/classifier/silva-138.1-ssu-nr99-319f-806r-classifier.qza \
	--o-classification  ${lib_name}/qiime/output/${lib_name}_ASV_taxonomy.qza \
	--p-n-jobs 20

	qiime deblur visualize-stats \
	--i-deblur-stats ${lib_name}/qiime/output/${lib_name}_stats_deblur.qza \
	--o-visualization ${lib_name}/qiime/output/${lib_name}_stats_deblur.qzv

	qiime feature-table tabulate-seqs \
	--i-data ${lib_name}/qiime/output/${lib_name}_ASV_seqs.qza \
	--o-visualization ${lib_name}/qiime/output/${lib_name}_ASV_seqs.qzv

	qiime feature-table summarize \
	--i-table ${lib_name}/qiime/output/${lib_name}_ASV_abundance.qza \
	--o-visualization ${lib_name}/qiime/output/${lib_name}_ASV_abundance.qzv \
	--m-sample-metadata-file data/mapping_${lib_name}.txt

	#qiime diversity core-metrics \
	#--i-table qiime2/ASV_abundance.qza \
	#--p-sampling-depth XXXX \
	#--m-metadata-file data/mapping_${lib_name}.txt
	#--output-dir qiime2/diversity-metrics

done

qiime feature-table merge \
    --i-tables lib20/qiime/output/lib20_ASV_abundance.qza \
    --i-tables lib21/qiime/output/lib21_ASV_abundance.qza \
    --i-tables lib27/qiime/output/lib27_ASV_abundance.qza \
    --i-tables lib29/qiime/output/lib29_ASV_abundance.qza \
    --i-tables lib30/qiime/output/lib30_ASV_abundance.qza \
    --i-tables lib32/qiime/output/lib32_ASV_abundance.qza \
    --i-tables lib48/qiime/output/lib48_ASV_abundance.qza \
    --o-merged-table ASV_abundance_merged.qza

qiime feature-table merge-seqs \
    --i-data lib20/qiime/output/lib20_ASV_seqs.qza \
    --i-data lib21/qiime/output/lib21_ASV_seqs.qza \
    --i-data lib27/qiime/output/lib27_ASV_seqs.qza \
    --i-data lib29/qiime/output/lib29_ASV_seqs.qza \
    --i-data lib30/qiime/output/lib30_ASV_seqs.qza \
    --i-data lib32/qiime/output/lib32_ASV_seqs.qza \
    --i-data lib48/qiime/output/lib48_ASV_seqs.qza \
    --o-merged-data ASV_seqs_merged.qza

qiime feature-table merge-taxa \
    --i-data lib20/qiime/output/lib20_ASV_taxonomy.qza \
    --i-data lib21/qiime/output/lib21_ASV_taxonomy.qza \
    --i-data lib27/qiime/output/lib27_ASV_taxonomy.qza \
    --i-data lib29/qiime/output/lib29_ASV_taxonomy.qza \
    --i-data lib30/qiime/output/lib30_ASV_taxonomy.qza \
    --i-data lib32/qiime/output/lib32_ASV_taxonomy.qza \
    --i-data lib48/qiime/output/lib48_ASV_taxonomy.qza \
    --o-merged-data ASV_taxonomy_merged.qza

qiime feature-table filter-samples \
	--i-table ASV_abundance_merged.qza \
	--m-metadata-file data/562_sampleids.tsv \
	--o-filtered-table ASV_abundance_562.qza

qiime taxa filter-table \
	--i-table ASV_abundance_562.qza \
	--i-taxonomy ASV_taxonomy_merged.qza \
	--p-include 'o__' \
	--p-mode 'contains' \
	--p-exclude 'mitochondria,chloroplast,o__;' \
	--o-filtered-table ASV_abundance_562_filtered.qza

qiime tools export \
  --input-path ASV_abundance_562_filtered.qza \
  --output-path results
qiime tools export \
  --input-path ASV_seqs_merged.qza \
  --output-path results
qiime tools export \
  --input-path ASV_taxonomy_merged.qza \
  --output-path results
