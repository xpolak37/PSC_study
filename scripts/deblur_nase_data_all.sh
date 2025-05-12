conda activate /home/povp/ppola/conda_envs/qiime2-amplicon-2024.2
export TMPDIR=/home/povp/data

path=/home/povp/ppola/runs/
path_vysledky=/home/povp/ppola/norsko_nase_data/
path_base=/home/povp/ppola/
path_bbmap=/home/povp/ppola/bbmap/


cd $path

# removing adapters and phiX
for lib_name in knihovna8
do
	cd ${path}/${lib_name}/
	mkdir ${path_vysledky}/${lib_name}
	mkdir ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk

	for read1 in *R1_001.fastq*
	do	
	read2=$(echo $read1| sed 's/R1_/R2_/')
	${path_bbmap}bbduk.sh \
	ref=${path_bbmap}resources/adapters.fa \
	in1=${read1} in2=${read2} out=stdout.fq \
	k=23 hdist=1 tbo cf=TRUE ftm=5 \
	2> ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/remove_adaptors.log | \
	${path_bbmap}bbduk.sh \
	ref=${path_bbmap}resources/phix174_ill.ref.fa.gz \
	in=stdin.fq int=t \
	out1=${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/noadaptors_nophix_${read1} \
	out2=${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/noadaptors_nophix_${read2} \
	k=31 hdist=1 2> ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/remove_phiX.log
	done
done

for lib_name in knihovna8
do
cd /home/povp/ppola/runs/${lib_name}/${lib_name}/ 
echo -e "sampleid\tbases_in_R1\treads_in_R1\tavg_read_length_in_R1" > /home/povp/ppola/runs/readcount_raw_${lib_name}.tsv 
for i in *_R1_* | cut -f 1)
do
echo -en "$i\t" >> /home/povp/ppola/runs/readcount_raw_${lib_name}.tsv; \
zcat ${i} | paste - - - - | cut -f 2 | perl -ne 'END {if ($. == 0) {print "0\t0\t0\n"} else {print "$c\t$.\t", sprintf("%.0f", $c/$.), "\n"}} $c+=(length($_)-1)' - >> /home/povp/ppola/runs/readcount_raw_${lib_name}.tsv; \
done 
done

# CUTADAPT

for lib_name in knihovna8
	do
	mkdir ${path_vysledky}/${lib_name}/trimmed

	cd ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk

	for read1 in *R1_001.fastq*
	do
		read2=$(echo $read1| sed 's/R1_/R2_/')
		sample1=$(echo $read1| sed 's/noadaptors_nophix_//')
		sample2=$(echo $read2| sed 's/noadaptors_nophix_//')
		command1="cutadapt "
		cd $path_base
		for item in $(cat mikrobiom_primery.csv) 
		do
			FWD=$(echo $item | cut -f2 -d ",")
			REV=$(echo $item | cut -f5 -d ",")

			command1+="-g ^${FWD} -G ^${REV} "
		done
		command1+="--discard-untrimmed -j 5 -o ${path_vysledky}/${lib_name}/trimmed/trimmed_${sample1} -p ${path_vysledky}/${lib_name}/trimmed/trimmed_${sample2} \
		${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/${read1} ${path_vysledky}/${lib_name}/remove_adapters_and_phiX_bbduk/${read2}"

		eval $command1
		
	done
done

# MERGING 
for lib_name in knihovna8
do
	cd ${path_vysledky}/${lib_name}/trimmed
	mkdir ${path_vysledky}/${lib_name}/merge_pairs_bbmerge

	for read1 in *R1_001.fastq*
	do

		read2=$(echo $read1| sed 's/R1/R2/')
		sample1=$(echo $read1| sed 's/trimmed_//')
		sample2=$(echo $read2| sed 's/trimmed_//')
		sample=$(echo $read1| sed 's/_L001_R1_001.fastq.gz//')
		${path_bbmap}bbmerge.sh \
		in1=${path_vysledky}/${lib_name}/trimmed/${read1} \
		in2=${path_vysledky}/${lib_name}/trimmed/${read2} \
	out=stdout.fq qtrim=r trimq=15 maxlength=500 mininsert=350 ihist=../merge_pairs_bbmerge/${sample}_insert_histogram.ihist \
	| ${path_bbmap}reformat.sh -Xmx20g -int=f -maxns=0 \
	in=stdin.fq out=${path_vysledky}/${lib_name}/merge_pairs_bbmerge/${sample}_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz \
	lhist=${path_vysledky}/${lib_name}/merge_pairs_bbmerge/${sample}_histogram_after_N_filter.txt;
	done
done

# ORIENTING
conda activate /home/povp/ppola/conda_envs/vsearch

for lib_name in knihovna8
do
	cd ${path_vysledky}/${lib_name}/merge_pairs_bbmerge
	mkdir ${path_vysledky}/${lib_name}/oriented
	for read in *mergedpairs.fastq.gz
	do
		sample=$(echo $read| sed 's/trimmed_//')
		sample=$(echo $sample| sed 's/_noadaptors_nophix_trimmedprimers_mergedpairs.fastq.gz//')
		vsearch --orient $read \
				--db /home/povp/ppola/silva_db/silva_nrr99_v138.1_orienting.udb \
				--fastqout ${path_vysledky}/${lib_name}/oriented/${sample}_oriented.fq \
				--tabbedout ${path_vysledky}/${lib_name}/oriented/orient.tx
	done
done

# 
conda activate /home/povp/ppola/conda_envs/qiime2-amplicon-2024.2

for lib_name in knihovna8
do
	cd ${path_vysledky}/${lib_name}/oriented
	mkdir ../qiime2
	mkdir ../qiime2/output
	echo -e "sample-id\tabsolute-filepath\tdirection" > ../qiime2/manifest_${lib_name}.tsv
	for read in *-oriented.fq; do
		sample=$(echo $read| sed 's/-oriented.fq//')
		echo -e "${sample}\t${path_vysledky}${lib_name}/oriented/${read}\tforward" >>  ../qiime2/manifest_${lib_name}.tsv
	done

	cd ../qiime2/

	qiime tools import \
    	--input-path manifest_${lib_name}.tsv \
    	--type 'SampleData[SequencesWithQuality]' \
    	--input-format SingleEndFastqManifestPhred33V2 \
    	--output-path ${lib_name}.qza

    	qiime demux summarize --i-data ${lib_name}.qza --o-visualization output/${lib_name}_imported_seqs.qzv --p-n 10000

    	qiime deblur denoise-16S \
	--i-demultiplexed-seqs ${lib_name}.qza \
	--p-trim-length 400 \
	--p-sample-stats \
	--o-representative-sequences output/${lib_name}_ASV_seqs.qza \
	--o-table output/${lib_name}_ASV_abundance.qza \
	--o-stats output/${lib_name}_stats_deblur.qza
	
	qiime feature-classifier classify-sklearn \
	--i-reads ${path_vysledky}/${lib_name}/qiime2/output/${lib_name}_ASV_seqs.qza \
	--i-classifier /home/povp/ppola/classifier/silva-138.1-ssu-nr99-319f-806r-classifier.qza \
	--p-n-jobs -1 \
	--o-classification ${path_vysledky}/${lib_name}/qiime2/output/${lib_name}_ASV_taxonomy.qza
	done


mkdir /home/povp/ppola/norsko_nase_data/with_lib8/

qiime feature-table merge \
    --i-tables /home/povp/ppola/norsko_nase_data/knihovna2/qiime2/output/knihovna2_ASV_abundance.qza \
    --i-tables /home/povp/ppola/norsko_nase_data/knihovna3/qiime2/output/knihovna3_ASV_abundance.qza\
    --i-tables /home/povp/ppola/norsko_nase_data/knihovna4/qiime2/output/knihovna4_ASV_abundance.qza\
    --i-tables /home/povp/ppola/norsko_nase_data/knihovna5/qiime2/output/knihovna5_ASV_abundance.qza\
    --i-tables /home/povp/ppola/norsko_nase_data/knihovna6/qiime2/output/knihovna6_ASV_abundance.qza\
    --i-tables /home/povp/ppola/norsko_nase_data/knihovna7/qiime2/output/knihovna7_ASV_abundance.qza\
    --i-tables /home/povp/ppola/norsko_nase_data/knihovna8/qiime2/output/knihovna8_ASV_abundance.qza\
    --o-merged-table /home/povp/ppola/norsko_nase_data/with_lib8/ASV_abundance_merged.qza

qiime feature-table merge-seqs \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna2/qiime2/output/knihovna2_ASV_seqs.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna3/qiime2/output/knihovna3_ASV_seqs.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna4/qiime2/output/knihovna4_ASV_seqs.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna5/qiime2/output/knihovna5_ASV_seqs.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna6/qiime2/output/knihovna6_ASV_seqs.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna7/qiime2/output/knihovna7_ASV_seqs.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna8/qiime2/output/knihovna8_ASV_seqs.qza \
    --o-merged-data /home/povp/ppola/norsko_nase_data/with_lib8/ASV_seqs_merged.qza

qiime feature-table merge-taxa \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna2/qiime2/output/knihovna2_ASV_taxonomy.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna3/qiime2/output/knihovna3_ASV_taxonomy.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna4/qiime2/output/knihovna4_ASV_taxonomy.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna5/qiime2/output/knihovna5_ASV_taxonomy.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna6/qiime2/output/knihovna6_ASV_taxonomy.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna7/qiime2/output/knihovna7_ASV_taxonomy.qza \
    --i-data /home/povp/ppola/norsko_nase_data/knihovna8/qiime2/output/knihovna8_ASV_taxonomy.qza \
    --o-merged-data /home/povp/ppola/norsko_nase_data/with_lib8/ASV_taxonomy_merged.qza

cd /home/povp/ppola/norsko_nase_data/with_lib8/

qiime taxa filter-table \
	--i-table ASV_abundance_merged.qza \
	--i-taxonomy ASV_taxonomy_merged.qza \
	--p-include 'o__' \
	--p-mode 'contains' \
	--p-exclude 'mitochondria,chloroplast,o__;' \
	--o-filtered-table ASV_abundance_filtered.qza

qiime tools export \
--input-path ASV_abundance_filtered.qza \
--output-path vysledky
qiime tools export \
--input-path ASV_seqs_merged.qza \
--output-path vysledky
qiime tools export \
--input-path ASV_taxonomy_merged.qza \
--output-path vysledky


qiime picrust2 full-pipeline \
   --i-table ASV_abundance_filtered.qza \
   --i-seq ASV_seqs_merged.qza \
   --output-dir picrust \
   --p-placement-tool sepp \
   --p-threads 10 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose

#cd picrust

#qiime tools export \
#  --input-path ec_metagenome.qza \
#  --output-path vysledky

#mv vysledky/feature-table.biom vysledky/ec_metagenome.biom

#qiime tools export \
#  --input-path ko_metagenome.qza \
#  --output-path vysledky

#mv vysledky/feature-table.biom vysledky/ko_metagenome.biom

#qiime tools export \
#  --input-path pathway_abundance.qza\
#  --output-path vysledky

#mv vysledky/feature-table.biom vysledky/pathway_abundance.biom