#!/bin/bash

# =============================================================================
# Advanced Bioinformatics 2026 Assessment
# NGS Pipeline Script - Alternative Pipeline (Bowtie2 + Freebayes)
# Author: K24121615
# Usage: bash pipeline_bowtie2.sh
# =============================================================================

# =============================================================================
# STEP 2.1 - TOOL INSTALLATION
# =============================================================================

cd ~/
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh
bash ./Anaconda3-2022.10-Linux-x86_64.sh
source ~/.bashrc

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -y samtools
conda install -y bwa
conda install -y freebayes
conda install -y picard
conda install -y bedtools
conda install -y trimmomatic
conda install -y fastqc
conda install -y snpeff

sudo apt install libvcflib-tools
sudo apt install bowtie2

tar -zxvf annovar.latest.tar.gz

cd ~/annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/

cd ~/anaconda3/share/snpeff-5.1-0/data/
wget https://snpeff-public.s3.amazonaws.com/databases/v5_0/snpEff_v5_0_hg19.zip
unzip snpEff_v5_0_hg19.zip
rm ~/anaconda3/share/snpeff-5.1-0/data/hg19/nextProt.bin

# =============================================================================
# STEP 2.1 - PROJECT SETUP AND INPUT FILE DOWNLOAD
# =============================================================================

mkdir -p ~/ngs_assignment_1/dnaseq/data/untrimmed_fastq
mkdir -p ~/ngs_assignment_1/dnaseq/data/trimmed_fastq
mkdir -p ~/ngs_assignment_1/dnaseq/data/aligned_data
mkdir -p ~/ngs_assignment_1/dnaseq/data/reference
mkdir -p ~/ngs_assignment_1/dnaseq/meta
mkdir -p ~/ngs_assignment_1/dnaseq/results
mkdir -p ~/ngs_assignment_1/dnaseq/logs

cd ~/ngs_assignment_1/dnaseq/data/untrimmed_fastq
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed \
  -P ~/ngs_assignment_1/dnaseq/data/

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz \
  -P ~/ngs_assignment_1/dnaseq/data/reference/

zcat ~/ngs_assignment_1/dnaseq/data/reference/hg19.fa.gz \
  > ~/ngs_assignment_1/dnaseq/data/reference/hg19.fa

# =============================================================================
# STEP 2.2 - PRE-ALIGNMENT QC
# =============================================================================

cd ~/ngs_assignment_1/dnaseq/data/untrimmed_fastq
fastqc -t 4 *.fastq.gz

mkdir -p ~/ngs_assignment_1/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_assignment_1/dnaseq/results/fastqc_untrimmed_reads/

cd ~/ngs_assignment_1/dnaseq/results/fastqc_untrimmed_reads/
for zip in *.zip
do
  unzip $zip
done
cat */summary.txt > ~/ngs_assignment_1/dnaseq/logs/fastqc_untrimmed_summaries.txt

trimmomatic PE \
  -threads 4 \
  -phred33 \
  ~/ngs_assignment_1/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
  ~/ngs_assignment_1/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
  -baseout ~/ngs_assignment_1/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.40-hdfd78af_0/share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

fastqc -t 4 \
  ~/ngs_assignment_1/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P \
  ~/ngs_assignment_1/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P

mkdir -p ~/ngs_assignment_1/dnaseq/results/fastqc_trimmed_reads
mv ~/ngs_assignment_1/dnaseq/data/trimmed_fastq/*fastqc* \
  ~/ngs_assignment_1/dnaseq/results/fastqc_trimmed_reads/

# =============================================================================
# STEP 2.3 - ALIGNMENT (BOWTIE2 - REPLACES BWA MEM)
# =============================================================================

samtools faidx ~/ngs_assignment_1/dnaseq/data/reference/hg19.fa

mkdir -p ~/ngs_assignment_1/dnaseq/data/reference/bowtie2_index

bowtie2-build \
  ~/ngs_assignment_1/dnaseq/data/reference/hg19.fa \
  ~/ngs_assignment_1/dnaseq/data/reference/bowtie2_index/hg19

bowtie2 -p 4 \
  --rg-id HWI-D0011.50.H7AP8ADXX.1.NGS0001 \
  --rg "SM:NGS0001" \
  --rg "PL:ILLUMINA" \
  --rg "LB:nextera-NGS0001-blood" \
  -x ~/ngs_assignment_1/dnaseq/data/reference/bowtie2_index/hg19 \
  -1 ~/ngs_assignment_1/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P \
  -2 ~/ngs_assignment_1/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P \
  | samtools view -h -b \
  | samtools sort -o ~/ngs_assignment_1/dnaseq/data/aligned_data/NGS0001_bowtie2_sorted.bam

samtools index ~/ngs_assignment_1/dnaseq/data/aligned_data/NGS0001_bowtie2_sorted.bam

cd ~/ngs_assignment_1/dnaseq/data/aligned_data

picard MarkDuplicates \
  I=NGS0001_bowtie2_sorted.bam \
  O=NGS0001_bowtie2_sorted_marked.bam \
  M=marked_dup_metrics_bowtie2.txt

samtools index NGS0001_bowtie2_sorted_marked.bam
rm NGS0001_bowtie2_sorted.bam

samtools view -F 1796 -q 20 \
  -o NGS0001_bowtie2_sorted_filtered.bam \
  NGS0001_bowtie2_sorted_marked.bam

samtools index NGS0001_bowtie2_sorted_filtered.bam
rm NGS0001_bowtie2_sorted_marked.bam

samtools flagstat NGS0001_bowtie2_sorted_filtered.bam \
  > ~/ngs_assignment_1/dnaseq/logs/NGS0001_bowtie2_flagstat.txt

samtools idxstats NGS0001_bowtie2_sorted_filtered.bam \
  > ~/ngs_assignment_1/dnaseq/logs/NGS0001_bowtie2_idxstats.txt

picard CollectInsertSizeMetrics \
  I=NGS0001_bowtie2_sorted_filtered.bam \
  O=~/ngs_assignment_1/dnaseq/logs/NGS0001_bowtie2_insert_size_metrics.txt \
  H=~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_insert_size_histogram.pdf

bedtools coverage \
  -a ~/ngs_assignment_1/dnaseq/data/annotation_sorted.bed \
  -b ~/ngs_assignment_1/dnaseq/data/aligned_data/NGS0001_bowtie2_sorted_filtered.bam \
  -sorted \
  -g ~/ngs_assignment_1/dnaseq/data/reference/hg19.genome \
  > ~/ngs_assignment_1/dnaseq/logs/NGS0001_bowtie2_coverage.txt

# =============================================================================
# STEP 2.4 - VARIANT CALLING AND FILTERING
# =============================================================================

freebayes \
  --bam ~/ngs_assignment_1/dnaseq/data/aligned_data/NGS0001_bowtie2_sorted_filtered.bam \
  --fasta-reference ~/ngs_assignment_1/dnaseq/data/reference/hg19.fa \
  --vcf ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2.vcf

bgzip ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2.vcf
tabix -p vcf ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2.vcf.gz

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
  ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2.vcf.gz \
  > ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered.vcf

bedtools intersect -header -wa \
  -a ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered.vcf \
  -b ~/ngs_assignment_1/dnaseq/data/annotation.bed \
  > ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered_annotated.vcf

bgzip ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered_annotated.vcf
tabix -p vcf ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered_annotated.vcf.gz

# =============================================================================
# STEP 2.5 - VARIANT ANNOTATION AND PRIORITISATION
# =============================================================================

cd ~/annovar

./convert2annovar.pl -format vcf4 \
  ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered_annotated.vcf.gz \
  > ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered_annotated.avinput

./table_annovar.pl \
  ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered_annotated.avinput humandb/ \
  -buildver hg19 \
  -out ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_annotated \
  -remove \
  -protocol refGene,clinvar_20180603,exac03,dbnsfp31a_interpro,avsnp150 \
  -operation g,f,f,f,f \
  -otherinfo -nastring . -csvout

/usr/lib/jvm/java-11-openjdk-amd64/bin/java \
  -jar ~/anaconda3/share/snpeff-5.1-0/snpEff.jar hg19 \
  ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_filtered_annotated.vcf.gz \
  > ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_snpeff.vcf

head -1 ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_annotated.hg19_multianno.csv \
  > ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_exonic_novel.csv

grep "exonic" ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_annotated.hg19_multianno.csv | \
  grep -v "ncRNA_exonic" | \
  awk -F',' '$25 == "." || $25 == "\".\""' \
  >> ~/ngs_assignment_1/dnaseq/results/NGS0001_bowtie2_exonic_novel.csv
