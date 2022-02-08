## software versions used in this pipeline:

# sratoolkit v2.9.2 https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-centos_linux64.tar.gz
# bwa-mem v0.7.15 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download
# picardtools v2.2.2 https://github.com/broadinstitute/picard/releases/tag/2.2.2
# GATK v3.7 https://console.cloud.google.com/storage/browser/gatk-software/package-archive !Google Cloud Account necessary!
# samtools v1.3.1 https://github.com/samtools/samtools/releases/tag/1.3.1
# varscan v2.4.2 https://github.com/dkoboldt/varscan/releases/tag/2.4.2
# annovar v2015Dec14 https://www.openbioinformatics.org/annovar/annovar_download_form.php FREE for academic use, registration required, ask the author for archived versions

## all software should be added to PATH

#!/bin/bash 
# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")

mkdir -p \
  $SCRIPTPATH/input_data/in_facility/{00_metrics,01_raw_reads,02_aligned_bam,03_inter_bam,04_processed_bam,05_mpileup,06_varscan,07_annovar,databases}

project_folder="/home/m3fga/ngs-filer/Franz/2022-01-26_shell_script_mutation_calling_paper/2022_mouse_mutation_calling"

## databases used in this pipeline:

# create bwa index of mm10.fa
# cd ${project_folder}/databases
# wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz
# unpigz mm10.fa.gz
# bwa index mm10.fa
# samtools faidx mm10.fa
# picard CreateSequenceDictionary R=mm10.fa

bwa_index=${project_folder}"/databases/mm10"
reference=${project_folder}"/databases/mm10.fa"

bwa_index="/m3-ngs/Archiv/databases/mm10/reference_bwa_Maria/mm10"
reference="/m3-ngs/Archiv/databases/mm10/reference_bwa_Maria/mm10.fa"
dbSNP137="/m3-ngs/Archiv/databases/mm10/snp_mm10_137/snp_order10-19_1-9_MXY_sort.vcf"

## sample IDs:

cd ${project_folder}

SRA_run=( $( cut -f 6 SRA_infos_anonymysed.txt ) )
bam_ID=( $( cut -f 5 SRA_infos_anonymysed.txt ) )

germ=( $( cut -f 1 matched_samples.txt ) )
tumor=( $( cut -f 2 matched_samples.txt ) )

## analysis pipeline:

for i in {1..51}; do
for i in 5 6 8; do
  
  cd ~/ngs-filer/programs/sratoolkit.2.9.2-centos_linux64/bin
  ./prefetch.2.9.2 ${SRA_run[i]}
  ./fasterq-dump.2.9.2 --outdir ${project_folder}/01_raw_reads --split-files ${SRA_run[i]}
  rm ~/ncbi/public/sra/${SRA_run[i]}.sra # this file might be in a different location on your system. You can remove it, it is not needed anymore.

  cd ${project_folder}
  pigz 01_raw_reads/${SRA_run[i]}*.fastq
  
  bwa mem \
    -M -t 8 \
    -R '@RG\tID:mouse\tPL:Illumina\tLB:Exome\tPU:none\tSM:'${bam_ID[i]} \
    ${reference} \
    01_raw_reads/${SRA_run[i]}"_1.fastq.gz" \
    01_raw_reads/${SRA_run[i]}"_2.fastq.gz" \
    | samtools sort -@8 -o 02_aligned_bam/${bam_ID[i]}.sorted.bam -

  samtools index 02_aligned_bam/${bam_ID[i]}.sorted.bam

  picard MarkDuplicates \
    INPUT=02_aligned_bam/${bam_ID[i]}.sorted.bam \
    OUTPUT=03_inter_bam/${bam_ID[i]}.sorted.marked.bam \
    METRICS_FILE=00_metrics/${bam_ID[i]}.sorted.marked.metrics \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

  eval /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.312.b07-1.el7_9.x86_64/jre/bin/java -Xmx32g -jar ~/ngs-filer/programs/GenomeAnalysisTK.jar -T RealignerTargetCreator \
    -R ${reference} \
    -I 03_inter_bam/${bam_ID[i]}.sorted.marked.bam \
    -o 03_inter_bam/${bam_ID[i]}.sorted.marked.list

  eval /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.312.b07-1.el7_9.x86_64/jre/bin/java -Xmx32g -jar ~/ngs-filer/programs/GenomeAnalysisTK.jar -T IndelRealigner \
    -targetIntervals 03_inter_bam/${bam_ID[i]}.sorted.marked.list \
    -R ${reference} \
    -I 03_inter_bam/${bam_ID[i]}.sorted.marked.bam \
    -o 03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.bam

  picard FixMateInformation \
    INPUT=03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.bam \
    OUTPUT=03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.bam \
    SO=coordinate \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true \
    TMP_DIR="temp/"

  picard AddOrReplaceReadGroups \
    INPUT=03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.bam \
    OUTPUT=03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.fixedRG.bam \
    SORT_ORDER=coordinate RGID=mouse RGLB=Exome RGPL=Illumina RGPU=none RGSM=${bam_ID[i]} \
    CREATE_INDEX=true

  eval /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.312.b07-1.el7_9.x86_64/jre/bin/java -Xmx32g -jar ~/ngs-filer/programs/GenomeAnalysisTK.jar -T BaseRecalibrator \
    -I 03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.fixedRG.bam \
    -R ${reference} \
    -knownSites ${dbSNP137} \
    -o 03_inter_bam/${bam_ID[i]}.recal_data.table

  eval /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.312.b07-1.el7_9.x86_64/jre/bin/java -Xmx32g -jar ~/ngs-filer/programs/GenomeAnalysisTK.jar -T PrintReads \
    -BQSR 03_inter_bam/${bam_ID[i]}.recal_data.table \
    -I 03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.fixedRG.bam \
    -R ${reference} \
    -o 04_processed_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.fixedRG.recalBam.bam

done

for i in {0 .. 41}; do
for i in 4 8; do
  samtools mpileup \
    -B -q 1 \
    -f ${reference} \
    04_processed_bam/${germ[$i]}".sorted.marked.realigned.fixed.fixedRG.recalBam.bam" \
    04_processed_bam/${tumor[$i]}".sorted.marked.realigned.fixed.fixedRG.recalBam.bam" \
    > 05_mpileup/${germ[$i]}"vs"${tumor[$i]}.mpileup

  varscan somatic \
    05_mpileup/${germ[$i]}"vs"${tumor[$i]}.mpileup \
    06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan \
    --mpileup 1 \
    --min-coverage-normal 5 \
    --min-coverage-tumor 5 \
    --min-var-freq 0.05 \
    --somatic-p-value 0.05 \
    --output-vcf 1 \
    --strand-filter 1

  varscan processSomatic \
    06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp.vcf
  varscan processSomatic \
    06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.indel.vcf

  varscan somaticFilter \
    06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp.Somatic.hc.vcf \
    -indel-file 06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.indel.vcf \
    -output-file 06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp.Somatic.hc.filter.vcf

  grep -v "#" 06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.indel.Somatic.hc.vcf \
    > 06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.indel.Somatic.hc.woheader.vcf
  cat 06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp.Somatic.hc.filter.vcf \
    06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.indel.Somatic.hc.woheader.vcf \
    > 06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf

  cd /m3-ngs/databases/mm10/annovar

  convert2annovar.pl \
    --format vcf4old \
    -includeinfo \
    06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf \
    > 07_annovar/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf.annovar

  table_annovar.pl \
    07_annovar/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf.annovar \
    /m3-ngs/Archiv/databases/mm10/annovar/mousedb/ \
    -buildver mm10 \
    -out 07_annovar/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf.annovar.out \
    -protocol RefGene,dbsnp137 \
    -operation g,f \
    -otherinfo \
    -csvout \
    -remove

done
