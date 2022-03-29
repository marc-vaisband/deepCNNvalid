## software versions used in this pipeline:

# sratoolkit v2.9.2 https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-centos_linux64.tar.gz
# bwa-mem v0.7.15 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.15.tar.bz2/download
# picardtools v2.2.2 https://github.com/broadinstitute/picard/releases/tag/2.2.2
# GATK v3.7 https://console.cloud.google.com/storage/browser/gatk-software/package-archive !Google Cloud Account necessary!
# samtools v1.3.1 https://github.com/samtools/samtools/releases/tag/1.3.1
# varscan v2.4.2 https://github.com/dkoboldt/varscan/releases/tag/2.4.2
# annovar v2015Dec14 https://www.openbioinformatics.org/annovar/annovar_download_form.php FREE for academic use, registration required, ask the author for archived versions

## all software should be added to PATH

## create a project folder structure similar to this:

mkdir -p \
  data/in_facility/{00_metrics,01_raw_reads,02_aligned_bam,03_inter_bam,04_processed_bam,05_mpileup,06_varscan,07_annovar,databases}

project_folder="data/in_facility"

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

# we used build 137 in our analysis, which is not available on the NCBI ftp server anymore
# if you detect large discrepancies in called mutations, we can provide you with our archived dbsnp137 vcf file
# please contact the authors in this case
#
# download mouse dbsnp database from ftp://ftp.ncbi.nih.gov/snp/.redesign/pre_build152/organisms/archive/mouse_10090/VCF
# you can either download the "00-All.vcf.gz" file and sort the vcf like this: chr1,chr2,...,chr19,chrX
# or download the single chromosomes and concatenate them in that order
# In any case, remove additional contigs e.g. chr_Un and chr_Alt, use only the 19 full chromosomes + chrX
# then save the unzipped vcf file in the databases folder

dbsnp_file=${project_folder}"/databases/snp_natural_sorted.vcf"

# download annovar databases
# annotate_variation.pl --buildver mm10 --downdb RefGene ${project_folder}"/databases/mm10/annovar/mousedb/"
# annotate_variation.pl --buildver mm10 --downdb dbsnp137 ${project_folder}"/databases/mm10/annovar/mousedb/"

annovar_db=${project_folder}"/databases/mm10/annovar/mousedb/"

## sample IDs:

cd ${project_folder}

SRA_run=( $( cut -d ";" -f 2 SRA_infos.csv ) )
bam_ID=( $( cut -d ";" -f 1 SRA_infos.csv ) )

germ=( $( cut -f 2 relative_annovarlist.csv ) )
tumor=( $( cut -f 3 relative_annovarlist.csv ) )

## analysis pipeline:

for i in {1..51}; do
  
  prefetch.2.9.2 ${SRA_run[i]}
  fasterq-dump.2.9.2 --outdir ${project_folder}/01_raw_reads --split-files ${SRA_run[i]}
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

  gatk -T RealignerTargetCreator \
    -R ${reference} \
    -I 03_inter_bam/${bam_ID[i]}.sorted.marked.bam \
    -o 03_inter_bam/${bam_ID[i]}.sorted.marked.list

  gatk -T IndelRealigner \
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

  gatk -T BaseRecalibrator \
    -I 03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.fixedRG.bam \
    -R ${reference} \
    -knownSites ${dbSNP137} \
    -o 03_inter_bam/${bam_ID[i]}.recal_data.table

  gatk -T PrintReads \
    -BQSR 03_inter_bam/${bam_ID[i]}.recal_data.table \
    -I 03_inter_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.fixedRG.bam \
    -R ${reference} \
    -o 04_processed_bam/${bam_ID[i]}.sorted.marked.realigned.fixed.fixedRG.recalBam.bam

done

for i in {0 .. 41}; do

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

  convert2annovar.pl \
    --format vcf4old \
    -includeinfo \
    06_varscan/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf \
    > 07_annovar/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf.annovar

  table_annovar.pl \
    07_annovar/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf.annovar \
    ${annovar_db} \
    -buildver mm10 \
    -out 07_annovar/${germ[$i]}"vs"${tumor[$i]}.VarScan.snp_indel.Somatic.hc.filter.vcf.annovar.out \
    -protocol RefGene,dbsnp137 \
    -operation g,f \
    -otherinfo \
    -csvout \
    -remove

done
