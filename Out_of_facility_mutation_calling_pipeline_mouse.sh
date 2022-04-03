## software versions used in this pipeline:

# picardtools v2.21.9 https://github.com/broadinstitute/picard/releases/tag/2.21.9
# GATK v3.7 https://console.cloud.google.com/storage/browser/gatk-software/package-archive !Google Cloud Account necessary!
# GATK v4.1.7.0 https://github.com/broadinstitute/gatk/releases/tag/4.1.7.0
# samtools v1.10 https://github.com/samtools/samtools/releases/tag/1.10
# varscan v2.4.4 https://github.com/dkoboldt/varscan/releases/tag/2.4.4
# annovar v2017July17 https://www.openbioinformatics.org/annovar/annovar_download_form.php FREE for academic use, registration required, ask the author for archived versions

## all software should be added to PATH
## except GATK v3.7:
## GATK="~/PATH/TO/GenomeAnalysisTK.jar" #GATK 3.7-0-gcfedb67 should be stored as a variable and needs java 1.8 to run
## java1="/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.282.b08-2.el8_3.x86_64/jre/bin/java" #java 1.8 should be in a similar location in your system

## create a project folder structure similar to this:

mkdir -p \
  data/kotani/{00_metrics,01_downloaded_bam,02_inter_bam,03_processed_bam,04_mpileup,05_varscan,06_annovar,databases}

project_folder="data/kotani"

## databases used in this pipeline:

# create bwa index of mm10.fa
# cd ${project_folder}/databases
# wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz
# unpigz mm10.fa.gz
# bwa index mm10.fa
# samtools faidx mm10.fa
# picard CreateSequenceDictionary R=mm10.fa

reference=${project_folder}"/databases/mm10.fa"

# we used dbsnp build 142 in our analysis, which is not available on the NCBI ftp server anymore
# if you detect large discrepancies in called mutations, we can provide you with our archived dbsnp142 vcf file
# please contact the authors in this case
#
# download mouse dbsnp database from ftp://ftp.ncbi.nih.gov/snp/.redesign/pre_build152/organisms/archive/mouse_10090/VCF
# you can either download the "00-All.vcf.gz" file and sort the vcf like this: chr1,chr2,...,chr19,chrX
# or download the single chromosomes and concatenate them in that order
# In any case, remove additional contigs e.g. chr_Un and chr_Alt, use only the 19 full chromosomes + chrX
# then save the unzipped vcf file in the databases folder

dbsnp_file=${project_folder}"/databases/snp_natural_sorted.vcf"
# grep commented header lines and "INDEL" and save in new file:
dbsnp_file_indels=${project_folder}"/databases/snp_natural_sorted.INDELS.vcf"

# get the SureSelect XT Mouse All Exon V2 target bed file from https://earray.chem.agilent.com/suredesign/ and liftover to mm10 using this webpage https://genome.ucsc.edu/cgi-bin/hgLiftOver
covered_lib_file=${project_folder}"/databases/covered_library_mm10.bed"
annovar_db="/home/m3fga/limcr-ngs/databases/mm10/annovar/mousedb/"

# download annovar databases
# annotate_variation.pl --buildver mm10 --downdb RefGene ${project_folder}"/databases/mm10/annovar/mousedb/"
# annotate_variation.pl --buildver mm10 --downdb dbsnp142 ${project_folder}"/databases/mm10/annovar/mousedb/"

annovar_db=${project_folder}"/databases/mm10/annovar/mousedb/"

## download sequence data from ENA

cd ${project_folder}/01_downloaded_bam
bash Out_of_facility_download_list.sh

mouse_id=( $( cut -f 1 ENA_mouse_list.txt ) )
sample_id=( $( cut -f 1 ENA_sample_list.txt ) )

for i in {0..83}; do

  samtools index ${project_folder}/01_downloaded_bam/${sample_id[i]}.bam

done

###################################################################################################
# INDEL REALIGNMENT
###################################################################################################

#make interval list for ALL files:
  
cd ${project_folder}/01_downloaded_bam

ls *.bam | sed s/\\s/\\n/g | sed "/^$/d" > bam_list
  
eval $java1 -Xmx32g -jar "$GATK" -T RealignerTargetCreator \
  -nt 48 \
  -R "${reference}" \
  -known "$dbsnp_file_indels" \
  -o sorted.bam.indels.list \
  -I bam.list

#apply IndelRealignment on matched germline-tumor pairs:

for j in 0 7 14 21 28 35; do

  for i in $(eval echo "{$j..$(expr $j + 6)}"); do

    eval $java1 -jar "$GATK" -T IndelRealigner \
    -I "${mouse_id[$i]}"_tail.bam \
    -I "${mouse_id[$i]}"_BM.bam \
    -R $reference \
    -known $dbsnp_file_indels \
    -targetIntervals sorted.bam.indels.list \
    -maxReads 1000000 \
    -nWayOut .realigned.bam
    
  done
  
  wait

done

mv *.realigned.ba* ${project_folder}/02_inter_bam

###################################################################################################
# BASE QUALITY SCORE RECALIBRATION (BQSR)
###################################################################################################

cd ${project_folder}

#Recalculate BQS:

for j in 0 14 28 42 56 70; do

  for i in $(eval echo "{$j..$(expr $j + 13)}"); do
  
    gatk BaseRecalibrator \
      -I 02_inter_bam/"${sample_id[$i]}".realigned.bam \
      -R $reference \
      --known-sites $dbsnp_file \
      -L $covered_lib_file \
      -O 02_inter_bam/"${sample_id[$i]}".realigned.recalBam.table
  
    done
	
  wait
  
done

#Apply recalculated BQS to new bam file:

for j in 0 14 28 42 56 70; do

  for i in $(eval echo "{$j..$(expr $j + 13)}"); do
  
    gatk ApplyBQSR \
      -I 02_inter_bam/"${sample_id[$i]}".realigned.bam \
      -R $reference \
      --bqsr-recal-file 02_inter_bam/"${sample_id[$i]}".realigned.recalBam.table \
      -O 03_processed_bam/"${sample_id[$i]}".realigned.recalBam.bam
 
    done
	
  wait

done

###################################################################################################
# ADDITIONAL BAM PROCeSSING
###################################################################################################

cd ${project_folder}

# repair NM Tags (messed by bwa mem)

for j in 0 14 28 42 56 70; do

  for i in $(eval echo "{$j..$(expr $j + 13)}"); do

    samtools calmd -b \
    03_processed_bam/"${sample_id[$i]}".realigned.recalBam.bam \
    $reference \
    > 03_processed_bam/"${sample_id[$i]}".realigned.recalBam.calmd.bam &
    
    done
	
  wait
  
done

# index new bam files

for j in 0 14 28 42 56 70; do

  for i in $(eval echo "{$j..$(expr $j + 13)}"); do
  
    samtools index 03_processed_bam/"${sample_id[$i]}".realigned.recalBam.calmd.bam
    
    done
	
  wait
  
done

###################################################################################################
# GENERATE MPILEUP FILES
###################################################################################################

cd ${project_folder}

for j in 0 14 28; do

  for i in $(eval echo "{$j..$(expr $j + 13)}"); do

  samtools mpileup -B -q 1 -f "$reference" \
    03_processed_bam/"${mouse_id[$i]}"_tail.realigned.recalBam.calmd.bam \
    03_processed_bam/"${mouse_id[$i]}"_BM.realigned.recalBam.calmd.bam \
    > 04_mpileup/"${mouse_id[$i]}".mpileup &
    
    done
	
  wait
  
done

cd 04_mpileup/
pigz -f -k *.mpileup

###################################################################################################
# VARSCAN MUTATION CALLING AND FILTERING
###################################################################################################

cd ${project_folder}

for j in 0 14 28; do

  for i in $(eval echo "{$j..$(expr $j + 13)}"); do
  
  # VarScan Mutation calling

  java -Xmx8g -jar "$VarScan" \
    somatic \
    04_mpileup/"${mouse_id[$i]}".mpileup \
    05_varscan/"${mouse_id[$i]}".VarScan2.4.4 \
    --mpileup 1 \
    --min-coverage-normal 10 \
    --min-coverage-tumor 10 \
    --min-var-freq 0.07 \
    --output-vcf 1 &
    
    done
	
  wait
  
done

for i in {0..41}; do

  # Define Somatic vs LOH vs Germline variants

  for k in snp indel; do
  
    java -Xmx8g -jar "$VarScan" \
      processSomatic \
      05_varscan/"${mouse_id[$i]}".VarScan2.4.4."${k}".vcf \
      --min-tumor-freq 0.07 \
      --max-normal-freq 0.05 \
      --p-value 0.05
      
  done

  # Filter somatic variants around indels

  java -Xmx8g -jar "$VarScan" \
    somaticFilter \
    05_varscan/"${mouse_id[$i]}".VarScan2.4.4.snp.Somatic.hc.vcf \
    --min-var-freq 0.07 \
    --indel-file 05_varscan/"${mouse_id[$i]}".VarScan2.4.4.indel.vcf \
    --output-file 05_varscan/"${mouse_id[$i]}".VarScan2.4.4.snp.Somatic.hc.filter.vcf
    
done

###################################################################################################
# ANNOVAR ANNOTATION
###################################################################################################

for i in {0..41}; do

  convert2annovar.pl \
    --format vcf4old \
    -includeinfo \
    "$project_folder"/05_varscan/"${mouse_id[$i]}".VarScan2.4.4.snp_indel.Somatic.hc.filter.vcf \
    > "$project_folder"/06_annovar/"${mouse_id[$i]}".VarScan2.4.4.snp_indel.Somatic.hc.filter.vcf.annovar
  
  table_annovar.pl \
    "$project_folder"/06_annovar/"${mouse_id[$i]}".VarScan2.4.4.snp_indel.Somatic.hc.filter.vcf.annovar \
    "$annovar_db" \
    -buildver mm10 \
    -out "$project_folder"/06_annovar/"${mouse_id[$i]}".VarScan2.4.4.snp_indel.Somatic.hc.filter.vcf.annovar.out \
    -protocol refGene,snp142 \
    -operation g,f \
    -otherinfo \
    -remove

done
