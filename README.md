# WGS-analysis

---
title: "WGS Round 2 April 2022"
author: "Troy LoBue"
date: "4/1/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
This document can be found in the following directory along with other related documents:

**/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs**

This document includes all scripts and processing steps used for the second round of WGS data received in March 2022. This data consists of 40 WGS from the AHCC. 

**Looking at the intial QC provided by HudsonAlpha, the file BC-CR-116-1 has a significantly lower "% of Perfect Index Reads" and significantly higher "% One Mismatch Reads".** Everything else looks good. This can be seen in: **/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/HudsonAlpha_Files/Raw_Data/Demultiplex_Stats_H7NWJDSX3_20220323_081137 (2)**

All of the computational processing done here is being performed on the Easley cluster in the following directory: 

**/hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/**

All job submission numbers and walltimes are referring to the jobs ran for Group1. 

This pipeline follows the GATK version 4.1.9.0 [data preprocessing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) and [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) workflows. 

## Step1: Downloading Fastq Files from HudsonAlpha: 1_HudsonAlpha_File_transfer.sh
This script utilizs the MacOS curl script provided by HudsonAlpha to download all Fastq files to Easley.
The curl script can be downloaded from the HudsonAlpha website  [here](https://gslweb.discoveryls.com/information/software/wget_curl_download). The HudsonAlpha files with the download links are found in the directory:

**/WGS2_Apr2022_TL/HudsonAlpha_downloaded_FILES/FastQ_download_FILES/**

#### 1_HudsonAlpha_File_transfer.sh
```
#!/bin/sh

## these do not work for the WGS2 dataset
#wget -i files_H7NWJDSX3_20220323_081137.txt
#wget -i files_H7NWJDSX3_20220323_081437.txt

## I ran the script twice, once for each file
#bash hadiscovery_gsl_curl_download.sh files_H7NWJDSX3_20220323_081437.txt
bash hadiscovery_gsl_curl_download.sh files_H7NWJDSX3_20220323_081137.txt
```

## Step2: Renaming Fastq Files: 2_Forward_File_Renaming.sh & 2_Reverse_File_Renaming.sh
These scripts change the names of the fastq files to reflect our nomenclature system. I used the following file to get the nomenclature key:

**/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/HudsonAlpha_Files/Experimental_records/WGS-Batch-0047-2022-02-28.xlsx/**

I generated the following excel file to construct the name-change commands for each file:

**/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/2_fastq_Renaming_WGS2_Apr2022.xlsx/**  

This excel file was copied over to Easley to create two separate scripts, one for the forward files and one for the reverse files. Scripts below only show the first lines of the script. 


#### 2_Forward_File_Renaming.sh
```
#!/bin/sh

cp H7NWJDSX3_s1_1_IDT8_UDI_008_i7-IDT8_UDI_008_i5_6887-NDM-0001.fastq.gz BC-CR-91-1_1_Apr2022.fastq.gz
cp H7NWJDSX3_s1_1_IDT8_UDI_009_i7-IDT8_UDI_009_i5_6887-NDM-0002.fastq.gz BC-CR-93-1_1_Apr2022.fastq.gz
cp H7NWJDSX3_s1_1_IDT8_UDI_010_i7-IDT8_UDI_010_i5_6887-NDM-0003.fastq.gz BC-CR-94-1_1_Apr2022.fastq.gz
cp H7NWJDSX3_s1_1_IDT8_UDI_011_i7-IDT8_UDI_011_i5_6887-NDM-0004.fastq.gz BC-CR-95-1_1_Apr2022.fastq.gz
cp H7NWJDSX3_s1_1_IDT8_UDI_012_i7-IDT8_UDI_012_i5_6887-NDM-0005.fastq.gz BC-CR-98-1_1_Apr2022.fastq.gz
```

#### 2_Reverse_File_Renaming.sh
```
#!/bin/sh

cp H7NWJDSX3_s1_2_IDT8_UDI_008_i7-IDT8_UDI_008_i5_6887-NDM-0001.fastq.gz BC-CR-91-1_2_Apr2022.fastq.gz
cp H7NWJDSX3_s1_2_IDT8_UDI_009_i7-IDT8_UDI_009_i5_6887-NDM-0002.fastq.gz BC-CR-93-1_2_Apr2022.fastq.gz
cp H7NWJDSX3_s1_2_IDT8_UDI_010_i7-IDT8_UDI_010_i5_6887-NDM-0003.fastq.gz BC-CR-94-1_2_Apr2022.fastq.gz
cp H7NWJDSX3_s1_2_IDT8_UDI_011_i7-IDT8_UDI_011_i5_6887-NDM-0004.fastq.gz BC-CR-95-1_2_Apr2022.fastq.gz
cp H7NWJDSX3_s1_2_IDT8_UDI_012_i7-IDT8_UDI_012_i5_6887-NDM-0005.fastq.gz BC-CR-98-1_2_Apr2022.fastq.gz
```


## Step3: FastQC on Untrimmed files: 3_Untrimmed_FastQC.sh 
This step serves to understand the quality of the raw, untrimmed sequences. From this output we can decide if we want to implement a head crop or tail crop. The fastQC outputs were transferred from easley to:

**/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/FastQC_Untrimmed/**

These fastQC outputs were used to run MultiQC on [Galaxy](https://usegalaxy.org/). MultiQC is not available on Easley which is why Galaxy was used. The MultiQC output is also saved in the same folder on Box Drive. Look to Step 5 to see more details on how MultiQC was used. The multiQC output can be found at **/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/Fastqc_untrimmed/Untrimmed_MultiQC_WGS2.html**

#### 3_Untrimmed_FastQC.sh
Job submission 169978: Time elapsed = 31:00:00
```
#!/bin/sh

##Change to the working directory if you are working on Easley
cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL

#load FastQC on easley
module load fastqc

########### Run FastQC on original forward fastq files.
########### What is FastQC analysis: A simple program used to assess the quality of the sequencing data generated from high throughput sequencing pipelines.

### This step should be performed before and after trimming.
#run fastqc on original forward fastq files

for UNTRIMMED_F_FASTQ in *_1_*.fastq.gz
do

fastqc $UNTRIMMED_F_FASTQ

done

mkdir Fastqc_untrimmed
mv *.html Fastqc_untrimmed
tar -czvf Fastqc_untrimmed.tar Fastqc_untrimmed
```


## Step4: Trim fastq  files: 4_Trimmomatic_WGS2_Apr2022.sh
After analyzing the fastQC files it was determined that the files had high quality reads with no need for a head or tail crop. One file, BC-CR-116-1 (same file with the low % perfect index reads), has a high adapter content which should be fixed by the Trim. Following this job unpaired files were moved to the /Unpaired_trimmed_FILES/ directory. The trim logs were moved to the /Trim_logs_FILES/ directory. Finally, the untrimmed fastq files were moved to the /Raw_renamed_fastq_FILES/. More information on Trimmomatic can be found [here](http://www.usadellab.org/cms/?page=trimmomatic).


#### 4_Trimmomatic_WGS2_Apr2022.sh
Job submission 170995: Time elapsed = 210:00:00
```
#!/bin/sh

cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL
## load the module

#load the module on Easley
module load trimmomatic/0.39

ls *_1_*.fastq.gz > FSamplesList.txt

FILELIST=`cat FSamplesList.txt`

for FILENAME in $FILELIST
do
#BC-CR-100-1_1_Apr2022.fastq.gz
SHORTER=`echo $FILENAME | awk -F "." '{print $1}'`
SHORT=`echo $SHORTER | awk -F "_" '{print $1}'`

#Make sure that the path to the trimmomatic.jar file is correct
java -jar /tools/trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 48 -phred33 -trimlog $SHORTER.trim.log "$SHORT"_1_Apr2022.fastq.gz "$SHORT"_2_Apr2022.fastq.gz -baseout $SHORT.trim.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:10

done


#mv *trim_*U.fastq.gz

#gunzip "$FILENAME"_trim_1P.fastq.gz
#gunzip "$FILENAME"_trim_2P.fastq.gz
```


## Step5: FastQC on Trimmed files: 5_Trimmed_FastQC.sh 
FastQC is now run on the paired-end outputs of the Trimmomatic job. These files have the ...1P.fastq.gz ending. The Unpaired outputs and trim logs were moved to the **Files_to_save/4_Trimmomatic_WGS2_Apr2022_FILES** directory. Paired end trimmed files were moved to the same directory. **there was node preemtion during the the gunzip process of the trimmmed files so make sure you check the sizes of these files before you use them again.**  MultiQC was also used again for the fastQC outputs. MultiQC uses the "fastqc_data.txt" file found in the zip output of FastQC. A script for extracting and renaming these files can be found in: 

 **MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/Fastqc_untrimmed/FastQC_zip_FILES** 

Also the script is printed below the Trimmomatic script below. It is titled "File_extract.sh". The MultiQC outputs are in **/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/Fastqc_trimmed/Trimmed_MultiQC_WGS2.html**

More information on FastQC can be found [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

#### 5_Trimmed_FastQC.sh
Job submission 176000: Time elapsed = 31:00:00
```
#!/bin/sh

##Change to the working directory if you are working on Easley
cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL

#load FastQC on easley
module load fastqc

########### Run FastQC on original forward fastq files.
########### What is FastQC analysis: A simple program used to assess the quality of the sequencing data generated from high throughput sequencing pipelines.

### This step should be performed before and after trimming.
#run fastqc on original forward fastq files

for TRIMMED_F_FASTQ in *1P.fastq.gz
do

fastqc $TRIMMED_F_FASTQ

done

mkdir Fastqc_trimmed
mv *.html Fastqc_trimmed
tar -czvf Fastqc_trimmed.tar Fastqc_trimmed
```

#### File_extract.sh
```
#!/bin/sh


ls -d *_Apr2022_fastqc > FSamplesList.txt
FILELIST=`cat FSamplesList.txt`

for i in $FILELIST
do
#BC-CR-100-1_1_Apr2022_fastqc
#SHORTER=`echo $i | awk -F "." '{print $1}'`
SHORT=`echo $i | awk -F "_" '{print $1}'`

cd $i
mv fastqc_data.txt "$SHORT"_fastqc_data.txt
cp "$SHORT"_fastqc_data.txt /Users/tml0023/Library/CloudStorage/Box-Box/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/Fastqc_untrimmed/FastQC_zip_FILES/Raw_fastQC_data_FILES
cd /Users/tml0023/Library/CloudStorage/Box-Box/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/Fastqc_untrimmed/FastQC_zip_FILES

done
```


## Step 6: Alignment to reference genome with Burrow-Wheeler Aligner: 6_BWA_hg38_WGS2_Apr2022.sh
Before performing this step I have split all of the files into four groups to speed up the following steps. There are now four directories with 10 files in each, which can all be found in the **/WGS2_Apr2022_TL** directory. Keep this in mind when looking at the 4 scripts because there will be minor differences in the scripts for each group. I have tried to point out where these changes are in the sample script below. 

This step aligns the sequencing data to the hg38 reference genome which can be found in:

**/hosted/cvmpt/archive/Human_Genome/genome** 

More information on BWA can be found [here](http://bio-bwa.sourceforge.net/bwa.shtml).
#### 6_BWA_hg38_WGS2_Apr2022_Group1.sh
Job submision 176678: Time elapsed = 120:00:00
```
#!/bin/sh
##Needs to be ran on fastq files that have already been trimmed
##Make sure to adjust the number of the threads requested to match the number of available processors per core on Easley

cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group<INSERT GROUP #>

############## Module Load #####################################
# bwa/0.7.17
# samtools/1.11

module load bwa
module load samtools

######### For Loop to align using BWA mem files #############################
##This requires a preexisting file containing in list form all of the basic file names used;
## only needs to have the first portion namesd to match - variables will be set below to 
##isolate portion of name needed for each step


FILELIST=`cat FSamplesList_Group<INSERT GROUP #>.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`

#Output = BC-CR-100-1
        #########################################################
### Aligning with 48 threads/cores/processors
##For reference Genome specify the directory and then the base part of the genome name
#####The part before the .fa

gzip -d "$SHORT".trim_1P.fastq.gz "$SHORT".trim_2P.fastq.gz
        
        bwa mem -t 48 -M /hosted/cvmpt/archive/Human_Genome/genome "$SHORT".trim_1P.fastq \
        "$SHORT".trim_2P.fastq > "$SHORT".mem.sam

        ##################################################################################

samtools view -Sb -@ 48 "$SHORT".mem.sam -o "$SHORT".mem.bam
samtools sort -@ 48 "$SHORT".mem.bam -o "$SHORT".memsorted.bam
samtools index -@ 48 "$SHORT".memsorted.bam

gzip "$SHORT".trim_1P.fastq "$SHORT".trim_2P.fastq

done
```

## Step 6a: Extract unmapped reads: 6a_Extract_Unmapped_reads.sh
There are projects going on in the lab working on identifying oncoviruses in the unmapped reads following Burrows-Wheeler Alignment. Because of this we have started extracting unmapped reads from the inital alignment (BAM) files. This is the script for carrying out this process. These unmapped files were then moved to:

**/hosted/cvmpt/archive/Files_to_Share/May2022_WGS2_UnmappedReads**

#### 6a_Extract_Unmapped_reads.sh
Job submission 184741: Time elapsed = 4:00:00
```
#!/bin/sh

cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group1

############## Module Load #####################################

# samtools/1.11

module load samtools

######### For Loop to align using BWA mem files #############################
##This requires a preexisting file containing in list form all of the basic file names used;
## only needs to have the first portion namesd to match - variables will be set below to 
##isolate portion of name needed for each step


FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`

samtools view -b -f 4 "$SHORT".memsorted.bam > "$SHORT".unmapped.bam

done
```

## Step 6b: Add or Replace Read Groups: 6b_AddReadGroups_WGS2_Apr2022_Group1.sh
In this step we are adding new read groups to each read in the bam files. This does not change the actual read data it actually just changes the metadata for each read. You must have read groups assigned in order to carry on through the GATK pipeline. I had initially forgot to sort and index this output, but I have fixed that on this document.

#### 6b_AddReadGroups_WGS2_Apr2022_Group1.sh
Job submision 184669: Time elapsed = 16:00:00
``` 
#!/bin/sh

module load samtools
module load gatk
module load picard

FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`

#Output = BC-CR-100-1

 java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
       I="$SHORT".memsorted.bam \
       O="$SHORT".rg.bam \
       RGID="$SHORT" \
       RGLB=WGS2 \
       RGPL=ILLUMINA \
       RGPU=<sequencer> \
       RGSM="$SHORT"_WGS2
 
## Now need to sort and index with samtools
samtools sort -@ 48 "$SHORT".rg.bam -o "$SHORT".rgsorted.bam
samtools index -@ 48 "$SHORT".rgsorted.bam

done
```

## Step 7: Marking Duplicate reads: 7_MarkedDuplicates_WGS2_Apr2022.sh
In this step we will be identifying and marking duplicate reads. This generates new bam files whitch includes this informaation. More information can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360051306171-MarkDuplicates-Picard-). 

#### 7_MarkedDuplicates_WGS2_Apr2022.sh
Job submission 199140: Time elapsed = 39:00:00
```
#!/bin/sh

module load picard
module load samtools

FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

##SHORT = BC-CR-100-1
#BC-CR-100-1.memsorted.bam 

java -jar /tools/picard-2.23.9/libs/picard.jar MarkDuplicates \
      I="$SHORT".rgsorted.bam \
      O="$SHORT".markdup.bam \
      M="$SHORT".marked_dup_metrics.txt

samtools sort -@ 48 "$SHORT".markdup.bam -o "$SHORT".markdup.sorted.bam
samtools index -@ 48 "$SHORT".markdup.sorted.bam

done
```

## Step7a: Analyzing Depth of Coverage: 7a_DepthOfCoverage_WGS2_Apr2022.sh
On this step we are analyzing the quality of the covered regions by determining the depth. This is a quality step that adds credibility to the files generated. This job will be performed before and after BQSR to compare the effectiveness. Documentation for this job can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360051307491-DepthOfCoverage-BETA-).

*wgs_calling_regions.hg38.interval_list* is a list of intervals provided by GATK that determine which regions of the genome to blacklist. These are centromeric regions and large regions with no genes. More information can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists). This file will be used again in GenomicDBimport. 

#### 7a_DepthOfCoverage_WGS2_Apr2022.sh
Job submission 206019: Time elapse = 95:00:00
```
#!/bin/sh

module load gatk 
module load samtools

FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

##SHORT = BC-CR-100-1
#BC-CR-100-1.memsorted.bam 

 gatk DepthOfCoverage \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -O "$SHORT".DOC_base \
   -I "$SHORT".markdup.sorted.bam \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   --summary-coverage-threshold 1 --summary-coverage-threshold 4 --summary-coverage-threshold 6 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --summary-coverage-threshold 30 --summary-coverage-threshold 35 --summary-coverage-threshold 40 --summary-coverage-threshold 45 --summary-coverage-threshold 50 \

done 
```

## Step7b: Use flagstat for more quality control: 7b_Flagstat_WGS2_Apr2022.sh 
This is another quality step that will be conducted before and after BQSR. I included the final output image below. FYI you must use a .png file if you want to include images in Markdown. The code used to create this image can be found in **WGS2_QC_plots.R**. More information on Flagstat can be found [here](http://www.htslib.org/doc/samtools-flagstat.html). 

#### 7b_Flagstat_WGS2_Apr2022.sh 
Job submission 206157: Time elapse = 3:00:00
```
#!/bin/sh

module load samtools

FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

##SHORT = BC-CR-100-1
#BC-CR-100-1.memsorted.bam 

samtools flagstat "$SHORT".markdup.sorted.bam > "$SHORT".flagstat.txt

done
```


```{r, echo=FALSE, out.width = '100%'}
knitr::include_graphics("~/Library/CloudStorage/Box-Box/MernerLab_General/Research_Projects/WGS1_Human_DataSets/WGS_2_AAs/Images/Flagstat_plot.png")
```


## Step 8: Base Recalibration, Application, and Plot Production: 8_BaseRecalibrator_WGS2_Apr2022.sh
This step uses variant databases such as dbSNP, gnomad, 1000 genome, HapMap, and Indel databases to recalibrate the base quality for the bam files. I downloaded updated dbSNP (version 155) and Indel databases from [NCBI](https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/) and [UCSC](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false), respectively. Unfortunately, these updated databases did not work as there was inconsistencies with contig naming adn the hg38 reference genome. So I went back to using older versions of of these databases, dbSNP_146 and Mills_and_1000G_gold_standard.indels.hg38. This script incorporates multiple jobs and is a little more complex. You can find information about these at [BaseRacalibration](https://gatk.broadinstitute.org/hc/en-us/articles/360050815072-BaseRecalibrator), [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR), and [Analyze Covariates](https://gatk.broadinstitute.org/hc/en-us/articles/360051304351-AnalyzeCovariates). Additional information regarding the databases can be found [here](https://gatk.broadinstitute.org/hc/en-us/community/posts/360075305092-Known-Sites-for-BQSR). 

#### 8_BaseRecalibrator_WGS2_Apr2022.sh
Job submission 205707: Time elapse = 72:00:00
```
#!/bin/sh

module load gatk 
module load samtools
module load R

FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

##SHORT = BC-CR-100-1
#BC-CR-100-1.memsorted.bam 

 gatk BaseRecalibrator \
   -I "$SHORT".rgsorted.bam \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   --known-sites /hosted/cvmpt/Human_Research/KnownSites/dbSNP150.hg38.vcf \
   --known-sites /hosted/cvmpt/Human_Research/KnownSites/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --known-sites /hosted/cvmpt/Human_Research/KnownSites/Homo_sapiens_assembly38.known_indels.vcf \
   --known-sites /hosted/cvmpt/Human_Research/KnownSites/hapmap_3.3.hg38.vcf \
   -O "$SHORT".recal_data.table
 
 gatk ApplyBQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -I "$SHORT".markdup.sorted.bam \
   --bqsr-recal-file "$SHORT".recal_data.table \
   -O "$SHORT".recal.bam
 
 gatk AnalyzeCovariates \
   -bqsr "$SHORT".recal_data.table \
   -plots "$SHORT".AnalyzeCovariates.pdf

samtools sort -@ 48 "$SHORT".recal.bam -o "$SHORT".recal.sorted.bam
samtools index -@ 48 "$SHORT".recal.sorted.bam

done
```

## Step 8a: Compiled Depth of Coverage for Recal bam Files: 8a_DepthOfCoverage_WGS2_Apr2022.sh
This depth of coverage is slightly different than the one I ran for the Marked Duplicate files. Here I used a file list as an input so that all files will be included in one output, rather than a separate output for each file. This made it easier to generate plots for the coverage data. The code for these plots can be found in **WGS2_QC_plots.R**. 

#### 8a_DepthOfCoverage_WGS2_Apr2022.sh
Job submission 230524: Time elapse = 110:00:00
```
#!/bin/sh

module load gatk 
module load samtools

FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

##SHORT = BC-CR-100-1
#BC-CR-100-1.memsorted.bam 

 gatk DepthOfCoverage \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -O Group1.DOC_base \
   -I Group1_bams.list \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   --summary-coverage-threshold 0 --summary-coverage-threshold 4 --summary-coverage-threshold 6 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --summary-coverage-threshold 30 --summary-coverage-threshold 35 --summary-coverage-threshold 40 --summary-coverage-threshold 45 --summary-coverage-threshold 50 \

done
```

## Step 9: Calling variants & generating VCF's: 9_HaplotypeCaller_WGS2_Apr2022.sh
Here we are finally calling variants with HaplotypeCaller which gives us our VCF files for each sample. These include the all variants for each sample but with none of the annotations. This is one of the longest steps so make sure to allocate plenty of time and memory. More information can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360050814612-HaplotypeCaller#--annotation-group).

#### 9_HaplotypeCaller_WGS2_Apr2022.sh
Job submission 223187: Time elapse = 250:00:00
```
#!/bin/sh

module load gatk 

FILELIST=`cat FSamplesList_Group1.txt`  ##Can be used if a file list is needed
for FILE in $FILELIST; do

        #########Creating Variable for File names#################
##Need to be changed based on naming system
#Input = BC-CR-100-1_1_Apr2022.fastq.gz

#SHORTER=`echo $FILE | awk -F "." '{print $1}'`
SHORT=`echo $FILE | awk -F "_" '{print $1}'`     

##SHORT = BC-CR-100-1
#BC-CR-100-1.memsorted.bam 

 gatk --java-options "-Xmx80g" HaplotypeCaller  \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -I "$SHORT".recal.sorted.bam \
   -O "$SHORT".g.vcf.gz \
   -ERC GVCF \
   -A AlleleFraction \
   -A BaseQuality \
   -A MappingQuality \
   --native-pair-hmm-threads 48 
   

done
```

#### Rename_VCF_sample.sh
This script is used to rename the sample name for each sample. This allows the files to be used in RVtests.

```
#!/bin/sh

module load picard

FILELIST=`cat OldSampleNames.txt`
for FILE in $FILELIST

do 
## input = BC-CR-100-1.g.vcf.gz

SHORT=`echo $FILE | awk -F "." '{print $1}'`

java -jar /tools/picard-2.23.9/libs/picard.jar RenameSampleInVcf \
      INPUT="$SHORT".g.vcf.gz \
      OUTPUT="$SHORT".rename.g.vcf.gz \
      NEW_SAMPLE_NAME="$SHORT"

done
```

#### Index_renamed.sh
Here we index the new renamed vcf files that were just generated for RVtests.

```
#!/bin/sh

module load gatk

FILELIST=`cat OldSampleNames.txt`
for FILE in $FILELIST

do 
## input = BC-CR-100-1.g.vcf.gz

SHORT=`echo $FILE | awk -F "." '{print $1}'`

 gatk IndexFeatureFile \
     -I "$SHORT".rename.g.vcf.gz

done
```

## Step 10: Merging VCf files: 10_GenomicsDBimport_WGS2_Apr2022.sh
After HaplotypeCaller we need to merge all the generated VCF files into a single workspace. This is done with GenomicsDBimport. All VCF files generated by HaplotypeCaller were moved to a shared directory named:

*/hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group_ALL*

Here I run Genomics Dbimport and generate one compiled workspace */Group_ALL/WGS2_GDBI_Workspace/*. This workspace is used in the next step to generate the compiled VCF file for all samples. More information on GenomicsDBimport can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360051305591-GenomicsDBImport).

**tmp-dir needs to be created prior to submitting the script** 

#### 10_GenomicsDBimport_WGS2_Apr2022.sh
Job submission 245790: Time elapse = 75:00:00
```
#!/bin/sh
module load samtools
module load picard
module load gatk
cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group_ALL
gatk --java-options "-Xmx550g -Xms550g" GenomicsDBImport \
-V BC-CR-100-1.g.vcf.gz \
-V BC-CR-114-1.g.vcf.gz \
-V BC-CR-116-1.g.vcf.gz \
-V BC-CR-120-1.g.vcf.gz \
-V BC-CR-123-1.g.vcf.gz \
-V BC-CR-91-1.g.vcf.gz \
-V BC-CR-93-1.g.vcf.gz \
-V BC-CR-94-1.g.vcf.gz \
-V BC-CR-95-1.g.vcf.gz \
-V BC-CR-98-1.g.vcf.gz \
-V BC-CR-99-1.g.vcf.gz \
-V BC-EAMC-124-1.g.vcf.gz \
-V BC-EAMC-125-1.g.vcf.gz \
-V BC-EAMC-130-1.g.vcf.gz \
-V BC-EAMC-131-1.g.vcf.gz \
-V BC-EAMC-133-1.g.vcf.gz \
-V BC-EAMC-134-1.g.vcf.gz \
-V BC-EAMC-137-1.g.vcf.gz \
-V BC-EAMC-139-1.g.vcf.gz \
-V BC-EAMC-147-1.g.vcf.gz \
-V BC-EAMC-149-1.g.vcf.gz \
-V BC-EAMC-159-1.g.vcf.gz \
-V BC-EAMC-163-1.g.vcf.gz \
-V BC-EAMC-164-1.g.vcf.gz \
-V BC-EAMC-173-1.g.vcf.gz \
-V BC-EAMC-176-1.g.vcf.gz \
-V BC-EAMC-188-1.g.vcf.gz \
-V BC-EAMC-191-1.g.vcf.gz \
-V BC-EAMC-193-1.g.vcf.gz \
-V BC-EAMC-194-1.g.vcf.gz \
-V BC-EAMC-199-1.g.vcf.gz \
-V BC-EAMC-201-1.g.vcf.gz \
-V BC-EAMC-202-1.g.vcf.gz \
-V BC-EAMC-204-1.g.vcf.gz \
-V BC-EAMC-208-1.g.vcf.gz \
-V BC-EAMC-209-1.g.vcf.gz \
-V BC-EAMC-213-1.g.vcf.gz \
-V BC-EAMC-220-1.g.vcf.gz \
-V BC-EAMC-223-1.g.vcf.gz \
-V BC-EAMC-87-1.g.vcf.gz \
-L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
--genomicsdb-workspace-path ./WGS2_GDBI_Workspace \
--tmp-dir ./WGS2_tmp \
--reader-threads 48 \
```

## Step 10a: Merginig 20 samples from WGS1 with 40 samples from WGS2: 10a_GenomicsDBimport_WGSALL_Apr2022.sh
I made a GDBI workspace for only the WGS2 files called *WGS2_GDBI_Workspace*. I copied this directory into *WGSALL_GDBI_Workspace* which is where I added the 20 WGS1 samples to. So at this point *WGSALL_GDBI_Workspace* consists of WGS1 (20 samples) and WGS2 (40samples).

*WGSALL_GDBI_Workspace should also be where all additional samples are added in the future so that all samples will be held in a single consolidated workspace. This is important to remember and continue to keep updated.*

#### 10a_GenomicsDBimport_WGSALL_Apr2022.sh
Job submission 250145: Time elapse = 55:00:00
```
#!/bin/sh
module load samtools
module load picard
module load gatk
cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group_ALL
gatk --java-options "-Xmx550g -Xms550g" GenomicsDBImport \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-14-6.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-17-2.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-23-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-28-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-32-2.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-33-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-36-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-38-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-4-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-50-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-51-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-CR-60-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-17-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-18-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-2-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-30-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-32-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-35-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-57-1.erc.g.vcf \
-V /hosted/cvmpt/archive/WGS_Human/BAM_files/vcf_files/BC-EAMC-86-1.erc.g.vcf \
-L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
--genomicsdb-update-workspace-path ./WGSALL_GDBI_Workspace \
--tmp-dir ./WGSALL_tmp \
--reader-threads 48 \
```

## Step 11: Generating single VCF: 11_GenotypeGVCF_WGS2_Apr2022.sh
This step utilizes the the GDBI workspace generated in the last step to create a single joint genotyped VCF file that incorporates all samples. I will do this for only WGS2 samples and then again using the *WGSALL_GDBI_Workspace* for all 60 samples. More info [here](https://gatk.broadinstitute.org/hc/en-us/articles/360050816072-GenotypeGVCFs).

#### 11_GenotypeGVCF_WGS2_Apr2022.sh
Job submission 248516: Time elapse 96:00:00
```
#!/bin/sh

module load gatk 

 gatk --java-options "-Xmx550g" GenotypeGVCFs \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V gendb://WGS2_GDBI_Workspace \
   -O GenotypeGVCF_output_WGS2.vcf.gz \
   --tmp-dir GenotypeGVCF_tmp
```

## Step 11a: Generating single VCF for WGS1 and WGS2 samples: 11a_GenotypeGVCF_WGSALL_Apr2022.sh
In this step I generate the master VCF for WGS1 samples and WGS2 samples. 

#### 11a_GenotypeGVCF_WGSALL_Apr2022.sh
Job submission 252876: Time elapse 200:00:00
```
#!/bin/sh

module load gatk 

 gatk --java-options "-Xmx550g" GenotypeGVCFs \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V gendb://WGSALL_GDBI_Workspace \
   -O GenotypeGVCF_output_WGSALL.vcf.gz \
   --tmp-dir GenotypeGVCF_ALL_tmp
```

## Step 12: Recalibrating variants: 12_VQSR_WGS2_Apr2022.sh
This step uses two jobs from GATK, [VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360050815872-VariantRecalibrator#--resource) and [ApplyVQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360051306591-ApplyVQSR). This script assigns values to the variants that acts as a quality metric. Filters can be applied in this step to control what variants are present in the final VCF file. More information can be found in the links above. 

#### 12_VQSR_WGS2_Apr2022.sh
Job submission 250362: Time elapse = 02:00:00
```
#!/bin/sh

module load gatk
module load R

## recalibrate SNPs
 gatk VariantRecalibrator \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V GenotypeGVCF_output_WGS2.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /hosted/cvmpt/Human_Research/KnownSites/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /hosted/cvmpt/Human_Research/KnownSites/dbSNP150.hg38.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /hosted/cvmpt/Human_Research/KnownSites/hapmap_3.3.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
   -mode SNP \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   -O VQSR_SNP_WGS2_Output.vcf.recal \
   --tranches-file VQSR_SNP_WGS2_Output.vcf.tranches \
   --rscript-file VQSR_SNP_WGS2_Output.vcf.plots.R

 gatk ApplyVQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V GenotypeGVCF_output_WGS2.vcf.gz \
   -O ApplyVQSR_SNPs_WGS2_Output.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file VQSR_SNP_WGS2_Output.vcf.tranches \
   --recal-file VQSR_SNP_WGS2_Output.vcf.recal \
   -mode SNP

## reaclibrate for INDELs
 gatk VariantRecalibrator \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V ApplyVQSR_SNPs_WGS2_Output.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /hosted/cvmpt/Human_Research/KnownSites/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /hosted/cvmpt/Human_Research/KnownSites/dbSNP150.hg38.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /hosted/cvmpt/Human_Research/KnownSites/hapmap_3.3.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   -O ApplyVQSR_INDEL_SNP_WGS2_Output.vcf.recal \
   --tranches-file VQSR_INDEL_SNP_WGS2_Output.vcf.tranches \
   --rscript-file VQSR_INDEL_SNP_WGS2_Output.vcf.plots.R

 gatk ApplyVQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V ApplyVQSR_SNPs_WGS2_Output.vcf.gz \
   -O ApplyVQSR_INDEL_SNP_WGS2_Output.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file VQSR_INDEL_SNP_WGS2_Output.vcf.tranches \
   --recal-file ApplyVQSR_INDEL_SNP_WGS2_Output.vcf.recal \
   -mode INDEL
```

## Step 12a: Recalibration for all 60 samples: 12a_VQSR_WGSALL_Apr2022.sh

#### 12a_VQSR_WGSALL_Apr2022.sh
Job submission 263018: Time elapse = 02:00:00
```
#!/bin/sh

module load gatk
module load R

## recalibrate SNPs
 gatk VariantRecalibrator \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V GenotypeGVCF_output_WGSALL.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /hosted/cvmpt/Human_Research/KnownSites/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /hosted/cvmpt/Human_Research/KnownSites/dbSNP150.hg38.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /hosted/cvmpt/Human_Research/KnownSites/hapmap_3.3.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
   -mode SNP \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   -O VQSR_SNP_WGSALL_Output.vcf.recal \
   --tranches-file VQSR_SNP_WGSALL_Output.vcf.tranches \
   --rscript-file VQSR_SNP_WGSALL_Output.vcf.plots.R

 gatk ApplyVQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V GenotypeGVCF_output_WGSALL.vcf.gz \
   -O ApplyVQSR_SNPs_WGSALL_Output.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file VQSR_SNP_WGSALL_Output.vcf.tranches \
   --recal-file VQSR_SNP_WGSALL_Output.vcf.recal \
   -mode SNP

## reaclibrate for INDELs
 gatk VariantRecalibrator \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V ApplyVQSR_SNPs_WGSALL_Output.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /hosted/cvmpt/Human_Research/KnownSites/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /hosted/cvmpt/Human_Research/KnownSites/dbSNP150.hg38.vcf \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /hosted/cvmpt/Human_Research/KnownSites/hapmap_3.3.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   -L /hosted/cvmpt/Human_Research/KnownSites/wgs_calling_regions.hg38.interval_list \
   -O ApplyVQSR_INDEL_SNP_WGSALL_Output.vcf.recal \
   --tranches-file VQSR_INDEL_SNP_WGSALL_Output.vcf.tranches \
   --rscript-file VQSR_INDEL_SNP_WGSALL_Output.vcf.plots.R

 gatk ApplyVQSR \
   -R /hosted/cvmpt/archive/Human_Genome/genome.fa \
   -V ApplyVQSR_SNPs_WGSALL_Output.vcf.gz \
   -O ApplyVQSR_INDEL_SNP_WGSALL_Output.vcf.gz \
   --truth-sensitivity-filter-level 90.0 \
   --tranches-file VQSR_INDEL_SNP_WGSALL_Output.vcf.tranches \
   --recal-file ApplyVQSR_INDEL_SNP_WGSALL_Output.vcf.recal \
   -mode INDEL

```

## Step 13: Annotate variants with ANNOVAR: 13_Annovar_WGS2_Apr2022.sh
This  step annotates all the variants with allele frequency, predicted mutation consequences, and other useful information. ANNOVAR is an open source program with more documentation found [here](https://annovar.openbioinformatics.org/en/latest/).

#### 13_Annovar_WGS2_Apr2022.sh
Job submission 250387: Time elapse = 11:00:00
```
#!/bin/sh

cd /hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group_ALL 

module load perl
module load annovar

table_annovar.pl ApplyVQSR_INDEL_SNP_WGS2_Output.vcf.gz /hosted/cvmpt/archive/Human_Genome/humandb \
-buildver hg38 \
-protocol refGene,dbnsfp42a,dbscsnv11,revel,avsnp150,clinvar_20220320,gnomad30_genome,intervar_20180118 \
--operation g,f,f,f,f,f,f,f --vcfinput \
--argument '--hgvs --exonicsplicing --splicing_threshold 5',,,,,,, \
--argument '--hgvs --indel_splicing_threshold 100',,,,,,, \
--outfile Annotated_WGS2_40_MASTER_Apr2022.vcf.gz
```

## Step 13a: Annotate all 60 samples with ANNOVAR: 13a_Annovar_WGSALL_Apr2022.sh 

#### 13a_Annovar_WGSALL_Apr2022.sh 
Job submission 263094: Time elapse = 20:00:00
```
#!/bin/sh

cd /hosted/cvmpt/archive/WGS_Human/z_ROUND_2_Mar2022_TL/Group_ALL 

module load perl
module load annovar

table_annovar.pl ApplyVQSR_INDEL_SNP_WGSALL_Output.vcf.gz /hosted/cvmpt/archive/Human_Genome/humandb \
-buildver hg38 \
-protocol refGene,dbnsfp42a,dbscsnv11,revel,avsnp150,clinvar_20220320,gnomad30_genome,intervar_20180118 \
--operation g,f,f,f,f,f,f,f --vcfinput \
--argument '--hgvs --exonicsplicing --splicing_threshold 5',,,,,,, \
--argument '--hgvs --indel_splicing_threshold 100',,,,,,, \
--outfile Annotated_WGSALL_60_MASTER_Apr2022.vcf.gz

```





## Important file locations

#### Raw fastq files
**/hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Files_to_save/2_File_renaming_FILES**

#### BAM files
Note that bam files will be found in the specific Group directories (Group1-4).

**/hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group1/Files_To_Save/6_BWA_hg38_WGS2_Apr2022_Group1_FILES**

#### Unmapped BAM files

**/hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group1/Files_To_Save/6a_Extract_Unmapped_reads_FILES**

#### GenomicsDB Workspace
The GenomicDB workspace that should be continued to be added to with all future sequencing data is located here: 

**/hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group_ALL/Files_to_save/10a_GenomicsDBimport_WGSALL_Apr2022_FILES**

This directory can be copied over to new directories that are used to process future data. 

#### MASTER VCF files
The MASTER vcf files for currently all 60 WGS samples is located here: 

**/hosted/cvmpt/archive/WGS_Human/WGS2_Apr2022_TL/Group_ALL/Annotated_MASTER_WGSALL_60_FILES**
