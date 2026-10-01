# WGS Germline Variant Calling Pipeline

A GATK-based pipeline for processing whole-genome sequencing (WGS) data from raw FASTQ files through joint-genotyped, annotated variant calls. Built around GATK 4.1.9.0 best practices for [data pre-processing](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) and [germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-).

Designed to run on an HPC cluster using environment modules (`module load ...`) and SLURM-style batch job submission, but the steps translate directly to any Linux cluster or workstation with the required tools installed.

## Pipeline Overview

1. [Download FASTQ files](#1-download-fastq-files)
2. [Rename FASTQ files](#2-rename-fastq-files)
3. [FastQC on raw reads](#3-fastqc-on-untrimmed-reads)
4. [Trim reads with Trimmomatic](#4-trim-reads)
5. [FastQC on trimmed reads](#5-fastqc-on-trimmed-reads)
6. [Align to reference genome (BWA-MEM)](#6-align-to-reference-genome)
   - 6a. [Extract unmapped reads](#6a-extract-unmapped-reads)
   - 6b. [Add/replace read groups](#6b-add-or-replace-read-groups)
7. [Mark duplicates](#7-mark-duplicate-reads)
   - 7a. [Depth of coverage (pre-BQSR)](#7a-depth-of-coverage)
   - 7b. [Flagstat QC](#7b-flagstat-qc)
8. [Base quality score recalibration (BQSR)](#8-base-quality-score-recalibration)
   - 8a. [Depth of coverage (post-BQSR)](#8a-depth-of-coverage-post-bqsr)
9. [Call variants per sample (HaplotypeCaller)](#9-call-variants-per-sample)
10. [Consolidate GVCFs (GenomicsDBImport)](#10-consolidate-gvcfs)
11. [Joint genotyping (GenotypeGVCFs)](#11-joint-genotyping)
12. [Variant quality score recalibration (VQSR)](#12-variant-quality-score-recalibration)
13. [Annotate variants (ANNOVAR)](#13-annotate-variants)

Each step below documents the tool used, what it does, and a generalized version of the script. Replace bracketed placeholders (`<...>`) with your own paths, sample names, and resource allocations.

---

## 1. Download FASTQ Files

Raw FASTQ files are retrieved from the sequencing provider using a vendor-supplied download script (e.g., a `curl`- or `wget`-based transfer utility).

```bash
#!/bin/sh

# Run once per manifest/batch file provided by the sequencing vendor
bash <vendor_download_script>.sh <manifest_file>.txt
```

## 2. Rename FASTQ Files

Raw FASTQ filenames from the sequencing provider are renamed to match an internal sample-naming convention (e.g., `<SAMPLE_ID>_1_<batch>.fastq.gz` / `<SAMPLE_ID>_2_<batch>.fastq.gz` for forward/reverse reads). The exact renaming logic depends on your own naming scheme and the provider's output format.

## 3. FastQC on Untrimmed Reads

Assesses raw read quality prior to trimming, to determine whether head/tail cropping is needed.

```bash
#!/bin/sh

module load fastqc

for FASTQ in *_1_*.fastq.gz; do
    fastqc "$FASTQ"
done

mkdir -p fastqc_untrimmed
mv *.html fastqc_untrimmed
tar -czvf fastqc_untrimmed.tar fastqc_untrimmed
```

FastQC outputs can be aggregated with [MultiQC](https://multiqc.info/) (e.g., via [Galaxy](https://usegalaxy.org/) if MultiQC isn't installed locally).

## 4. Trim Reads

Trims adapters and low-quality bases with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Adjust `ILLUMINACLIP`, `SLIDINGWINDOW`, `HEADCROP`, etc. based on your FastQC results.

```bash
#!/bin/sh

module load trimmomatic

ls *_1_*.fastq.gz > forward_samples.txt

for FASTQ in $(cat forward_samples.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')

    java -jar <path_to>/trimmomatic.jar PE -threads <N> -phred33 \
        -trimlog "$SAMPLE".trim.log \
        "$SAMPLE"_1_<batch>.fastq.gz "$SAMPLE"_2_<batch>.fastq.gz \
        -baseout "$SAMPLE".trim.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        HEADCROP:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:10
done
```

## 5. FastQC on Trimmed Reads

Re-run FastQC on the paired trimmed output (`*_1P.fastq.gz` / `*_2P.fastq.gz`) to confirm trimming improved read quality, then re-aggregate with MultiQC.

```bash
#!/bin/sh

module load fastqc

for FASTQ in *1P.fastq.gz; do
    fastqc "$FASTQ"
done

mkdir -p fastqc_trimmed
mv *.html fastqc_trimmed
tar -czvf fastqc_trimmed.tar fastqc_trimmed
```

A helper script can extract and rename each sample's `fastqc_data.txt` from its FastQC zip output for downstream MultiQC aggregation:

```bash
#!/bin/sh

ls -d *_fastqc > fastqc_dirs.txt

for DIR in $(cat fastqc_dirs.txt); do
    SAMPLE=$(echo "$DIR" | awk -F "_" '{print $1}')
    cd "$DIR"
    mv fastqc_data.txt "$SAMPLE"_fastqc_data.txt
    cp "$SAMPLE"_fastqc_data.txt <output_directory>/
    cd ..
done
```

## 6. Align to Reference Genome

Aligns trimmed reads to the reference genome (e.g., hg38) using [BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml), then sorts and indexes with samtools.

For large cohorts, splitting samples into groups and running them in parallel batches can significantly speed up this and subsequent steps.

```bash
#!/bin/sh

module load bwa
module load samtools

cd <working_directory>

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')

    gzip -d "$SAMPLE".trim_1P.fastq.gz "$SAMPLE".trim_2P.fastq.gz

    bwa mem -t <N> -M <reference_genome_prefix> \
        "$SAMPLE".trim_1P.fastq "$SAMPLE".trim_2P.fastq > "$SAMPLE".mem.sam

    samtools view -Sb -@ <N> "$SAMPLE".mem.sam -o "$SAMPLE".mem.bam
    samtools sort -@ <N> "$SAMPLE".mem.bam -o "$SAMPLE".sorted.bam
    samtools index -@ <N> "$SAMPLE".sorted.bam

    gzip "$SAMPLE".trim_1P.fastq "$SAMPLE".trim_2P.fastq
done
```

### 6a. Extract Unmapped Reads

Optionally extract unmapped reads from the alignment for downstream analyses (e.g., pathogen/oncovirus screening).

```bash
#!/bin/sh

module load samtools

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')
    samtools view -b -f 4 "$SAMPLE".sorted.bam > "$SAMPLE".unmapped.bam
done
```

### 6b. Add or Replace Read Groups

GATK requires read group metadata on every read. This step adds/replaces read groups without altering the underlying alignment data, then re-sorts and indexes.

```bash
#!/bin/sh

module load samtools
module load gatk
module load picard

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')

    java -jar <path_to>/picard.jar AddOrReplaceReadGroups \
        I="$SAMPLE".sorted.bam \
        O="$SAMPLE".rg.bam \
        RGID="$SAMPLE" \
        RGLB=<library> \
        RGPL=ILLUMINA \
        RGPU=<sequencer_id> \
        RGSM="$SAMPLE"

    samtools sort -@ <N> "$SAMPLE".rg.bam -o "$SAMPLE".rgsorted.bam
    samtools index -@ <N> "$SAMPLE".rgsorted.bam
done
```

## 7. Mark Duplicate Reads

Identifies and flags PCR/optical duplicate reads with [Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360051306171-MarkDuplicates-Picard-).

```bash
#!/bin/sh

module load picard
module load samtools

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')

    java -jar <path_to>/picard.jar MarkDuplicates \
        I="$SAMPLE".rgsorted.bam \
        O="$SAMPLE".markdup.bam \
        M="$SAMPLE".marked_dup_metrics.txt

    samtools sort -@ <N> "$SAMPLE".markdup.bam -o "$SAMPLE".markdup.sorted.bam
    samtools index -@ <N> "$SAMPLE".markdup.sorted.bam
done
```

### 7a. Depth of Coverage

Runs [DepthOfCoverage](https://gatk.broadinstitute.org/hc/en-us/articles/360051307491-DepthOfCoverage-BETA-) before BQSR as a baseline QC metric, restricted to a standard set of callable intervals (GATK's [`wgs_calling_regions`](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) interval list, which excludes centromeric/low-complexity regions).

```bash
#!/bin/sh

module load gatk
module load samtools

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')

    gatk DepthOfCoverage \
        -R <reference_genome>.fa \
        -O "$SAMPLE".depth_of_coverage \
        -I "$SAMPLE".markdup.sorted.bam \
        -L <wgs_calling_regions>.interval_list \
        --summary-coverage-threshold 1 --summary-coverage-threshold 4 \
        --summary-coverage-threshold 6 --summary-coverage-threshold 10 \
        --summary-coverage-threshold 15 --summary-coverage-threshold 20 \
        --summary-coverage-threshold 25 --summary-coverage-threshold 30 \
        --summary-coverage-threshold 35 --summary-coverage-threshold 40 \
        --summary-coverage-threshold 45 --summary-coverage-threshold 50
done
```

### 7b. Flagstat QC

Runs [samtools flagstat](http://www.htslib.org/doc/samtools-flagstat.html) as an additional QC checkpoint, run both before and after BQSR for comparison.

```bash
#!/bin/sh

module load samtools

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')
    samtools flagstat "$SAMPLE".markdup.sorted.bam > "$SAMPLE".flagstat.txt
done
```

## 8. Base Quality Score Recalibration

Recalibrates base quality scores using known-site variant databases (dbSNP, Mills/1000 Genomes gold-standard indels, GATK known indels, HapMap), applies the recalibration, and generates before/after covariate plots. See [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360050815072-BaseRecalibrator), [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR), and [AnalyzeCovariates](https://gatk.broadinstitute.org/hc/en-us/articles/360051304351-AnalyzeCovariates) for details, and [this GATK discussion](https://gatk.broadinstitute.org/hc/en-us/community/posts/360075305092-Known-Sites-for-BQSR) on choosing known-sites resources.

```bash
#!/bin/sh

module load gatk
module load samtools
module load R

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')

    gatk BaseRecalibrator \
        -I "$SAMPLE".rgsorted.bam \
        -R <reference_genome>.fa \
        --known-sites <dbSNP>.vcf \
        --known-sites <Mills_and_1000G_gold_standard_indels>.vcf \
        --known-sites <known_indels>.vcf \
        --known-sites <hapmap>.vcf \
        -O "$SAMPLE".recal_data.table

    gatk ApplyBQSR \
        -R <reference_genome>.fa \
        -I "$SAMPLE".markdup.sorted.bam \
        --bqsr-recal-file "$SAMPLE".recal_data.table \
        -O "$SAMPLE".recal.bam

    gatk AnalyzeCovariates \
        -bqsr "$SAMPLE".recal_data.table \
        -plots "$SAMPLE".analyze_covariates.pdf

    samtools sort -@ <N> "$SAMPLE".recal.bam -o "$SAMPLE".recal.sorted.bam
    samtools index -@ <N> "$SAMPLE".recal.sorted.bam
done
```

### 8a. Depth of Coverage (post-BQSR)

Re-run DepthOfCoverage on the recalibrated BAMs. Using a BAM list as input (rather than looping per-sample) produces a single combined output across all samples, which is convenient for cohort-level coverage plots.

```bash
#!/bin/sh

module load gatk
module load samtools

gatk DepthOfCoverage \
    -R <reference_genome>.fa \
    -O <cohort>.depth_of_coverage \
    -I <bam_list>.list \
    -L <wgs_calling_regions>.interval_list \
    --summary-coverage-threshold 0 --summary-coverage-threshold 4 \
    --summary-coverage-threshold 6 --summary-coverage-threshold 10 \
    --summary-coverage-threshold 15 --summary-coverage-threshold 20 \
    --summary-coverage-threshold 25 --summary-coverage-threshold 30 \
    --summary-coverage-threshold 35 --summary-coverage-threshold 40 \
    --summary-coverage-threshold 45 --summary-coverage-threshold 50
```

## 9. Call Variants Per Sample

Calls variants per sample with [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360050814612-HaplotypeCaller) in GVCF mode, producing one `g.vcf.gz` per sample containing all sites (variant and reference) with no annotation/filtering applied yet. This is typically the most resource- and time-intensive step.

```bash
#!/bin/sh

module load gatk

for FASTQ in $(cat sample_list.txt); do
    SAMPLE=$(echo "$FASTQ" | awk -F "_" '{print $1}')

    gatk --java-options "-Xmx<heap_size>g" HaplotypeCaller \
        -R <reference_genome>.fa \
        -I "$SAMPLE".recal.sorted.bam \
        -O "$SAMPLE".g.vcf.gz \
        -ERC GVCF \
        -A AlleleFraction \
        -A BaseQuality \
        -A MappingQuality \
        --native-pair-hmm-threads <N>
done
```

Sample names embedded in each GVCF can be renamed and re-indexed (e.g., for compatibility with downstream tools like RVTests):

```bash
#!/bin/sh

module load picard

for SAMPLE in $(cat sample_names.txt); do
    java -jar <path_to>/picard.jar RenameSampleInVcf \
        INPUT="$SAMPLE".g.vcf.gz \
        OUTPUT="$SAMPLE".rename.g.vcf.gz \
        NEW_SAMPLE_NAME="$SAMPLE"
done
```

```bash
#!/bin/sh

module load gatk

for SAMPLE in $(cat sample_names.txt); do
    gatk IndexFeatureFile -I "$SAMPLE".rename.g.vcf.gz
done
```

## 10. Consolidate GVCFs

Merges all per-sample GVCFs into a single [GenomicsDB](https://gatk.broadinstitute.org/hc/en-us/articles/360051305591-GenomicsDBImport) workspace, which is used as the input for joint genotyping. The workspace directory should be created once and updated incrementally (via `--genomicsdb-update-workspace-path`) as new samples are added, rather than rebuilt from scratch each time.

**Note:** the `--tmp-dir` path must exist before submitting the job.

```bash
#!/bin/sh

module load samtools
module load picard
module load gatk

gatk --java-options "-Xmx<heap_size>g -Xms<heap_size>g" GenomicsDBImport \
    -V <sample_1>.g.vcf.gz \
    -V <sample_2>.g.vcf.gz \
    -V <sample_N>.g.vcf.gz \
    -L <wgs_calling_regions>.interval_list \
    --genomicsdb-workspace-path ./genomicsdb_workspace \
    --tmp-dir ./tmp \
    --reader-threads <N>
```

To add new samples to an existing workspace:

```bash
gatk --java-options "-Xmx<heap_size>g -Xms<heap_size>g" GenomicsDBImport \
    -V <new_sample_1>.g.vcf.gz \
    -V <new_sample_2>.g.vcf.gz \
    -L <wgs_calling_regions>.interval_list \
    --genomicsdb-update-workspace-path ./genomicsdb_workspace \
    --tmp-dir ./tmp \
    --reader-threads <N>
```

## 11. Joint Genotyping

Performs joint genotyping across all samples in the GenomicsDB workspace, producing a single cohort-level VCF. See [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360050816072-GenotypeGVCFs).

```bash
#!/bin/sh

module load gatk

gatk --java-options "-Xmx<heap_size>g" GenotypeGVCFs \
    -R <reference_genome>.fa \
    -V gendb://genomicsdb_workspace \
    -O joint_genotyped.vcf.gz \
    --tmp-dir ./tmp
```

## 12. Variant Quality Score Recalibration

Assigns a quality metric to each variant using machine-learning-based recalibration ([VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360050815872-VariantRecalibrator)), then filters the callset based on a sensitivity threshold ([ApplyVQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360051306591-ApplyVQSR)). Run separately for SNPs and indels.

```bash
#!/bin/sh

module load gatk
module load R

## Recalibrate SNPs
gatk VariantRecalibrator \
    -R <reference_genome>.fa \
    -V joint_genotyped.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 <Mills_and_1000G_gold_standard_indels>.vcf \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 <dbSNP>.vcf \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 <hapmap>.vcf \
    -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
    -mode SNP \
    -L <wgs_calling_regions>.interval_list \
    -O vqsr_snp.recal \
    --tranches-file vqsr_snp.tranches \
    --rscript-file vqsr_snp.plots.R

gatk ApplyVQSR \
    -R <reference_genome>.fa \
    -V joint_genotyped.vcf.gz \
    -O vqsr_snp_filtered.vcf.gz \
    --truth-sensitivity-filter-level 90.0 \
    --tranches-file vqsr_snp.tranches \
    --recal-file vqsr_snp.recal \
    -mode SNP

## Recalibrate indels
gatk VariantRecalibrator \
    -R <reference_genome>.fa \
    -V vqsr_snp_filtered.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 <Mills_and_1000G_gold_standard_indels>.vcf \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 <dbSNP>.vcf \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 <hapmap>.vcf \
    -an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum -an MQRankSum \
    -mode INDEL \
    -L <wgs_calling_regions>.interval_list \
    -O vqsr_indel.recal \
    --tranches-file vqsr_indel.tranches \
    --rscript-file vqsr_indel.plots.R

gatk ApplyVQSR \
    -R <reference_genome>.fa \
    -V vqsr_snp_filtered.vcf.gz \
    -O vqsr_snp_indel_filtered.vcf.gz \
    --truth-sensitivity-filter-level 90.0 \
    --tranches-file vqsr_indel.tranches \
    --recal-file vqsr_indel.recal \
    -mode INDEL
```

## 13. Annotate Variants

Annotates the final filtered callset with gene/transcript consequences, population allele frequencies, functional prediction scores, and clinical significance using [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/).

```bash
#!/bin/sh

module load perl
module load annovar

table_annovar.pl vqsr_snp_indel_filtered.vcf.gz <humandb_path> \
    -buildver hg38 \
    -protocol refGene,dbnsfp42a,dbscsnv11,revel,avsnp150,clinvar_20220320,gnomad30_genome,intervar_20180118 \
    --operation g,f,f,f,f,f,f,f --vcfinput \
    --argument '--hgvs --exonicsplicing --splicing_threshold 5',,,,,,, \
    --argument '--hgvs --indel_splicing_threshold 100',,,,,,, \
    --outfile annotated_final.vcf.gz
```

---

## Notes

- Known-sites databases (dbSNP, Mills/1000G gold-standard indels, HapMap, GATK known indels) should be matched to the same reference build (here, hg38/GRCh38) used for alignment. Newer database releases may use different contig-naming conventions than the reference FASTA — verify compatibility before swapping versions.
- Steps 6–9 are well suited to parallelization by splitting the sample cohort into batches/groups.
- For cohorts that grow over time, maintain a single running GenomicsDB workspace and update it incrementally (Step 10) rather than rebuilding it for each new batch, so all samples remain jointly genotyped together.
- QC checkpoints (FastQC, flagstat, DepthOfCoverage) are run both before and after key processing steps (trimming, BQSR) to track data quality and confirm genuine improvement.
