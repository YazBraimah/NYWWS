"""
Author: Y. Ahmed-Braimah
--- Map RNA-seq reads with Hisat2 and annotate genome
--- with Trinotate.
"""

#########
############
## fix the kraken2/perl incompatibility problem (use conda directive for all rules except kraken2.
## REMAINING RULES
#1. LCSvariantCaller
#2. kallistonVariantCaller


import json
import os
import re
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

## define environment variables

##--------------------------------------------------------------------------------------##
## Global config files:
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']
INDEX = config['INDEX']
script_dir = config['SCRIPT_DIR']
covidRefSequences = config['covidRefSequences']
k2db = config['K2_STD_DB_PATH']
pri_bed = config['PRIMER_BED']
allCOVdb = config['ALL_COVID']
majCOVdb = config['MAJOR_COVID']
VAR_DEF = config['VAR_DEF']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
# HOME_DIR = config['HOME_DIR']  # the "launch_snakemake.sh" and "config.yml" files should be here

# Samples and their corresponding filenames.
# single-end
seFILES = json.load(open(config['SE_SAMPLES_JSON']))
seSAMPLES = sorted(seFILES.keys())
# paired-end:
peFILES = json.load(open(config['PE_SAMPLES_JSON']))
peSAMPLES = sorted(peFILES.keys())

combinedSam = [peSAMPLES, seSAMPLES]
SAMPLES = [y for x in combinedSam for y in x]

## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

##--------------------------------------------------------------------------------------##
#
# _____ _             _               _               _
#|  ___(_)_ __   __ _| |   ___  _   _| |_ _ __  _   _| |_ ___
#| |_  | | '_ \ / _` | |  / _ \| | | | __| '_ \| | | | __/ __|
#|  _| | | | | | (_| | | | (_) | |_| | |_| |_) | |_| | |_\__ \
#|_|   |_|_| |_|\__,_|_|  \___/ \__,_|\__| .__/ \__,_|\__|___/
#                                        |_|
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all:
    input:
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html'),
        expand(join(OUT_DIR, 'Variants', 'iVar', '{sample}.rawVarCalls.tsv'), sample = SAMPLES),
        # expand(join(OUT_DIR, 'Kraken', '{sample}', '{sample}.majCovid.bracken'), sample = SAMPLES),
        expand(join(OUT_DIR, 'iVar', '{sample}', 'consensus', 'consensus.fa'), sample = SAMPLES),
        expand(join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_bootstrap.png'), sample = SAMPLES),
        expand(join(OUT_DIR, 'QC', '{sample}', 'pos-coverage-quality.tsv'), sample = SAMPLES)

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##
#  ___   ____
# / _ \ / ___|
#| | | | |
#| |_| | |___
# \__\_\\____|
#
##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# Rule to check raw SE read quality
rule fastqcSE:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.benchmark.tsv')
    message:
        """--- Checking read quality of SE sample "{wildcards.sample}" with FastQC """
    run:
        shell('/home/yahmed/software/FastQC/fastqc'
                ' -o ' + join(OUT_DIR, 'fastQC') +
                ' {input.r1}'
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to map reads with BWA
rule bwa_mem:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    params:
        index = INDEX
    output:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.csorted.bwa.bam')
    log:
        join(OUT_DIR, 'BWA', 'bwa_{sample}.log')
    benchmark:
        join(OUT_DIR, 'BWA', '{sample}', 'bwa_map_se.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=8000
    message:
        """--- Mapping SE reads for sample {wildcards.sample} to genome with Bowtie2 """
    run:
        shell('(bwa mem'
                ' -t 4'
                ' {params.index}'
                ' {input.r1}) 2>{log}'
                ' | samtools sort -@ 8 -o {output.bam} -')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to map PE reads with HISAT2
rule iVar_trimming: ## Need to include the primer bed files later
    input:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.csorted.bwa.bam')
    params:
        pri_bed = pri_bed
    output:
        bam = join(OUT_DIR, 'iVar', '{sample}', '{sample}.resorted.bwa.bam'),
        stats_cs = join(OUT_DIR, 'iVar', '{sample}', '{sample}.cssorted.stats'),
        stats_re = join(OUT_DIR, 'iVar', '{sample}', '{sample}.resorted.stats')
    log:
        join(OUT_DIR, 'iVar', 'bwa_{sample}.log')
    benchmark:
        join(OUT_DIR, 'iVar', '{sample}', 'iVar.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=8000
    message:
        """--- Trimming round of inital BAM files for sample {wildcards.sample} """
    run:
        shell('ivar'
                ' trim'
                ' -e -b {params.pri_bed}'
                ' -p trimmed.{wildcards.sample}'
                ' -i {input.bam}')
        shell('samtools'
                ' sort'
                ' trimmed.{wildcards.sample}.bam'
                ' -o {output.bam}'
                ' -@ 8')
        shell('samtools'
                ' stats'
                ' {input.bam}'
                ' | grep ^SN'
                ' | cut -f 2-'
                ' > {output.stats_cs}')
        shell('samtools'
                ' stats'
                ' {output.bam}'
                ' | grep ^SN'
                ' | cut -f 2-'
                ' > {output.stats_re}')
        shell('rm trimmed.{wildcards.sample}.bam')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule BAM_to_FastQ: ## May need to incorporate a subsampling step if the read numbers are high
    input:
        bam = join(OUT_DIR, 'iVar', '{sample}', '{sample}.resorted.bwa.bam')
    output:
        fastq = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.resorted.fq.gz')
    log:
        join(OUT_DIR, 'FastQs', '{sample}', 'bam2fq.log')
    benchmark:
        join(OUT_DIR, 'FastQs', '{sample}', 'bam2fq.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=8000
    message:
        """--- Converting BAM file to FastQ for sample {wildcards.sample} """
    run:
        shell('samtools'
                ' view'
                ' --threads 2'
                ' {input.bam}'
                ' -o {wildcards.sample}.resorted.bwa.sam')
        shell('sam2fastq.py'
                ' {wildcards.sample}.resorted.bwa.sam'
                ' {wildcards.sample}.resorted.fq')
        shell('gzip {wildcards.sample}.resorted.fq')
        shell('mv {wildcards.sample}.resorted.fq.gz ' + join(OUT_DIR, 'FastQs', '{wildcards.sample}'))
        shell('rm {wildcards.sample}.resorted.bwa.sam')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule generate_pielup:
    input:
        bam = join(OUT_DIR, 'iVar', '{sample}', '{sample}.resorted.bwa.bam')
    params:
        dna = DNA
    output:
        pileup = join(OUT_DIR, 'Pileups', '{sample}.up')
    log:
        join(OUT_DIR, 'Pileups', '{sample}.pielup.log')
    benchmark:
        join(OUT_DIR, 'Pileups', '{sample}.pielup.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=16000
    message:
        """--- Generating pielups for sample {wildcards.sample} """
    run:
        shell ('samtools'
                ' mpileup'
                ' -aa'
                ' -A'
                ' -d 10000'
                ' -B'
                ' -Q 0'
                ' --reference {params.dna}'
                ' -o {output.pileup}'
                ' {input.bam}')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## measure transcript abundances with Stringtie
rule call_variants:
    input:
        pileup = join(OUT_DIR, 'Pileups', '{sample}.up')
    params:
        crs = covidRefSequences,
        dna = DNA
    output:
        countsSt = join(OUT_DIR, 'Variants', 'iVar', '{sample}.rawVarCalls.tsv')
    log:
        join(OUT_DIR, 'Variants', 'iVar', '{sample}.log')
    benchmark:
        join(OUT_DIR, 'Variants', 'iVar', '{sample}.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=16000
    message:
        """--- Calling variants for sample {wildcards.sample}. """
    run:
        shell('cat {input.pileup} |'
                ' ivar variants'
                ' -p {wildcards.sample}.rawVarCalls'
                ' -g {params.crs}'
                ' -r {params.dna}'
                ' -m 10'
                ' > {log} 2>&1'
                ' && mv {wildcards.sample}.rawVarCalls.tsv {output.countsSt}')



##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##


## Rule for mapping PE reads to the genome with Bowtie2
rule kraken:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    params:
        k2db = k2db
    output:
        join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_std.out')
    log:
        join(OUT_DIR, 'Kraken', '{sample}', 'k2_std.log')
    benchmark:
        join(OUT_DIR, 'Kraken', '{sample}', 'k2_std_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Kraken2 search for sample "{wildcards.sample}"."""
    conda:
        "kraken2"
    shell:
        'kraken2 --db {params.k2db} --threads 8 --report {output} {input.r1}> /dev/null'
    # run:
    #     shell('/home/yahmed/GitHub_Repositories/kraken2/kraken2'
    #             ' {input.r1}'
    #             ' --db {params.k2db}'
    #             ' --threads 8'
    #             ' --report {output}'
    #             ' > /dev/null')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##


## Rule for mapping PE reads to the genome with Bowtie2
rule QCplots:
    input:
        pileup = join(OUT_DIR, 'Pileups', '{sample}.up')
    params:
        pri_bed = pri_bed
    output:
        join(OUT_DIR, 'QC', '{sample}', 'pos-coverage-quality.tsv')
    log:
        join(OUT_DIR, 'QC', '{sample}', 'qc.log')
    benchmark:
        join(OUT_DIR, 'QC', '{sample}', 'qc_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=8000
    message:
        """--- QC for variant calls "{wildcards.sample}"."""
    run:
        shell('python3 ' + script_dir + '/plotQC.py'
                ' {input.pileup}'
                ' {params.pri_bed}'
                ' ' + join(OUT_DIR, 'QC', '{wildcards.sample}') +
                ' > {log} 2>&1')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##


rule krakenVariantCaller:
    input:
        fastq = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.resorted.fq.gz')
    params:
        allCOVdb = allCOVdb,
        majCOVdb = majCOVdb
    output:
        allCov = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_allCovid.out'),
        allCovid_bracken = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.allCovid.bracken'),
        majCov = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_majCovid.out'),
        majCovid_bracken = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.majCovid.bracken')
    log:
        all_brak = join(OUT_DIR, 'Kraken', '{sample}', 'k2_std.log'),
        maj_brak = join(OUT_DIR, 'Kraken', '{sample}', 'k2_std.log')
    benchmark:
        join(OUT_DIR, 'Kraken', '{sample}', 'k2_std_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Kraken2 search for sample "{wildcards.sample}"."""
    run:
        shell('/home/yahmed/GitHub_Repositories/kraken2/kraken2'
                ' {input.fastq}'
                ' --db {params.allCOVdb}'
                ' --threads 4'
                ' --report {output.allCov}'
                ' > /dev/null')
        shell('/home/yahmed/GitHub_Repositories/Bracken/bracken'
                ' -d {params.allCOVdb}'
                ' -i {output.allCov}'
                ' -o {output.allCovid_bracken}'
                ' -l P'
                ' > {log.all_brak} 2>&1')
        shell('/home/yahmed/GitHub_Repositories/kraken2/kraken2'
                ' {input.fastq}'
                ' --db {params.majCOVdb}'
                ' --threads 4'
                ' --report {output.majCov}'
                ' > /dev/null')
        shell('/home/yahmed/GitHub_Repositories/Bracken/bracken'
                ' -d {params.majCOVdb}'
                ' -i {output.majCov}'
                ' -o {output.majCovid_bracken}'
                ' -l C'
                ' > {log.maj_brak} 2>&1')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule consensusSequence:
    input:
        bam = join(OUT_DIR, 'iVar', '{sample}', '{sample}.resorted.bwa.bam')
    params:
        dna = DNA
    output:
        vcf = join(OUT_DIR, 'iVar', '{sample}', 'consensus', 'calls.vcf.gz'),
        consensus = join(OUT_DIR, 'iVar', '{sample}', 'consensus', 'consensus.fa')
    log:
        join(OUT_DIR, 'iVar', '{sample}', 'consensus', 'calls.log')
    benchmark:
        join(OUT_DIR, 'iVar', '{sample}', 'consensus', 'consensus_benchmark.tsv')
    message:
        """--- Running MultiQC """
    run:
        shell('bcftools mpileup'
                ' -d 10000'
                ' -Ou'
                ' -f {params.dna}'
                ' {input.bam} |'
                ' bcftools call'
                ' --ploidy 1'
                ' -mv'
                ' -Oz'
                ' -o {output.vcf}'
                ' > {log} 2>&1')
        shell('bcftools index {output.vcf}')
        shell('cat {params.dna} | bcftools consensus {output.vcf} > {output.consensus}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule pangolinVariantCaller:
    input:
        consensus = join(OUT_DIR, 'iVar', '{sample}', 'consensus', 'consensus.fa')
    output:
        join(OUT_DIR, 'Pangolin', '{sample}', 'lineage_report.csv')
    log:
        join(OUT_DIR, 'Pangolin', '{sample}', 'pangolin.log')
    benchmark:
        join(OUT_DIR, 'Pangolin', '{sample}', 'pangolin_benchmark.tsv')
    message:
        """--- Running Pangolin """
    run:
        shell('pangolin'
                ' {input.consensus}'
                ' --outfile {output}'
                ' --threads 4'
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule linearDeconVariantCaller:
    input:
        countsSt = join(OUT_DIR, 'Variants', 'iVar', '{sample}.rawVarCalls.tsv')
    params:
        var_def = VAR_DEF
    output:
        join(OUT_DIR, 'LDVC', '{sample}', 'ldvc_report.csv')
    log:
        join(OUT_DIR, 'LDVC', '{sample}', 'ldvc.log')
    benchmark:
        join(OUT_DIR, 'LDVC', '{sample}', 'ldvc_benchmark.tsv')
    message:
        """--- Running LDVC """
    run:
        shell('deconvolveVariants.py'
                ' {input.countsSt}'
                ' ' + join(OUT_DIR, 'LDVC', '{wildcards.sample}') +
                ' {params.var_def}'
                ' > {output}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule freyjaVariantCaller:
    input:
        bam = join(OUT_DIR, 'iVar', '{sample}', '{sample}.resorted.bwa.bam')
    params:
        dna = DNA
    output:
        variants = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.variants.tsv'),
        depths = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.depths.tsv')
    log:
        join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_variant.log')
    benchmark:
        join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_variant_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja variants"""
    run:
        shell('freyja variants'
                ' {input.bam} '
                ' --variants {output.variants}'
                ' --depths {output.depths}'
                ' --ref {params.dna}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule freyjaDemix:
    input:
        variants = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.variants.tsv'),
        depths = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.depths.tsv')
    output:
        demix = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.demix')
    log:
        join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_demix.log')
    benchmark:
        join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_demi_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja """
    run:
        shell('freyja demix'
                ' {input.variants}'
                ' {input.depths}'
                ' --output {output.demix}'
                ' --confirmedonly')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule freyjaBoot:
    input:
        variants = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.variants.tsv'),
        depths = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.depths.tsv'),
        demix = join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja.demix')
    output:
        boot = join(OUT_DIR, 'iVar', '{sample}', 'freyja', '{sample}_freyja_boot_lineages.csv'),
        png =  join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_bootstrap.png')
    log:
        join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_boot.log')
    benchmark:
        join(OUT_DIR, 'iVar', '{sample}', 'freyja', 'freyja_boot_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja boot"""
    run:
        shell('freyja boot'
                ' {input.variants}'
                ' {input.depths}'
                ' --nt 4'
                ' --nb 10'
                ' --output_base {wildcards.sample}_freyja_boot'
                ' &&'
                ' mv {wildcards.sample}_freyja_boot*csv ' + join(OUT_DIR, 'iVar', '{wildcards.sample}', 'freyja'))
        shell('parseFreyjaBootstraps.py {input.demix} {output.boot} {output.png}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##


## Rule for mapping PE reads to the genome with Bowtie2
rule qualiMap:
    input:
        bam = join(OUT_DIR, 'iVar', '{sample}', '{sample}.resorted.bwa.bam'),
        gtf = covidRefSequences
    output:
        join(OUT_DIR, 'iVar', '{sample}', '{sample}.qualimap', 'qualimapReport.html')
    log:
        join(OUT_DIR, 'iVar', '{sample}', 'qualmap.log')
    benchmark:
        join(OUT_DIR, 'iVar', '{sample}', 'qualmap_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Evaluating mapping quality with QualiMap for sample "{wildcards.sample}"."""
    run:
        shell('qualimap bamqc'
                ' -bam {input.bam}'
                ' --java-mem-size=32G'
                # ' -gff {input.gtf}'
                ' -outdir ' + join(OUT_DIR, 'iVar', '{wildcards.sample}', '{wildcards.sample}.qualimap') +
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_std.out'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Pangolin', '{sample}', 'lineage_report.csv'), sample = SAMPLES),
        expand(join(OUT_DIR, 'LDVC', '{sample}', 'ldvc_report.csv'), sample = SAMPLES),
        expand(join(OUT_DIR, 'BWA', '{sample}', '{sample}.csorted.bwa.bam'), sample = SAMPLES),
        expand(join(OUT_DIR, 'iVar', '{sample}', '{sample}.qualimap', 'qualimapReport.html'), sample = SAMPLES)
    output:
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC', 'multiQC.benchmark.tsv')
    message:
        """--- Running MultiQC """
    run:
        shell('ls -1 ' + join(OUT_DIR) + '/fastQC/*fastqc.zip >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/BWA/*log > ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/iVar/*/*.qualimap | grep ":" | sed "s/://g" >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -l -dd 2 ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt') +
                ' > {log} 2>&1')
