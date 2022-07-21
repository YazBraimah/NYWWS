"""
Author: Y. Ahmed-Braimah
--- NY Waste Water Surveillance pipeline
--- adapted from https://github.com/CFSAN-Biostatistics/C-WAP
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
GENEXUS_DNA = config['GENEXUS_DNA']
INDEX = config['INDEX']
script_dir = config['SCRIPT_DIR']
covidRefSequences = config['covidRefSequences']
k2db = config['K2_STD_DB_PATH']
pri_bed = config['PRIMER_BED']
allCOVdb = config['ALL_COVID']
majCOVdb = config['MAJOR_COVID']
VAR_DEF = config['VAR_DEF']
BAMS_DIR = config['BAMS_DIR']

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
        # expand(join(OUT_DIR, 'iVar', '{sample}.rawVarCalls.tsv'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_allCovid.out'), sample = SAMPLES),
        # expand(join(OUT_DIR, 'Sequences', '{sample}', 'consensus', 'consensus.fa'), sample = SAMPLES),
        join(OUT_DIR, 'Freyja', 'Aggregate', 'Genexus', 'freyja_stacked_barplots.png')
        # expand(join(OUT_DIR, 'QC', '{sample}', 'pos-coverage-quality.tsv'), sample = SAMPLES)

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# Rule to check raw SE read quality
rule FastQC:
    input:
        r1 = lambda wildcards: seFILES[wildcards.sample]['R1']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}' + '.fastQC_se.benchmark.tsv')
    message:
        """--- Checking read quality of sample "{wildcards.sample}" with FastQC """
    conda:
        'envs/fastqc_env.yml'
    shell:
        'fastqc -o ' + join(OUT_DIR, 'fastQC') + ' {input.r1} > {log} 2>&1'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule BWA_MEM:
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
        """--- Mapping reads for sample {wildcards.sample} with BWA MEM """
    run:
        shell('(bwa mem'
                ' -t 4'
                ' {params.index}'
                ' {input.r1}) 2>{log}'
                ' | samtools sort -@ 8 -o {output.bam} -')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule iVar_trimming:
    input:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.csorted.bwa.bam')
    params:
        pri_bed = pri_bed
    output:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.bwa.bam'),
        stats_cs = join(OUT_DIR, 'BWA', '{sample}', '{sample}.cssorted.stats'),
        stats_re = join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.stats')
    log:
        join(OUT_DIR, 'BWA', 'iVar_trimming_{sample}.log')
    benchmark:
        join(OUT_DIR, 'BWA', '{sample}', 'iVar_trimming.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=8000
    message:
        """--- Trimming inital BAM files for sample {wildcards.sample} """
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

rule BAM_to_FastQ:
    input:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.bwa.bam')
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
        """--- Converting trimmed BAM to FastQ for sample {wildcards.sample} """
    run:
        shell('samtools'
                ' view'
                ' --threads 2'
                ' {input.bam}'
                ' -o {wildcards.sample}.resorted.bwa.sam')
        shell(script_dir + '/sam2fastq.py'
                ' {wildcards.sample}.resorted.bwa.sam'
                ' {wildcards.sample}.resorted.fq')
        shell('gzip {wildcards.sample}.resorted.fq')
        shell('mv {wildcards.sample}.resorted.fq.gz ' + join(OUT_DIR, 'FastQs', '{wildcards.sample}'))
        shell('rm {wildcards.sample}.resorted.bwa.sam')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Samtools_pielup:
    input:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.bwa.bam')
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

rule iVar_variants:
    input:
        pileup = join(OUT_DIR, 'Pileups', '{sample}.up')
    params:
        crs = covidRefSequences,
        dna = DNA
    output:
        countsSt = join(OUT_DIR, 'iVar', '{sample}.rawVarCalls.tsv')
    log:
        join(OUT_DIR, 'iVar', '{sample}.log')
    benchmark:
        join(OUT_DIR, 'iVar', '{sample}.benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=16000
    message:
        """--- Calling variants for sample {wildcards.sample} with iVar. """
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

rule Kraken:
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
        "envs/kraken_env.yml"
    shell:
        'kraken2 --db {params.k2db} --threads 8 --report {output} {input.r1} > /dev/null'


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

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
        """--- QC for iVar variant calls of sample "{wildcards.sample}"."""
    run:
        shell('python3 ' + script_dir + '/plotQC.py'
                ' {input.pileup}'
                ' {params.pri_bed}'
                ' ' + join(OUT_DIR, 'QC', '{wildcards.sample}') +
                ' > {log} 2>&1')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Kraken_variants:
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
    conda:
        "envs/kraken_env.yml"
    message:
        """--- Discovering variants with Kraken for "{wildcards.sample}"."""
    shell:
        'kraken2 {input.fastq} --db {params.allCOVdb} --threads 4 --report {output.allCov} > /dev/null && '

        '/home/software/Bracken/bracken'
                ' -d {params.allCOVdb}'
                ' -i {output.allCov}'
                ' -o {output.allCovid_bracken}'
                ' -l P'
                ' > {log.all_brak} 2>&1 && '
        'kraken2'
                ' {input.fastq}'
                ' --db {params.majCOVdb}'
                ' --threads 4'
                ' --report {output.majCov}'
                ' > /dev/null && '
        '/home/software/Bracken/bracken'
                ' -d {params.majCOVdb}'
                ' -i {output.majCov}'
                ' -o {output.majCovid_bracken}'
                ' -l C'
                ' > {log.maj_brak} 2>&1'


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Consensus_sequence:
    input:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.bwa.bam')
    params:
        dna = DNA
    output:
        vcf = join(OUT_DIR, 'Sequences', '{sample}', 'consensus', 'calls.vcf.gz'),
        consensus = join(OUT_DIR, 'Sequences', '{sample}', 'consensus', 'consensus.fa')
    log:
        join(OUT_DIR, 'Sequences', '{sample}', 'consensus', 'calls.log')
    benchmark:
        join(OUT_DIR, 'Sequences', '{sample}', 'consensus', 'consensus_benchmark.tsv')
    message:
        """--- Outputting consensus sequence for sample "{wildcards.sample}" """
    conda:
        "envs/bcftools_env.yml"
    shell:
        'bcftools mpileup'
                ' -d 10000'
                ' -Ou'
                ' -f {params.dna}'
                ' {input.bam} |'
                ' bcftools call'
                ' --ploidy 1'
                ' -mv'
                ' -Oz'
                ' -o {output.vcf}'
                ' > {log} 2>&1 && '
        'bcftools index {output.vcf} && '
        'cat {params.dna} | bcftools consensus {output.vcf} > {output.consensus}'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Panglin_variants:
    input:
        consensus = join(OUT_DIR, 'Sequences', '{sample}', 'consensus', 'consensus.fa')
    output:
        join(OUT_DIR, 'Pangolin', '{sample}', 'lineage_report.csv')
    log:
        join(OUT_DIR, 'Pangolin', '{sample}', 'pangolin.log')
    benchmark:
        join(OUT_DIR, 'Pangolin', '{sample}', 'pangolin_benchmark.tsv')
    message:
        """--- Running Pangolin on sample "{wildcards.sample}" """
    run:
        shell('pangolin'
                ' {input.consensus}'
                ' --outfile {output}'
                ' --threads 4'
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule LDVC_variants:
    input:
        countsSt = join(OUT_DIR, 'iVar', '{sample}.rawVarCalls.tsv')
    params:
        var_def = VAR_DEF
    output:
        join(OUT_DIR, 'LDVC', '{sample}', 'ldvc_report.csv')
    log:
        join(OUT_DIR, 'LDVC', '{sample}', 'ldvc.log')
    benchmark:
        join(OUT_DIR, 'LDVC', '{sample}', 'ldvc_benchmark.tsv')
    message:
        """--- Running LDVC for sample "{wildcards.sample}".  """
    run:
        shell(script_dir + '/deconvolveVariants.py'
                ' {input.countsSt}'
                ' ' + join(OUT_DIR, 'LDVC', '{wildcards.sample}') +
                ' {params.var_def}'
                ' > {output}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_variants:
    input:
        bwa_bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.bwa.bam'),
        gen_bam = join(BAMS_DIR, '{sample}.ptrim.bam')
    params:
        dna = DNA,
        gdna = GENEXUS_DNA
    output:
        bwa_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'BWA', '{sample}.freyja.variants.tsv'),
        bwa_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'BWA', '{sample}.freyja.depths.tsv'),
        gen_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'Genexus', '{sample}.freyja.variants.tsv'),
        gen_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'Genexus', '{sample}.freyja.depths.tsv')
    log:
        join(OUT_DIR, 'Freyja', 'Variants', '{sample}_freyja_var-dep.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Variants', '{sample}_freyja_var-dep_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja variants for sample "{wildcards.sample}". """
    run:
        shell('freyja variants'
                ' {input.bwa_bam}'
                ' --variants {output.bwa_variants}'
                ' --depths {output.bwa_depths}'
                ' --ref {params.dna}')
        shell('freyja variants'
                ' {input.gen_bam}'
                ' --variants {output.gen_variants}'
                ' --depths {output.gen_depths}'
                ' --ref {params.gdna}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_demix:
    input:
        bwa_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'BWA', '{sample}.freyja.variants.tsv'),
        bwa_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'BWA', '{sample}.freyja.depths.tsv'),
        gen_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'Genexus', '{sample}.freyja.variants.tsv'),
        gen_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'Genexus', '{sample}.freyja.depths.tsv')
    output:
        bwa_demix = join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'BWA', '{sample}_freyja.demix'),
        gen_demix = join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'Genexus', '{sample}_freyja.demix')
    log:
        join(OUT_DIR, 'Freyja', 'Demix', '{sample}_freyja_demix.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Demix', '{sample}_freyja_demix_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja demix for sample "{wildcards.sample}". """
    run:
        shell('freyja demix'
                ' {input.bwa_variants}'
                ' {input.bwa_depths}'
                ' --output {output.bwa_demix}'
                ' --confirmedonly')
        shell('freyja demix'
                ' {input.gen_variants}'
                ' {input.gen_depths}'
                ' --output {output.gen_demix}'
                ' --confirmedonly')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_boot:
    input:
        bwa_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'BWA', '{sample}.freyja.variants.tsv'),
        bwa_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'BWA', '{sample}.freyja.depths.tsv'),
        gen_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'Genexus', '{sample}.freyja.variants.tsv'),
        gen_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', 'Genexus', '{sample}.freyja.depths.tsv'),
        bwa_demix = join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'BWA', '{sample}_freyja.demix'),
        gen_demix = join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'Genexus', '{sample}_freyja.demix')
    output:
        bwa_boot = join(OUT_DIR, 'Freyja', 'Boot', 'Results', 'BWA', '{sample}_freyja_boot_lineages.csv'),
        bwa_png =  join(OUT_DIR, 'Freyja', 'Boot', 'Results', 'BWA', '{sample}_freyja_bootstrap.png'),
        gen_boot = join(OUT_DIR, 'Freyja', 'Boot', 'Results', 'Genexus', '{sample}_freyja_boot_lineages.csv'),
        gen_png =  join(OUT_DIR, 'Freyja', 'Boot', 'Results', 'Genexus', '{sample}_freyja_bootstrap.png')
    log:
        join(OUT_DIR, 'Freyja', 'Boot', '{sample}_freyja_boot.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Boot', '{sample}_freyja_boot_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja boot for sample "{wildcards.sample}". """
    run:
        shell('freyja boot'
                ' {input.bwa_variants}'
                ' {input.bwa_depths}'
                ' --nt 4'
                ' --nb 10'
                ' --output_base {wildcards.sample}_bwa_freyja_boot'
                ' &&'
                ' mv {wildcards.sample}_bwa_freyja_boot*csv ' + join(OUT_DIR, 'Freyja', 'Boot', 'Results', 'BWA'))
        shell(script_dir + 'parseFreyjaBootstraps.py {input.bwa_demix} {output.bwa_boot} {output.bwa_png}')

        shell('freyja boot'
                ' {input.gen_variants}'
                ' {input.gen_depths}'
                ' --nt 4'
                ' --nb 10'
                ' --output_base {wildcards.sample}_gen_freyja_boot'
                ' &&'
                ' mv {wildcards.sample}_gen_freyja_boot*csv ' + join(OUT_DIR, 'Freyja', 'Boot', 'Results', 'Genexus'))
        shell(script_dir + 'parseFreyjaBootstraps.py {input.gen_demix} {output.gen_boot} {output.gen_png}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_aggregate:
    input:
        expand(join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'BWA', '{sample}_freyja.demix'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'Genexus', '{sample}_freyja.demix'), sample = SAMPLES)
    output:
        bwa_agg = join(OUT_DIR, 'Freyja', 'Aggregate', 'BWA', 'aggregated_results.tsv'),
        bwa_png = join(OUT_DIR, 'Freyja', 'Aggregate', 'BWA', 'freyja_stacked_barplots.png'),
        gen_agg = join(OUT_DIR, 'Freyja', 'Aggregate', 'Genexus', 'aggregated_results.tsv'),
        gen_png = join(OUT_DIR, 'Freyja', 'Aggregate', 'Genexus', 'freyja_stacked_barplots.png')
    log:
        join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_agg.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_agg_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja aggregate and outputting barplots. """
    run:
        shell('freyja aggregate'
                ' --output {output.bwa_agg}'
                ' ' + join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'BWA') + '/')
        shell('freyja plot'
                ' {output.bwa_agg}'
                ' --output {output.bwa_png}')
        shell('freyja aggregate'
                ' --output {output.gen_agg}'
                ' ' + join(OUT_DIR, 'Freyja', 'Demix', 'Results', 'Genexus') + '/')
        shell('freyja plot'
                ' {output.gen_agg}'
                ' --output {output.gen_png}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##


rule QualiMap:
    input:
        bam = join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.bwa.bam'),
        gtf = covidRefSequences
    output:
        join(OUT_DIR, 'BWA', '{sample}', '{sample}.qualimap', 'qualimapReport.html')
    log:
        join(OUT_DIR, 'BWA', '{sample}', 'qualmap.log')
    benchmark:
        join(OUT_DIR, 'BWA', '{sample}', 'qualmap_benchmark.tsv')
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
                ' -outdir ' + join(OUT_DIR, 'BWA', '{wildcards.sample}', '{wildcards.sample}.qualimap') +
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Rule to collate fastQC and HISAT2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '.R1_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_std.out'), sample = SAMPLES),
        expand(join(OUT_DIR, 'BWA', '{sample}', '{sample}.resorted.stats'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Pangolin', '{sample}', 'lineage_report.csv'), sample = SAMPLES),
        expand(join(OUT_DIR, 'LDVC', '{sample}', 'ldvc_report.csv'), sample = SAMPLES),
        expand(join(OUT_DIR, 'BWA', '{sample}', '{sample}.csorted.bwa.bam'), sample = SAMPLES),
        expand(join(OUT_DIR, 'BWA', '{sample}', '{sample}.qualimap', 'qualimapReport.html'), sample = SAMPLES)
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
        # shell('ls -1 ' + join(OUT_DIR) + '/BWA/*log > ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/BWA/*/*.qualimap | grep ":" | sed "s/://g" >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -l -dd 2 ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt') +
                ' > {log} 2>&1')
