"""
Author: Y. Ahmed-Braimah
--- NY Waste Water Surveillance pipeline
--- adapted from https://github.com/CFSAN-Biostatistics/C-WAP
"""

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
KalIdx = config['KalIdx']
LCS_DIR = config['LCS_DIR']
SAM_HEADER = config['SAM_HEADER']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
HOME_DIR = config['HOME_DIR']

# Samples and their corresponding filenames.
# single-end
# FILES = json.load(open(config['SE_SAMPLES_JSON']))
# seSAMPLES = sorted(FILES.keys())
# paired-end:
# peFILES = json.load(open(config['PE_SAMPLES_JSON']))
# peSAMPLES = sorted(peFILES.keys())

# combinedSam = [peSAMPLES, seSAMPLES]
# SAMPLES = [y for x in combinedSam for y in x]

FILES = json.load(open(config['BAM_SAMPLES_JSON']))
SAMPLES = sorted(FILES.keys())

## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all:
    input:
        join(OUT_DIR, 'MultiQC', 'multiqc_report.html'),
        join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_stacked_barplots.png'),
        expand(join(OUT_DIR, 'Freyja', 'Boot', 'Results', '{sample}_freyja_bootstrap.png'), sample = SAMPLES),
        # expand(join(OUT_DIR, 'SummaryPie', '{sample}', 'pieChart_lcs.png'), sample = SAMPLES),
        expand(join(OUT_DIR, 'QC', '{sample}', 'pos-coverage-quality.tsv'), sample = SAMPLES)

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule BAM_to_FastQ:
    input:
        bam = join(BAMS_DIR, '{sample}.ptrim.bam')
    output:
        proc_sam = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.sam'),
        proc_bam = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.bam'),
        fastq = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.fq.gz'),
        stats_cs = join(OUT_DIR, 'BAM', '{sample}', '{sample}.cssorted.stats'),
        stats_re = join(OUT_DIR, 'BAM', '{sample}', '{sample}.resorted.stats')
    params:
        header = SAM_HEADER,
        dna = DNA
    log:
        join(OUT_DIR, 'FastQs', '{sample}', 'bam2fq.log')
    benchmark:
        join(OUT_DIR, 'FastQs', '{sample}', 'bam2fq.benchmark.tsv')
    threads:
        4
    resources:
        mem_mb=8000
    message:
        """--- Converting Genexus BAM to FastQ for sample {wildcards.sample} """
    run:
        shell('samtools index'
                ' {input.bam}')
        shell('samtools'
                ' view'
                ' --threads 2'
                ' {input.bam}'
                ' 2019-nCoV'
                ' -b'
                ' | samtools reheader'
                ' {params.header}'
                ' -'
                ' | samtools view'
                ' -b'
                ' - |'
                ' samtools sort'
                ' -@ 8'
                ' -'
                ' -o {output.proc_bam}')
        shell('samtools view'
                ' {output.proc_bam}'
                ' -o {output.proc_sam}')
        shell('samtools'
                ' stats'
                ' {input.bam}'
                ' | grep ^SN'
                ' | cut -f 2-'
                ' > {output.stats_cs}')
        shell('samtools'
                ' stats'
                ' {output.proc_bam}'
                ' | grep ^SN'
                ' | cut -f 2-'
                ' > {output.stats_re}')
        shell(script_dir + '/sam2fastq.py'
                ' {output.proc_sam}'
                ' {wildcards.sample}.fq')
        shell('gzip {wildcards.sample}.fq')
        shell('mv {wildcards.sample}.fq.gz ' + join(OUT_DIR, 'FastQs', '{wildcards.sample}'))

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule FastQC:
    input:
        r1 = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.fq.gz')
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}' + '_fastqc.html')
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

rule Samtools_pielup:
    input:
        bam = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.bam')
    params:
        dna = DNA
    output:
        pileup = join(OUT_DIR, 'Pileups', '{sample}.up')
    log:
        join(OUT_DIR, 'Pileups', '{sample}.pielup.log')
    benchmark:
        join(OUT_DIR, 'Pileups', '{sample}.pielup.benchmark.tsv')
    threads:
        4
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
        4
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

# rule Kraken:
#     input:
#         r1 = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.fq.gz')
#     params:
#         k2db = k2db
#     output:
#         join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_std.out')
#     log:
#         join(OUT_DIR, 'Kraken', '{sample}', 'k2_std.log')
#     benchmark:
#         join(OUT_DIR, 'Kraken', '{sample}', 'k2_std_benchmark.tsv')
#     threads:
#         4
#     resources:
#         mem_mb=32000
#     message:
#         """--- Kraken2 search for sample "{wildcards.sample}"."""
#     conda:
#         "envs/kraken_env.yml"
#     shell:
#         'kraken2 --db {params.k2db} --threads 8 --report {output} {input.r1} > {log} 2>&1'

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
        4
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

# rule Kraken_variants:
#     input:
#         fastq = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.fq.gz')
#     params:
#         allCOVdb = allCOVdb,
#         majCOVdb = majCOVdb
#     output:
#         allCov = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_allCovid.out'),
#         allCovid_bracken = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_allCovid_bracken_phylums.out'),
#         majCov = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_majCovid.out'),
#         majCovid_bracken = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_majCovid_bracken_classes.out')
#     log:
#         all_brak = join(OUT_DIR, 'Kraken', '{sample}', 'k2_std_all.log'),
#         maj_brak = join(OUT_DIR, 'Kraken', '{sample}', 'k2_std_maj.log')
#     benchmark:
#         join(OUT_DIR, 'Kraken', '{sample}', 'k2_std_benchmark.tsv')
#     threads:
#         8
#     resources:
#         mem_mb=32000
#     conda:
#         "envs/kraken_env.yml"
#     message:
#         """--- Discovering variants with Kraken for "{wildcards.sample}"."""
#     shell:
#         '''
#             readNum=$(zcat {input.fastq} | grep read | wc -l);
#             if [[ $readNum -gt 5000 ]]
#             then
#             kraken2 {input.fastq} --db {params.allCOVdb} --threads 4 --report {output.allCov} > /dev/null &&
#             /home/software/Bracken/bracken \
#                      -d {params.allCOVdb} \
#                      -i {output.allCov} \
#                      -o {output.allCovid_bracken} \
#                      -l P \
#                      > {log.all_brak} 2>&1 &&
#             kraken2 \
#                      {input.fastq} \
#                      --db {params.majCOVdb} \
#                      --threads 4 \
#                      --report {output.majCov} \
#                      > /dev/null &&
#             /home/software/Bracken/bracken \
#                      -d {params.majCOVdb} \
#                      -i {output.majCov} \
#                      -o {output.majCovid_bracken} \
#                      -l C \
#                      > {log.maj_brak} 2>&1
#             else
#             echo -e "100.00\t0\t0\tR\t1\troot" > {output.allCovid_bracken}
#             echo -e "100.00\t0\t0\tR\t1\tError" >> {output.allCovid_bracken}
#             cp {output.allCovid_bracken} {output.majCovid_bracken}
#             cp {output.allCovid_bracken} {output.allCov}
#             cp {output.allCovid_bracken} {output.majCov}
#             fi
#             '''

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Consensus_sequence:
    input:
        bam = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.bam')
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

# rule LDVC_variants:
#     input:
#         countsSt = join(OUT_DIR, 'iVar', '{sample}.rawVarCalls.tsv')
#     params:
#         var_def = VAR_DEF
#     output:
#         join(OUT_DIR, 'LDVC', '{sample}', 'linearDeconvolution_abundance.csv')
#     log:
#         join(OUT_DIR, 'LDVC', '{sample}', 'ldvc.log')
#     benchmark:
#         join(OUT_DIR, 'LDVC', '{sample}', 'ldvc_benchmark.tsv')
#     message:
#         """--- Running LDVC for sample "{wildcards.sample}".  """
#     run:
#         shell(script_dir + '/deconvolveVariants.py'
#                 ' {input.countsSt}'
#                 ' ' + join(OUT_DIR, 'LDVC', '{wildcards.sample}') +
#                 ' {params.var_def}'
#                 ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# rule Kallisto_variants:
#     input:
#         fastq = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.fq.gz')
#     params:
#         kalidx = KalIdx
#     output:
#         join(OUT_DIR, 'Kallisto', '{sample}', 'abundance.tsv')
#     log:
#         join(OUT_DIR, 'Kallisto', '{sample}', 'kallisto.log')
#     benchmark:
#         join(OUT_DIR, 'Kallisto', '{sample}', 'kallisto_benchmark.tsv')
#     message:
#         """--- Running Kallisto for sample "{wildcards.sample}".  """
#     conda:
#         "envs/kallisto_env.yml"
#     shell:
#         '''
#             readNum=$(zcat {input.fastq} | grep read | wc -l);
#             if [[ $readNum -gt 5000 ]]
#             then
#             kallisto quant \
#             --index {params.kalidx} \
#             --output-dir + join(OUT_DIR, 'Kallisto', '{wildcards.sample}')  \
#             --plaintext -t 2 \
#             --single \
#             -l 300 \
#             -s 50 \
#             {input.fastq} \
#             > {log} 2>&1
#             else
#             echo -e "target_id\tlength\teff_length\test_counts\ttpm" > {output}
#             echo -e "Error\t29903\t29903\t100\t100" >> {output}
#             fi
#             '''


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# rule LCS_variants:
#     input:
#         fastq = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.fq.gz')
#     params:
#         lcs_dir = LCS_DIR
#     output:
#         join(OUT_DIR, 'LCS', '{sample}', 'outputs', 'decompose', 'lcs.out')
#     log:
#         join(OUT_DIR, 'LCS', '{sample}', 'lcs.log')
#     benchmark:
#         join(OUT_DIR, 'LCS', '{sample}', 'lcs_benchmark.tsv')
#     message:
#         """--- Running LCS for sample "{wildcards.sample}". """
#     conda:
#         "envs/lcs_env.yml"
#     shell:
#         'cp -r {params.lcs_dir} {wildcards.sample}_LCS &&'
#         ' cd {wildcards.sample}_LCS &&'
#         ' mkdir -p outputs/variants_table &&'
#         ' zcat data/pre-generated-marker-tables/pango-designation-markers-v1.2.124.tsv.gz > outputs/variants_table/pango-markers-table.tsv &&'
#         ' mkdir data/fastq &&'
#         ' cp {input.fastq} data/fastq/resorted.fastq.gz &&'
#         ' echo "resorted" > data/tags_pool_lcs &&'
#         ' snakemake --config markers=pango dataset=lcs --cores 2 &&'
#         ' sed -i "s/resorted/{wildcards.sample}/g" outputs/decompose/lcs.out &&'
#         ' cp -r outputs ' + join(OUT_DIR, 'LCS', '{wildcards.sample}') + '/'
#         ' && rm -rf ' + join(HOME_DIR, '{wildcards.sample}_LCS')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_update:
    output:
        done = join(OUT_DIR, 'Freyja', 'Update', 'update.ok')
    log:
        join(OUT_DIR, 'Freyja', 'Update', 'update.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Update', 'update_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=16000
    message:
        """--- Updating Freyja database. """
    run:
        shell('freyja update')
        shell('touch {output.done}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_variants:
    input:
        done = join(OUT_DIR, 'Freyja', 'Update', 'update.ok'),
        bt2_bam = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.bam')
    params:
        dna = DNA
    output:
        bt2_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', '{sample}.freyja.variants.tsv'),
        bt2_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', '{sample}.freyja.depths.tsv')
    log:
        join(OUT_DIR, 'Freyja', 'Variants', '{sample}_freyja_bt2_var-dep.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Variants', '{sample}_freyja_bt2_var-dep_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja Variants for sample "{wildcards.sample}". """
    run:
        shell('freyja variants'
                ' {input.bt2_bam}'
                ' --variants {output.bt2_variants}'
                ' --depths {output.bt2_depths}'
                ' --ref {params.dna}')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_demix:
    input:
        bt2_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', '{sample}.freyja.variants.tsv'),
        bt2_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', '{sample}.freyja.depths.tsv')
    output:
        bt2_demix = join(OUT_DIR, 'Freyja', 'Demix', 'Results', '{sample}_freyja.demix')
    log:
        join(OUT_DIR, 'Freyja', 'Demix', '{sample}_freyja_bt2_demix.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Demix', '{sample}_freyja_bt2_demix_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja Demix for sample "{wildcards.sample}". """
    run:
        shell('freyja demix'
                ' {input.bt2_variants}'
                ' {input.bt2_depths}'
                ' --output {output.bt2_demix}'
                # ' --confirmedonly'
                )

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_boot:
    input:
        bt2_variants = join(OUT_DIR, 'Freyja', 'Variants', 'Results', '{sample}.freyja.variants.tsv'),
        bt2_depths = join(OUT_DIR, 'Freyja', 'Variants', 'Results', '{sample}.freyja.depths.tsv'),
        bt2_demix = join(OUT_DIR, 'Freyja', 'Demix', 'Results', '{sample}_freyja.demix')
    output:
        bt2_boot = join(OUT_DIR, 'Freyja', 'Boot', 'Results', '{sample}_freyja_boot_lineages.csv'),
        bt2_png =  join(OUT_DIR, 'Freyja', 'Boot', 'Results', '{sample}_freyja_bootstrap.png')
    log:
        join(OUT_DIR, 'Freyja', 'Boot', '{sample}_freyja_bt2_boot.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Boot', '{sample}_freyja_bt2_boot_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja Boot for sample "{wildcards.sample}". """
    run:
        shell('freyja boot'
                ' {input.bt2_variants}'
                ' {input.bt2_depths}'
                ' --nt 4'
                ' --nb 10'
                ' --output_base {wildcards.sample}_freyja_boot'
                ' &&'
                ' mv {wildcards.sample}_freyja_boot*csv ' + join(OUT_DIR, 'Freyja', 'Boot', 'Results'))
        shell(script_dir + 'parseFreyjaBootstraps.py {input.bt2_demix} {output.bt2_boot} {output.bt2_png}')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule Freyja_aggregate:
    input:
        expand(join(OUT_DIR, 'Freyja', 'Demix', 'Results', '{sample}_freyja.demix'), sample = SAMPLES)
    output:
        bt2_agg = join(OUT_DIR, 'Freyja', 'Aggregate', 'aggregated_results.tsv'),
        bt2_png = join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_stacked_barplots.png')
    log:
        join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_agg.log')
    benchmark:
        join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_agg_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Running Freyja Aggregate and outputting barplots. """
    run:
        shell('freyja aggregate'
                ' --output {output.bt2_agg}'
                ' ' + join(OUT_DIR, 'Freyja', 'Demix', 'Results') + '/')
        shell('freyja plot'
                ' {output.bt2_agg}'
                ' --output {output.bt2_png}')


##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule QualiMap_BAM:
    input:
        bam = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.bam'),
        gtf = covidRefSequences
    output:
        join(OUT_DIR, 'BAM', '{sample}', '{sample}.qualimap', 'qualimapReport.html')
    log:
        join(OUT_DIR, 'BAM', '{sample}', 'qualmap.log')
    benchmark:
        join(OUT_DIR, 'BAM', '{sample}', 'qualmap_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Evaluating BAM mapping quality with QualiMap for sample "{wildcards.sample}"."""
    run:
        shell('qualimap bamqc'
                ' -bam {input.bam}'
                ' --java-mem-size=32G'
                # ' -gff {input.gtf}'
                ' -outdir ' + join(OUT_DIR, 'BAM', '{wildcards.sample}', '{wildcards.sample}.qualimap') +
                ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# rule QualiMap_Genexus:
#     input:
#         gen_bam = join(OUT_DIR, 'FastQs', '{sample}', '{sample}.bam')
#     output:
#         join(OUT_DIR, 'Genexus_BAM', '{sample}', '{sample}.qualimap', 'qualimapReport.html')
#     log:
#         join(OUT_DIR, 'Genexus_BAM', '{sample}', 'qualmap.log')
#     benchmark:
#         join(OUT_DIR, 'Genexus_BAM', '{sample}', 'qualmap_benchmark.tsv')
#     threads:
#         8
#     resources:
#         mem_mb=32000
#     message:
#         """--- Evaluating Genexus BAM with QualiMap for sample "{wildcards.sample}"."""
#     run:
#         shell('qualimap bamqc'
#                 ' -bam {input.bam}'
#                 ' -nt 12'
#                 ' --java-mem-size=32G'
#                 ' -outdir ' + join(OUT_DIR, 'Genexus_BAM', '{wildcards.sample}', '{wildcards.sample}.qualimap') +
#                 ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

# rule PieChart_summaries:
#     input:
#         ldvc = join(OUT_DIR, 'LDVC', '{sample}', 'linearDeconvolution_abundance.csv'),
#         kallisto = join(OUT_DIR, 'Kallisto', '{sample}', 'abundance.tsv'),
#         allCovid_bracken = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_allCovid_bracken_phylums.out'),
#         majCovid_bracken = join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_majCovid_bracken_classes.out'),
#         bt2_demix = join(OUT_DIR, 'Freyja', 'Demix', 'Results', '{sample}_freyja.demix'),
#         lcs = join(OUT_DIR, 'LCS', '{sample}', 'outputs', 'decompose', 'lcs.out')
#     params:
#         var_def = VAR_DEF
#     output:
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'kallisto.out'),
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'pieChart_deconvolution.png'),
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'pieChart_freyja.png'),
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'pieChart_k2_allCovid.png'),
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'pieChart_k2_majorCovid.png'),
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'pieChart_kallisto.png'),
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'pieChart_lcs.png')
#     log:
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'qualmap.log')
#     benchmark:
#         join(OUT_DIR, 'SummaryPie', '{sample}', 'qualmap_benchmark.tsv')
#     threads:
#         8
#     resources:
#         mem_mb=32000
#     message:
#         """--- Outputting sumamry variant pie charts from all methods for sample "{wildcards.sample}"."""
#     run:
#         shell(script_dir + 'plotPieChartsforAbundance.py ' + join(OUT_DIR, 'SummaryPie', '{wildcards.sample}') +
#                 ' {params.var_def}'
#                 ' {input.ldvc}'
#                 ' {input.kallisto}'
#                 ' {input.allCovid_bracken}'
#                 ' {input.majCovid_bracken}'
#                 ' {input.bt2_demix}'
#                 ' {input.lcs}'
#                 ' > {log} 2>&1')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule multiQC:
    input:
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '_fastqc.html'), sample = SAMPLES),
        # expand(join(OUT_DIR, 'iVar_trimming_{sample}.log'), sample = SAMPLES),
        # expand(join(OUT_DIR, 'Kallisto', '{sample}', 'abundance.tsv'), sample = SAMPLES),
        # expand(join(OUT_DIR, 'Kraken', '{sample}', '{sample}.k2_std.out'), sample = SAMPLES),
        expand(join(OUT_DIR, 'BAM', '{sample}', '{sample}.resorted.stats'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Pangolin', '{sample}', 'lineage_report.csv'), sample = SAMPLES),
        # expand(join(OUT_DIR, 'LDVC', '{sample}', 'linearDeconvolution_abundance.csv'), sample = SAMPLES),
        expand(join(OUT_DIR, 'FastQs', '{sample}', '{sample}.bam'), sample = SAMPLES),
        expand(join(OUT_DIR, 'BAM', '{sample}', '{sample}.qualimap', 'qualimapReport.html'), sample = SAMPLES)
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
        # shell('ls -1 ' + join(OUT_DIR) + '/*log >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        # shell('ls -1 ' + join(OUT_DIR) + '/Kallisto/*/kallisto.log >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/BAM/*/*sorted.stats >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/*/* | grep ":" | sed "s/://g" >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        # shell('ls -1 ' + join(OUT_DIR) + '/Kraken/* | grep ":" | sed "s/://g" >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))

        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -l -dd 2 ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt') +
                ' > {log} 2>&1')
