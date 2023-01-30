"""
Author: Y. Ahmed-Braimah
--- NY Waste Water Surveillance pipeline
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
BAMS_DIR = config['BAMS_DIR']
SAM_HEADER = config['SAM_HEADER']
SAMPLE_INFO = config['SAMPLE_INFO']
covidRefSequences = config['covidRefSequences']
r_script = config['R_script']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
HOME_DIR = config['HOME_DIR']

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
        join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_parse.csv'),
        join(OUT_DIR, 'Coverage', 'coverageReport.txt'),
        expand(join(OUT_DIR, 'iVar', '{sample}.rawVarCalls.tsv'), sample = SAMPLES),
        join(OUT_DIR, 'Summary', 'variant_tables', 'varTables_ok')

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
        shell(HOME_DIR + 'scripts/sam2fastq.py'
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

rule Pangolin_variants:
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
    conda:
        "envs/pangolin_env.yml"
    shell:
        '/home/yahmed/miniconda3/envs/pangolin/bin/pangolin'
                ' {input.consensus}'
                ' --outfile {output}'
                ' --threads 4'
                ' > {log} 2>&1'

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

rule Freyja_parse:
    input:
        bt2_agg = join(OUT_DIR, 'Freyja', 'Aggregate', 'aggregated_results.tsv')
    output:
        parse = join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_parse.csv')
    threads:
        2
    resources:
        mem_mb=4000
    message:
        """--- Parsing Freyija output. """
    run:
        shell('python ' + join(HOME_DIR, 'scripts', 'parse_demix.py') + ' {input.bt2_agg} > {output.parse}')

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

rule Coverage_summary:
    input:
        expand(join(OUT_DIR, 'BAM', '{sample}', '{sample}.qualimap', 'qualimapReport.html'), sample = SAMPLES)
    output:
        join(OUT_DIR, 'Coverage', 'coverageReport.txt')
    log:
        join(OUT_DIR, 'Coverage', 'cov.log')
    benchmark:
        join(OUT_DIR, 'Coverage', 'cov_benchmark.tsv')
    threads:
        2
    resources:
        mem_mb=8000
    message:
        """--- Outputting genome wide coverage stats."""
    run:
        shell('cd ' + join(OUT_DIR, 'BAM') + ' &&'
                ' ls -1 | for file in `cat -`; \
                do sed 1d $file/$file.qualimap/raw_data_qualimapReport/coverage_across_reference.txt | \
                sed "s/^/$file\t/g" >> {output} ;\
                done')

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule report_summary:
    input:
        cov = join(OUT_DIR, 'Coverage', 'coverageReport.txt'),
        parse = join(OUT_DIR, 'Freyja', 'Aggregate', 'freyja_parse.csv')
    params:
        samInfo = SAMPLE_INFO,
        r_script = r_script
    output:
        coveragePointPlot = join(OUT_DIR, 'Summary', 'coveragePointPlot_by_date.pdf'),
        coverageCovClass = join(OUT_DIR, 'Summary', 'coveragePointPlot_by_covClass.pdf'),
        varTables = join(OUT_DIR, 'Summary', 'variant_tables', 'varTables_ok'),
        varFreqBarPlotsAll = join(OUT_DIR, 'Summary', 'Variant_frequency_barplots_by_county_all_samples.pdf'),
        varFreqBarPlotsMin20X = join(OUT_DIR, 'Summary', 'Variant_frequency_barplots_by_county_min_20X.pdf')
    log:
        join(OUT_DIR, 'Summary', 'sum.log')
    benchmark:
        join(OUT_DIR, 'Summary', 'sum_benchmark.tsv')
    threads:
        8
    resources:
        mem_mb=32000
    message:
        """--- Outputting summary tables and plots."""
    conda:
        'envs/r_env.yml'
    shell:
        'Rscript {params.r_script} {input.parse} {input.cov} {params.samInfo} {output.coveragePointPlot} {output.coverageCovClass} {output.varFreqBarPlotsAll} {output.varFreqBarPlotsMin20X} && '
        'mv *_County_variant_table.pdf ' + join(OUT_DIR, 'Summary', 'variant_tables') +
        ' && touch {output.varTables}'

##--------------------------------------------------------------------------------------##
##--------------------------------------------------------------------------------------##

rule multiQC:
    input:
        expand(join(OUT_DIR, 'fastQC', '{sample}' + '_fastqc.html'), sample = SAMPLES),
        expand(join(OUT_DIR, 'BAM', '{sample}', '{sample}.resorted.stats'), sample = SAMPLES),
        expand(join(OUT_DIR, 'Pangolin', '{sample}', 'lineage_report.csv'), sample = SAMPLES),
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
        shell('ls -1 ' + join(OUT_DIR) + '/BAM/*/*sorted.stats >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/*/* | grep ":" | sed "s/://g" >> ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt'))

        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -l -dd 2 ' + join(OUT_DIR, 'MultiQC', 'summary_files.txt') +
                ' > {log} 2>&1')
