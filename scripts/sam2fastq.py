#! /bin/env python

# A custom script to obtain a fastq file from a sam file via CIGAR strings
# This is useful when there is a bam file with soft clipping
# that needs to be converted into a fastq file
# Soft clipped primers will be hard clipped in the output fastq.

import csv
import re
import gzip


samfilename = snakemake.input[0]
fastqfilename = snakemake.output[0]


with open(samfilename, 'r') as infile, gzip.open(fastqfilename, 'wb') as outfile:
    reader = csv.reader(infile, delimiter="\t", quoting=csv.QUOTE_NONE)
    read_counter=1
    for row in reader:
        if row[0][0] != '@':
            CIGAR = row[5]
            seq = row[9]
            quality = row[10]
            
            # Process the CIGAR string
            CIGAR_letters = [x for x in re.split('\d', CIGAR) if x]
            CIGAR_numbers = [x for x in re.split('\D', CIGAR) if x]
            
            if CIGAR_letters[0] == 'S':
                start_pos = int(CIGAR_numbers[0])
            else:
                start_pos = 0
            
            if CIGAR_letters[-1] == 'S':
                end_pos = -int(CIGAR_numbers[-1])
            else:
                end_pos = len(seq)
    
            trimmed_seq = seq[start_pos:end_pos]
            trimmed_quality = quality[start_pos:end_pos]
            
            # Output to the fastq file
            lines = ["@read%d\n" % read_counter, trimmed_seq, "\n+\n", trimmed_quality, "\n"]
            outfile.writelines([line.encode("utf-8") for line in lines])
            read_counter += 1

