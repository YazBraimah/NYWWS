#!/usr/bin/env python3
'''
Make a samples.json file with sample names and file names.
'''
def msg(name=None):                                                            
    return ''' make_json_samples.py <samples_files>

    fastq file names should have the following format:
        paired-end: <sample_name>.R1.fastq.gz
                    <sample_name>.R2.fastq.gz

        single-end: <sample_name>.R1.fastq.gz
        '''

import json
from glob import glob
from sys import argv
import sys
import argparse
import glob
parser = argparse.ArgumentParser(description='Make a samples.json file with sample names and file names.', usage=msg())


fqFiles = argv[1:]
fastqs = []
for fq in fqFiles:
    fastqs.extend(glob.glob(fq))

peFILES = {}
seFILES = {}

# Change this line to extract a sample name from each filename.
SAMPLES = [fastq.split('/')[-1].split('.')[0] for fastq in fastqs]

for sample in SAMPLES:
    mate1 = lambda fastq: sample in fastq and 'R1' in fastq
    mate2 = lambda fastq: sample in fastq and 'R2' in fastq
    if any('R2' in s for s in sorted(filter(mate2, fastqs))):
       
        peFILES[sample] = {}
        peFILES[sample]['R1'] = sorted(filter(mate1, fastqs))
        peFILES[sample]['R2'] = sorted(filter(mate2, fastqs))
    else:
        
        seFILES[sample] = {}
        seFILES[sample]['R1'] = sorted(filter(mate1, fastqs))



js_pe = json.dumps(peFILES, indent = 4, sort_keys=True)
js_se = json.dumps(seFILES, indent = 4, sort_keys=True)
open('samples_pe.json', 'w').writelines(js_pe)
open('samples_se.json', 'w').writelines(js_se)