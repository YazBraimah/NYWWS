#!/usr/bin/env python3
'''
Make a samples.json file with sample names and file names.
'''
def msg(name=None):
    return ''' make_bam_samples.py <samples_files>

        '''

import json
from glob import glob
from sys import argv
import sys
import argparse
import glob
parser = argparse.ArgumentParser(description='Make a samples.json file with sample names and file names.', usage=msg())


fqFiles = argv[1:]
bams = []
for fq in fqFiles:
    bams.extend(glob.glob(fq))

bamFILES = {}

# Change this line to extract a sample name from each filename.
SAMPLES = [bam.split('/')[-1].split('.')[0] for bam in bams]

for sample in SAMPLES:
    mate1 = lambda bam: sample in bam and 'trim' in bam
    bamFILES[sample] = {}
    bamFILES[sample]['trim'] = sorted(filter(mate1, bams))



js_bam = json.dumps(bamFILES, indent = 4, sort_keys=True)
open('samples_bam.json', 'w').writelines(js_bam)
