#!/usr/bin/env python
from __future__ import print_function
import sys
USAGE = """USAGE is:
./cleanDatasets.py <STR_file> <blood_DNA>
Pass in the file with the STR first, should be called something like "MASTERAnonymizedClinicalDataset23Feb2012.txt". Pass in the name of the file with the microarray data, wchich should be called something like "genotyping blood DNA Glasgow...".
"""
try:
    clinicalDataset = sys.argv[1]
except IndexError as e:
    print(USAGE)
    sys.exit(1)
with open(clinicalDataset, "r") as f:
    line = f.readlines()[int(sys.argv[2])]
    for index, entry in enumerate(line.split("\t")):
        print(index, entry)
