#!/bin/env python
from __future__ import print_function

import scipy.stats
import numpy
import random

norm = scipy.stats.norm

genes = 20000

individuals = 35

genLoc = 10
genSca = 3

samSca = 2

print("IDs", end="\t")
for i in range(individuals):
    print(i, end="\t")
print()

print("ModalAllele", end="\t")
for i in range(individuals):
    modalAlleleLength = random.randint(0, 1300)
    print(modalAlleleLength, end="\t")
print()

print("Affymetrix")
for i in range(genes):
    print(random.randint(20000000, 30000000), end="\t")
    samLoc = norm.rvs(size=1, loc=genLoc, scale=genSca)[0]
    genes = norm.rvs(size=individuals, loc=samLoc, scale=samSca)
    for gene in genes:
        print(gene, end="\t")
    print()
