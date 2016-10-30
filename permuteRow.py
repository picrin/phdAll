#!/bin/env python
from __future__ import print_function
import random
import sys

modalAllele = sys.stdin.readlines()[0]
modalAllele = modalAllele.rstrip()
modalAllele = modalAllele.split("\t")
head = modalAllele[0]
tail = modalAllele[1:]
random.shuffle(tail)

print(head, end="\t")
for elem in tail:
    print(elem, end="\t")
print("")
