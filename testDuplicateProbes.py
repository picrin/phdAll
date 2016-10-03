#!/bin/env python
from __future__ import print_function

# test with smallDuplicateProbes file, which has 3 duplicates.different and 1 duplicates.same
import sys

probes = {}

sameDups = "duplicates.same"
differentDups = "duplicates.different"

print("Puts duplicates in", sameDups, "and", differentDups)

with open(sameDups, "w") as f:
    f.write("")

with open(differentDups, "w") as f:
    f.write("")

filename = sys.argv[1]
with open(filename, "r") as f:
    lines = f.readlines()
    nextLineProbe = False
    for line in lines:
        columns = line.rstrip().split("\t")
        if len(columns) == 0:
            break
        if nextLineProbe:
            if columns[0] not in probes:
                probes[columns[0]] = columns[1:]
            else:
                if probes[columns[0]] == columns[1:]:
                    with open(sameDups, "a") as f:
                        print(columns, file=f)
                else:
                    with open(differentDups, "a") as f:
                        print(columns, file=f)
        if "Affymetrix" in columns[0]:
                nextLineProbe = True
