#!/bin/env python
from __future__ import print_function
import sys
import json

filename = sys.argv[1]
with open(sys.argv[2]) as f:
    idsMapping = json.load(f)
with open(filename, "r") as f:
    lines = f.readlines()
    swapped = [idsMapping[line] for line in lines[0].rstrip().split("\t")[1:]]
print("IDs\t" + "\t".join(swapped))

