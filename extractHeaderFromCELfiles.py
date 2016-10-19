#!/bin/env python
from __future__ import print_function

import sys
import os, os.path
import re
import json

CELfiles = os.listdir(sys.argv[1])

def processCELFile(basePath, lastPath):
    relPath = os.path.join(basePath, lastPath)
    result = ""
    with open(relPath, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines[0:20]):
            try:
                result = re.search("(\d)*HUEX\dA(\d)*", line).group()
            except AttributeError:
                pass
    return result
returna = {}
for celfile in CELfiles:
    returna[celfile] = processCELFile(sys.argv[1], celfile)
print(json.dumps(returna))
