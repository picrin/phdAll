#!/usr/bin/env python
from __future__ import print_function
import sys
USAGE = """USAGE is:
./cleanDatasets.py <STR_file> <Affymetrix_file>

STR_file
A file with the STR data, should be called something like "MASTERAnonymizedClinicalDataset23Feb2012.txt". Needs to contain Gene Logic ID

Affymetrix_file
A file with the microarray data, wchich should be called something like "genotyping blood DNA Glasgow...".

id_mapping
A file which maps ids from STR file to ids in Affymetrix_file.
"""

# Variables relating to the STR_file
firstLineWithDataSTR = 6
columnWithIDIndex = 109
columnWithSTRProgeny = 112
columnWithSTRModal = 113


def countDuplicates(list1):
    set1 = set(list1)
    return len(list1) - len(set1)

def countIncludedExcluded(froma, to):
    fromSet = set(froma)
    toSet = set(to)
    included = fromSet.intersection(toSet)
    excluded = len(fromSet) - included
    return included, excluded

try:
    STR_file = sys.argv[1]
    Affymetrix_file = sys.argv[2]
    id_mapping = sys.argv[3]
except IndexError as e:
    print(USAGE)
    sys.exit(1)

STRIDs = {}
patientsWithoutSTR = 0

with open(STR_file, "r") as f:
    for line in f.readlines()[firstLineWithDataSTR:]:
        line = line.split("\t")
        ID = line[columnWithIDIndex]
        progeny = line[columnWithSTRProgeny]
        modal = line[columnWithSTRModal]
        # If a patient is missing this, then they are a control, even thought they're
        # included in this file (sic)!
        if progeny == '' or modal == '':
            patientsWithoutSTR += 1
        else:
            STRIDs[ID] = {"progenySTR": progeny, "modalSTR": modal}

# Variables relating to the id_mapping file
firstLineWithDataIDMap = 3
columnWithBloodID = 13
columnWithGeneLogicID = 2

bloodIDToGeneLogic = {}

with open(id_mapping, "r") as f:
    for line in f.readlines()[firstLineWithDataIDMap:]:
        line = line.split("\t")
        geneLogicID = line[columnWithGeneLogicID]
        bloodID = line[columnWithBloodID]
        bloodIDToGeneLogic[bloodID] = geneLogicID

AffymetrixIDs = None
with open(Affymetrix_file, "r") as f:
    lines = f.readlines()
    AffymetrixIDs = [id for id in lines[0].rstrip().split("\t")[1:]]


controls = 0
DM1 = 0
print("ModalAllele", end="\t")
for bloodID in AffymetrixIDs:
        try:
            geneLogicID = bloodIDToGeneLogic[bloodID]
            print(STRIDs[geneLogicID]["modalSTR"], end="\t")
        except KeyError:
            if bloodID in ["189821HUEX1A11"]:
                print(-1, end="\t")
            else:
                print(0, end="\t")

print("\nEstProgAllele", end="\t")
for bloodID in AffymetrixIDs:
        try:
            geneLogicID = bloodIDToGeneLogic[bloodID]
            print(STRIDs[geneLogicID]["progenySTR"], end="\t")
            DM1 += 1
        except KeyError:
            if bloodID in ["189821HUEX1A11"]:
                print(-1, end="\t")
            else:
                print(0, end="\t")
                controls += 1
import sys
print(DM1, controls, file=sys.stderr)
#print(len(controls), "controls", controls)
#print(len(DM1), "DM1", DM1)

#print("Patients with STR measured", len(STRIDs))
#print("Patients without STR data:", patientsWithoutSTR)
#print("Duplicate STRIDs: ", countDuplicates(STRIDs))
#print("zip(STRIDs, STRProgeny, STRModal)")
#print(zip(STRIDs, STRProgeny, STRModal))
