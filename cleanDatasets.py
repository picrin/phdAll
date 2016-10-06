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
columnWithMuscleID = 14
columnWithGeneLogicID = 2
columnWithCase = 16

bloodIDToGeneLogic = {}
muscleIDToGeneLogic = {}
geneLogicIDToCase = {}

with open(id_mapping, "r") as f:
    for line in f.readlines()[firstLineWithDataIDMap:]:
        line = line.split("\t")
        geneLogicID = line[columnWithGeneLogicID]
        bloodID = line[columnWithBloodID]
        muscleID = line[columnWithMuscleID]
        bloodIDToGeneLogic[bloodID] = geneLogicID
        muscleIDToGeneLogic[muscleID] = geneLogicID
        case = line[columnWithCase].rstrip()
        if case == "DM1":
            case = True
        elif case == "Control":
            case = False
        geneLogicIDToCase[geneLogicID] = case
AffymetrixIDs = None
with open(Affymetrix_file, "r") as f:
    lines = f.readlines()
    AffymetrixIDs = [id for id in lines[0].rstrip().split("\t")[1:]]

for X in [["ModalAllele", "modalSTR"],["\nEstProgAllele", "progenySTR"]]:
    print(X[0], end="\t")
    for bloodID in AffymetrixIDs:
        if bloodID in bloodIDToGeneLogic:
            geneLogicID = bloodIDToGeneLogic[bloodID]
            try:
                print(STRIDs[geneLogicID][X[1]], end="\t")
            except KeyError:
                if not geneLogicIDToCase[geneLogicID]:
                    print(0, end="\t")
                else:
                    raise(ValueError("should be a control"))
            else:
                if not geneLogicIDToCase[geneLogicID]:
                    raise(ValueError("should be a case"))
        elif bloodID in muscleIDToGeneLogic:
            geneLogicID = muscleIDToGeneLogic[bloodID]
            try:
                print(STRIDs[geneLogicID][X[1]], end="\t")
            except KeyError:
                if not geneLogicIDToCase[geneLogicID]:
                    print(0, end="\t")
                else:
                    print("missing data", end="\t")
            else:
                if not geneLogicIDToCase[geneLogicID]:
                    raise(ValueError("should be a case"))
        else:
            raise(ValueError("insufficientData"))
print("")

print("BloodOrMuscle", end="\t")
for bloodID in AffymetrixIDs:
    if bloodID in bloodIDToGeneLogic:
        print("blood", end="\t")
    elif bloodID in muscleIDToGeneLogic:
        print("muscle", end="\t")
    else:
        raise(ValueError("insufficientData"))
print("")

print("CaseOrControl", end="\t")
for bloodID in AffymetrixIDs:
    if bloodID in bloodIDToGeneLogic:
        geneLogicID = bloodIDToGeneLogic[bloodID]
    elif bloodID in muscleIDToGeneLogic:
        geneLogicID = muscleIDToGeneLogic[bloodID]
    else:
        raise(ValueError("insufficientData"))
    if geneLogicIDToCase[geneLogicID]:
        print("case", end="\t")
    else:
        print("control", end="\t")
print("")

print("1 DM1 without modal data", file=sys.stderr)
