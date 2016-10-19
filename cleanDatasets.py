#!/usr/bin/env python
from __future__ import print_function
import sys
USAGE = """USAGE is:
./cleanDatasets.py <STR_file> <Affymetrix_file> <id_mapping>

STR_file
A file with the STR data, should be called something like "MASTERAnonymizedClinicalDataset23Feb2012.txt". Needs to contain Gene Logic ID

Affymetrix_file
A file with the microarray data, wchich should be called something like "genotyping blood DNA Glasgow...".

id_mapping
A file which maps ids from STR file to ids in Affymetrix_file.
"""


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
    print(USAGE, file=sys.stderr)
STRIDs = {}
patientsWithoutSTR = 0

# Variables relating to the STR_file
firstLineWithDataSTR = 6
columnWithIDIndex = 109
columnWithSTRProgeny = 112
columnWithSTRModal = 113
columnWithSTRNormal1 = 110
columnWithSTRNormal2 = 111

with open(STR_file, "r") as f:
    for i, line in enumerate(f.readlines()[firstLineWithDataSTR:]):
        line = line.split("\t")
        ID = line[columnWithIDIndex]
        progeny = line[columnWithSTRProgeny]
        modal = line[columnWithSTRModal]
        if progeny == '' or modal == '':
            try:
                al1 = int(line[columnWithSTRNormal1])
            except ValueError:
                al1 = None
            try:
                al2 = int(line[columnWithSTRNormal2])
            except ValueError:
                al2 = None
            if al1 is None and al2 is not None:
                STRIDs[ID] = {"progenySTR": al2, "modalSTR": al2}
            if al2 is None and al1 is not None:
                STRIDs[ID] = {"progenySTR": al1, "modalSTR": al1}
            if al1 is None and al2 is None:
                print(i, al1, al2, file=sys.stderr)
                #raise ValueError("Insufficient Data")
            STRIDs[ID] = {"progenySTR": max(al1, al2), "modalSTR": max(al1, al2)}
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
                #if not geneLogicIDToCase[geneLogicID]:
                #    raise(ValueError("should be a case"))
                pass
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
                #if not geneLogicIDToCase[geneLogicID]:
                #    raise(ValueError("should be a case"))
                pass
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
