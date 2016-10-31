#!/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib as mpl 
mpl.use('agg')
import matplotlib.pyplot as plt 
import scipy.stats

cleaned = sys.argv[1]

with open(cleaned) as f:
    rowID = None
    skipCount = 0
    while rowID != "Affymetrix":
        line = f.readline().rstrip().split()
        skipCount += 1
        rowID = line[0]
        if rowID == "IDs":
            columns = range(len(line))[1:]
            ids = line[1:]
        if rowID == "ModalAllele":
            modalAllele = [int(i) for i in line[1:]]
            

with open(cleaned) as f:
    dataArray =  np.loadtxt(f, delimiter="\t", skiprows=skipCount, usecols=columns)

residualsArray = np.zeros_like(dataArray)

rowNo, columnNo = dataArray.shape

#print(modalAllele)
#print(dataArray[0])
for i in range(rowNo):
    slope, intercept, _, _, _ = scipy.stats.linregress(modalAllele, list(dataArray[i]))
    expected = [slope * value + intercept for value in modalAllele]
    residuals = [actual - expected for actual, expected in zip(dataArray[i], expected)]
    for j, residual in enumerate(residuals):
        residualsArray[i][j] = residual
    #break
averageResiduals = {}
rejected = 0
for i in range(columnNo):
    averageResiduals[ids[i]] = sum(residualsArray[:,i])/len(residualsArray[:,i])
    if scipy.stats.ttest_1samp(residualsArray[:,i], 0)[1] < 0.05:
        rejected += 1
print(rejected)
import json

# Test results the same with other residuals.
#with open("residuals") as f:
#    otherResiduals = json.load(f)
#for key in otherResiduals:
#    if averageResiduals[key] - otherResiduals[key] > 0.0000001:
#        raise(ValueError("Programmatic Error"))


# Boxplot
data_to_plot = [[dataArray[:, i]] for i in range(columnNo)]
fig = plt.figure(1, figsize=(9, 6))
ax = fig.add_subplot(111)
bp = ax.boxplot(data_to_plot, showfliers=False)
fig.savefig(sys.argv[1] + 'Boxplot.png', bbox_inches='tight')
