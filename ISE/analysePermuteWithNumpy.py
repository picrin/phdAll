#!/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib as mpl 
mpl.use('agg')
import matplotlib.pyplot as plt 
import scipy.stats
import random

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
    random.shuffle(residuals)
    for j, residual in enumerate(residuals):
        residualsArray[i][j] = residual
    #break
averageResiduals = {}
p_values = []
for i in range(columnNo):
    #averageResiduals[ids[i]] = sum(residualsArray[:,i])/len(residualsArray[:,i])
    p_value = scipy.stats.ttest_1samp(residualsArray[:,i], 0)[1]
    p_values.append(p_value)

import json

# Test results the same with other residuals.
#with open("residuals") as f:
#    otherResiduals = json.load(f)
#for key in otherResiduals:
#    if averageResiduals[key] - otherResiduals[key] > 0.0000001:
#        raise(ValueError("Programmatic Error"))


# Boxplot
data_to_plot = [[dataArray[:, i]] for i in range(columnNo)]
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111)
bp = ax.boxplot(data_to_plot, showfliers=False)
fig.savefig(sys.argv[1] + 'PermutedBoxplot.png', bbox_inches='tight')

# p-values
# data_to_plot = []
#print(zip(ids,p_values))

import math
fig = plt.figure(figsize=(9,6))
ax = fig.add_subplot(111)
x = range(len(ids))
#loggedPValues = [math.log(i, 10) for i in p_values]
#loggedPValues = p_values
#ax.semilogx(x, p_values)
ax.semilogy(range(len(p_values)), p_values)
fig.savefig(sys.argv[1] + "PermutedPValues.png", bbox_inches='tight')


