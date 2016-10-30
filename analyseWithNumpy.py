#!/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib as mpl 
mpl.use('agg')
import matplotlib.pyplot as plt 

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

#print(rowID)
#print(columns)
#print(skipCount)
with open(cleaned) as f:
    dataArray =  np.loadtxt(f, delimiter="\t", skiprows=skipCount, usecols=columns)

rowNo, columnNo = dataArray.shape


## combine these different collections into a list    
data_to_plot = [[dataArray[:, i]] for i in range(columnNo)]

# Create a figure instance
fig = plt.figure(1, figsize=(9, 6))

# Create an axes instance
ax = fig.add_subplot(111)

# Create the boxplot
bp = ax.boxplot(data_to_plot, showfliers=False)

# Save the figure
fig.savefig(sys.argv[1] + 'Boxplot.png', bbox_inches='tight')
