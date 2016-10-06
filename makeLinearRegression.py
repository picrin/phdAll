#!/bin/env python
import scipy, scipy.stats
import sys
import matplotlib.pyplot as plt
datafilename = sys.argv[1]

def lowerUpperBoundary(varsy, slackFraction):
    varsyMax = max(varsy)
    varsyMin = min(varsy)
    slack = (varsyMax - varsyMin) * slackFraction
    return varsyMin - slack, varsyMax + slack

def bestFitLine(leftBorder, rightBorder, intercept, slope):
    lowY = leftBorder * slope + intercept
    highY = rightBorder * slope + intercept
    return [leftBorder, rightBorder], [lowY, highY]

def makePlot(independentVar, dependentVar, independentVarName, dependentVarName, plotname):
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(independentVar, dependentVar)

    leftBorder, rightBorder = lowerUpperBoundary(independentVar, 0.1)

    leftBorderLine, rightBorderLine = lowerUpperBoundary(independentVar, 0.05)

    bottomBorder, topBorder = lowerUpperBoundary(dependentVar, 0.3)

    xs, ys = bestFitLine(leftBorderLine, rightBorderLine, intercept, slope)



    fig = plt.figure(figsize=(10, 5))

    verticalBar = scipy.std(dependentVar)

    horizontalBar = scipy.std(independentVar)    

    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.15, bottom=0.4, top=0.85)
    #ax.set_title(descriptions[figureNo] + '\n', fontsize=16)
    #figureNo += 1
    d1 = 'Best fit line slope: ' + '%.3g' % slope + "\n"
    d2 = "Pearson's coefficient-value: " + '%.3g' % r_value + "\n"
    d3 = 'Null hypothesis p-value: ' + '%.3g' % p_value + "\n"
    
    fig.text(0.12, 0.02, d1 + d2 + d3 , fontsize=16)
    
    ax.set_xlabel(independentVarName, fontsize = 18)
    ax.set_ylabel(dependentVarName, fontsize = 18)
    
    ax.plot(xs, ys, label='line of best fit') # shows line of best fit
    #fig.suptitle('sample correlation', fontsize=12)
    ax.errorbar(independentVar, dependentVar, fmt='go', label='individuals')
    ax.legend(loc=(0.74, -0.72), shadow=True, fancybox=True)
    ax.axis([leftBorder, rightBorder, bottomBorder, topBorder]) # creates canvas

    plt.savefig(plotname)


with open(datafilename, "r") as f:
    lines = f.readlines()
    nextLineProbe = False
    for index, line in enumerate(lines):
        columns = line.rstrip().split("\t")
        if len(columns) == 0:
            break
        if "ModalAllele" in columns[0]:
            independentVar = [int(i) for i in columns[1:]]
        if index == 4:
            dependentVar = [float(i) for i in columns[1:]]
        if "Affymetrix" in columns[0]:
            nextLineProbe = True

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(independentVar, dependentVar)
print(slope, intercept, r_value, p_value, std_err)
makePlot(independentVar, dependentVar, "ModalAllele", "first probe", "firstProbePlot")
