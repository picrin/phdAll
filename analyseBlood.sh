#!/bin/sh
set -e
tail CELfiles/blood/core -n +2 > normalisedLastLines
cat normalisedFirstLine normalisedLastLines > normalised
affymetrixFilename=normalised
echo cleaning the dataset $affymetrixFilename.
head $affymetrixFilename -n 1 > cleaned
python cleanDatasets.py MASTERAnonymizedClinicalDataset23Feb2012.txt $affymetrixFilename DMBDI\ genotyping\ blood\ DNA\ Glasgow\ data\ with\ DMBDI\ ID\ for\ John.txt >> cleaned
echo Affymetrix >> cleaned
tail $affymetrixFilename -n +2 >> cleaned
echo file cleaned written to $PWD
echo file containes $(wc -l cleaned | awk '{print $1}') lines
