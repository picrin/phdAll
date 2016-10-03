#!/bin/sh
set -xe
affymetrixFilename=$1
echo !!!
echo pass in a filename containing Affymetrix data you want analysed as the first argument.
echo !!!
echo cleaning the dataset $affymetrixFilename.
head $affymetrixFilename -n 1 > cleaned
python cleanDatasets.py MASTERAnonymizedClinicalDataset23Feb2012.txt $affymetrixFilename DMBDI\ genotyping\ blood\ DNA\ Glasgow\ data\ with\ DMBDI\ ID\ for\ John.txt >> cleaned
echo "Affymetrix" >> cleaned
tail $affymetrixFilename -n +2 >> cleaned
./testDuplicateProbes.py cleaned
echo repeated probes
wc -l duplicates.same
echo conflicting probes
wc -l duplicates.different
./makeLinearRegression.py cleaned
