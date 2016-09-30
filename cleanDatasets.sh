#!/bin/sh
affymetrixFilename=AffymetrixSmall.txt
head $affymetrixFilename -n 1 > STR_to_append
python cleanDatasets.py MASTERAnonymizedClinicalDataset23Feb2012.txt $affymetrixFilename DMBDI\ genotyping\ blood\ DNA\ Glasgow\ data\ with\ DMBDI\ ID\ for\ John.txt >> STR_to_append
tail $affymetrixFilename -n +2 >> STR_to_append
