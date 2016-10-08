#/bin/sh
set -xev
#Rscript normalise.r
rm -f self-normalised.txt self-normalised2.txt
cp analyses/$1 self-normalised.txt
head -n 1 self-normalised.txt > self-normalised2.txt
echo Affymetrix >> self-normalised2.txt
tail -n +2 self-normalised.txt >> self-normalised2.txt
./testDuplicateProbes.py self-normalised2.txt 
wc -l duplicates.same
wc -l duplicates.different
