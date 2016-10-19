#/bin/sh

mkdir -p analyses

ls CELfiles/blood/core; coreExists=$?

set -e
if [[ $coreExists == 0 ]]; then
	echo core file exists.
else
	Rscript coreNormalise.r
fi
cp CELfiles/blood/core analyses

./extractHeaderFromCELfiles.py CELfiles/blood > idsMapping
./swapIDsCELToHUEX.py analyses/core idsMapping > normalisedFirstLine
tail analyses/core -n +2 > normalisedLastLines
cat normalisedFirstLine normalisedLastLines > normalised
