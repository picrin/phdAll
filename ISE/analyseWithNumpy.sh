set -xe
./analyseWithNumpy.py cleaned
./generatePermuted.sh
./analyseWithNumpy.py permuted
./generateRandomData.py > randomData
./analyseWithNumpy.py randomData
./generateRandomTranscriptionFactors.py > randomTranscriptionFactors
./analyseWithNumpy.py randomTranscriptionFactors
./generateIndividualTranscriptionFactors.py > randomIndividualTranscriptionFactors
./analyseWithNumpy.py randomIndividualTranscriptionFactors
./analysePermuteWithNumpy.py cleaned
