#!/bin/sh
./generateData.py > randomData
./computeResiduals.py randomData randomAvgResidual > randomResiduals
