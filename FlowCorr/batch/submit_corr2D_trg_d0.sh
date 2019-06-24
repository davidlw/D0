#!/bin/bash

echo "setup cmssw"
cd /home/yz144/pPb2016
source setupfile8xy
cd CMSSW_8_0_31/src
eval `scramv1 runtime -sh`
cd D0/FlowCorr/batch
echo PWD: $PWD

./corr2D_trg_d0 list/data.list.$1
