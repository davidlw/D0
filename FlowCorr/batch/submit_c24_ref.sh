#!/bin/bash

echo "setup cmssw"
cd /home/wl33
source setupfile8xy
eval `scramv1 runtime -sh`
cd /home/wl33/D0/FlowCorr/batch
echo PWD: $PWD

../bin/c24_ref prefix$1
