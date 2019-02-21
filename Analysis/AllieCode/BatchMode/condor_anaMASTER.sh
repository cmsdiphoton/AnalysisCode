#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 project CMSSW CMSSW_8_0_24_patch1`
mv anaNNNN.C CMSSW_8_0_24_patch1/src/
mv Analysis.tar.gz CMSSW_8_0_24_patch1/src/
cd CMSSW_8_0_24_patch1/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
echo "CMSSW: "$CMSSW_BASE
tar -xf Analysis.tar.gz 
rm Analysis.tar.gz
mv anaNNNN.C Analysis/macro
cd Analysis/macro
root -q -b anaNNNN.C
xrdcp -f OOOO FINALDESTINATION
rm OOOO
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_8_0_24_patch1