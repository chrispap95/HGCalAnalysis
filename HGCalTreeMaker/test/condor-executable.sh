#!/bin/bash

# requires 4 argument inputs:
# 1: UNIQUE_ID - any unique string identifier
# 2: CONDOR_PROCESS - condor process number
# RUN_DIR - running directory (CMSSW_X_Y_Z/subdir)
# mode.  should be 0 for background or 1 for signal

UNIQUE_ID=$1
CONDOR_PROCESS=$2
RUN_DIR=$3
OUT_DIR=$4

START_TIME=`/bin/date`
echo "started at $START_TIME"

#
# setup CMSSW software environment at UMD
#
cd $RUN_DIR

export PATH=/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh


cd /data/users/chpapage/CMSSW_10_6_3_patch1/src
export USERBASE=`pwd`
export ARCH=slc7_amd64_gcc820
export SCRAM_ARCH=slc7_amd64_gcc820

FINAL_PREFIX_NAME=`echo ${UNIQUE_ID}_${CONDOR_PROCESS}`
FINAL_LOG=`echo $FINAL_PREFIX_NAME.log`

cd $RUN_DIR
cmsDriver.py Configuration/GenProduction/python/SingleGammaPt100Eta1p7_pythia8_cfi.py \
    --conditions auto:phase2_realistic -n 100 --era Phase2C8_timing_layer_bar --eventcontent FEVTDEBUG \
    -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC --geometry Extended2023D41 --fileout file:step1_${UNIQUE_ID}_${CONDOR_PROCESS}.root

cmsDriver.py step2  --conditions auto:phase2_realistic \
    -s DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,HLT:@fake2 \
    --datatier GEN-SIM-DIGI-RAW -n 100 --geometry Extended2023D41 --era Phase2C8_timing_layer_bar \
    --eventcontent FEVTDEBUGHLT --filein  file:step1_${UNIQUE_ID}_${CONDOR_PROCESS}.root  --fileout file:step2_${UNIQUE_ID}_${CONDOR_PROCESS}.root

rm -f step1_${UNIQUE_ID}_${CONDOR_PROCESS}.root

cmsDriver.py step3  --conditions auto:phase2_realistic -n 100 --era Phase2C8_timing_layer_bar \
    --eventcontent FEVTDEBUGHLT --runUnscheduled  \
    -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT \
    --datatier GEN-SIM-RECO --geometry Extended2023D41 --filein  file:step2_${UNIQUE_ID}_${CONDOR_PROCESS}.root  --fileout file:step3_${UNIQUE_ID}_${CONDOR_PROCESS}.root

rm -f step2_${UNIQUE_ID}_${CONDOR_PROCESS}.root

cmsRun /data/users/chpapage/CMSSW_10_6_3_patch1/src/full-sim/HGCalAnalysis/HGCalTreeMaker/test/run_HGCalTupleMaker_2023.py inputFiles='file:step3_${UNIQUE_ID}_${CONDOR_PROCESS}.root' outputFile='${UNIQUE_ID}_${CONDOR_PROCESS}.root'

rm -f step3_${UNIQUE_ID}_${CONDOR_PROCESS}.root
mv ${UNIQUE_ID}_${CONDOR_PROCESS}.root ${OUT_DIR}

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
