#!/bin/bash

# requires 7 argument inputs:
#
# 1: UNIQUE_ID - any unique string identifier
# 2: CONDOR_PROCESS - condor process number
# 3: OUT_DIR - output directory
# 4: RUN_DIR - running directory (CMSSW_X_Y_Z/subdir)
# 5: PARAMETER_SET - the filename of the template parameter set
#                    (must be in RUN_DIR)
# 6: RANDOM_NUMBER - starting random seed (provide a dummy
#                    number if the job doesn't need a seed)
# 7: NUM_EVENTS_PER_JOB - number of events per job
#
# Creates a new configuration file where key phrases
# in the template PAREMETER_SET file are replaced with values
# modified for each job.  If the key phrase isn't
# present in the template file, no errors ensue,
# the phrase is simply ignored.
#
# CONDOR_RANDOMLOG : Modifies the name of the random
#                    number generator log file, if one
#                    exists.  Naming syntax is of the form
#                    ${OUT_DIR}/${UNIQUE_ID}_${CONDOR_PROCESS}.random
# CONDOR_RANDOMNUMBER : Modifies the random seed for each
#                       job, if needed.  Sets the seed for
#                       each job to:
#                       $RANDOM_NUMBER + 100*$CONDOR_PROCESS
# CONDOR_MAXEVENTS : Set to $NUM_EVENTS_PER_JOB
# CONDOR_SKIPEVENTS : Set to $CONDOR_PROCESS*$NUM_EVENTS_PER_JOB
# CONDOR_OUTPUTFILENAME : ${OUT_DIR}/${UNIQUE_ID}_${CONDOR_PROCESS}.output.root
# CONDOR_HISTOFILENAME : Needed if the job produces two output
#                        files (such as flat histo and ED).
#                        ${OUT_DIR}/${UNIQUE_ID}_${CONDOR_PROCESS}.histo.root
#
# Then executes cmsRun inside $RUN_DIR, sending output to
# $OUT_DIR.

UNIQUE_ID=$1
CONDOR_PROCESS=$2
OUT_DIR=$3
RUN_DIR=$4
PARAMETER_SET=$5
RANDOM_NUMBER=$6
NUM_EVENTS_PER_JOB=$7

#
# header
#

echo ""
echo "CMSSW on Condor"
echo ""

START_TIME=`/bin/date`
echo "started at $START_TIME"

echo ""
echo "parameter set:"
echo "UNIQUE_ID: $UNIQUE_ID"
echo "CONDOR_PROCESS: $CONDOR_PROCESS"
echo "OUT_DIR: $OUT_DIR"
echo "RUN_DIR: $RUN_DIR"
echo "PARAMETER_SET: $PARAMETER_SET"
echo "RANDOM_NUMBER: $RANDOM_NUMBER"
echo "NUM_EVENTS_PER_JOB: $NUM_EVENTS_PER_JOB"

#
# setup CMSSW software environment at UMD
#
export VO_CMS_SW_DIR=/sharesoft/cmssw
. $VO_CMS_SW_DIR/cmsset_default.sh
cd $RUN_DIR
eval `scramv1 runtime -sh`

#
# modify parameter-set
#

FINAL_PREFIX_NAME=`echo ${OUT_DIR}/${UNIQUE_ID}_${CONDOR_PROCESS}`
FINAL_PARAMETER_SET=`echo $FINAL_PREFIX_NAME.py`
RANDOM_LOG=`echo ${UNIQUE_ID}_${CONDOR_PROCESS}.random`
FINAL_LOG=`echo $FINAL_PREFIX_NAME.log`
FINAL_FILENAME=`echo $FINAL_PREFIX_NAME.output.root`
FINAL_HISTONAME=`echo $FINAL_PREFIX_NAME.histo.root`
let "FINAL_RANDOMNUMBER = $CONDOR_PROCESS*100 + $RANDOM_NUMBER"
let "SKIP_EVENTS = $CONDOR_PROCESS*$NUM_EVENTS_PER_JOB"
let "PROCESS_PLUS_ONE = $CONDOR_PROCESS+1"
echo ""
echo "Writing final parameter-set: $FINAL_PARAMETER_SET to OUT_DIR: $OUT_DIR"
echo ""

cat $PARAMETER_SET \
| sed -e s~CONDOR_RANDOMLOG~$RANDOM_LOG~ \
| sed -e s/CONDOR_RANDOMNUMBER/$FINAL_RANDOMNUMBER/ \
| sed -e s/CONDOR_MAXEVENTS/$NUM_EVENTS_PER_JOB/ \
| sed -e s/CONDOR_SKIPEVENTS/$SKIP_EVENTS/ \
| sed -e s~CONDOR_OUTPUTFILENAME~$FINAL_FILENAME~ \
| sed -e s~CONDOR_HISTOFILENAME~$FINAL_HISTONAME~ \
> $FINAL_PARAMETER_SET

#
# run cmssw
#
cmsDriver.py Configuration/GenProduction/python/SingleGammaPt100Eta1p7_pythia8_cfi.py \
    --conditions auto:phase2_realistic -n $NUM_EVENTS_PER_JOB --era Phase2C8_timing_layer_bar --eventcontent FEVTDEBUG \
    -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC --geometry Extended2023D41 --fileout file:step1_${UNIQUE_ID}_${CONDOR_PROCESS}.root

cmsDriver.py step2  --conditions auto:phase2_realistic \
    -s DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,HLT:@fake2 \
    --datatier GEN-SIM-DIGI-RAW -n $NUM_EVENTS_PER_JOB --geometry Extended2023D41 --era Phase2C8_timing_layer_bar \
    --eventcontent FEVTDEBUGHLT --filein  file:step1_${UNIQUE_ID}_${CONDOR_PROCESS}.root  --fileout file:step2_${UNIQUE_ID}_${CONDOR_PROCESS}.root

rm -f step1_${UNIQUE_ID}_${CONDOR_PROCESS}.root

cmsDriver.py step3  --conditions auto:phase2_realistic -n $NUM_EVENTS_PER_JOB --era Phase2C8_timing_layer_bar \
    --eventcontent FEVTDEBUGHLT --runUnscheduled  \
    -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT \
    --datatier GEN-SIM-RECO --geometry Extended2023D41 --filein  file:step2_${UNIQUE_ID}_${CONDOR_PROCESS}.root  --fileout file:step3_${UNIQUE_ID}_${CONDOR_PROCESS}.root

rm -f step2_${UNIQUE_ID}_${CONDOR_PROCESS}.root

cmsRun /data/users/chpapage/CMSSW_10_6_3_patch1/src/full-sim/HGCalAnalysis/HGCalTreeMaker/test/run_HGCalTupleMaker_2023.py inputFiles='file:step3_${UNIQUE_ID}_${CONDOR_PROCESS}.root' outputFile='${UNIQUE_ID}_${CONDOR_PROCESS}.root'

rm -f step3_${UNIQUE_ID}_${CONDOR_PROCESS}.root
mv ${UNIQUE_ID}_${CONDOR_PROCESS}.root ${OUT_DIR}

exitcode=$?

# The RandomNumberServiceHelper wants the log placed in the running dir.
# Transfer the random # log to the out_dir once job completes.
echo "Moving random log, if it exists.  Ignore error if file does not exist, as not all jobs have random seed logs."
mv $RUN_DIR/$RANDOM_LOG $OUT_DIR/$RANDOM_LOG

#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
