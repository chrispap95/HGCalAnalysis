#!/bin/tcsh
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
setenv PATH /usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
printenv
source /cvmfs/cms.cern.ch/cmsset_default.csh  ## if a bash script, use .sh instead of .csh
echo "source done"
#
# TO BE CHECKED before submission
#
setenv MYCMSSW CMSSW_10_4_0_pre2 ## <========= TO BE CHECKED
setenv VERSION v01                        ## <========= TO BE CHECKED
#
#
echo ${MYCMSSW} ${VERSION}
#
printf "Start time: "; /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd
voms-proxy-info
voms-proxy-info -fqan
#
### for case 1. EOS have the following line, otherwise remove this line in case 2.
xrdcp -s root://kodiak-se.baylor.edu//store/user/bcaraway/condor/tarballs/${MYCMSSW}_HGCAL_condor.tgz .
tar -xf ${MYCMSSW}_HGCAL_condor.tgz
rm ${MYCMSSW}_HGCAL_condor.tgz
#setenv SCRAM_ARCH slc6_amd64_gcc530
ls -R
cd ${MYCMSSW}/src
scramv1 b ProjectRename
eval `scramv1 runtime -csh` # cmsenv is an alias not on the workers
cmsRun ../../run_HGCalTupleMaker_2023_PF_full.py maxEvents=1000 skipEvents=`echo ${1}\*1000|bc`
foreach f (`ls *10_4_D28*.root`) # # <=== to be checked
   echo $f
   set name=`basename $f .root`
   echo $name
   gfal-copy --just-copy ${f} gsiftp://kodiak-se.baylor.edu/cms/data/store/user/bcaraway/condor/outputs/${name}_${MYCMSSW}_${1}_${VERSION}.root
end
#gfal-copy --just-copy relval_minbias_2018_MCfull.root gsiftp://kodiak-se.baylor.edu/cms/data/store/user/hatake/condor/outputs/relval_minbias_2018_MCfull_${MYCMSSW}_${1}_${VERSION}.root
### remove the output file if you don't want it automatically transferred when the job ends
### rm relval_minbias_2018_MCfull.root
cd ${_CONDOR_SCRATCH_DIR}
rm -rf ${MYCMSSW}
