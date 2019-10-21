# HGCalAnalysis
## Instructions for initial setup (UMD)
First set up a CMSSW release. In this study we need to use geometry D41, so you need at least a version ```CMSSW_10_6_X``` and up. To run that you need a ```slc7_amd64_gcc820``` architecture. You can either find that by logging into ```lxplus7.cern.ch``` or ```siab-1.umd.edu```.

In your work area (say ```/data/users/$USER``` or ```/afs/cern.ch/work/${your_inital}/$USER```) do
```bash
cmsrel CMSSW_10_6_3_patch1
cd CMSSW_10_6_3_patch1/src
cmsenv
git clone https://github.com/chrispap95/HGCalAnalysis.git
scramv1 b -j 8
```

## Instructions for creating samples (UMD)
You can find correct workflows by using ```runTheMatrix.py```. Specifically, doing
```bash
runTheMatrix.py -w upgrade -n
```
is going to give all possible predetermined workflows. You can grep on them using ```2023D41_``` to look for the correct geometry and no PU.
One interesting example is 29005.0 (SingleGammaPt35). To get the workflow:
```bash
runTheMatrix.py -w upgrade -l 29005.0 --command="--no_exec"
```
Then you can execute step by step the commands listed under ```29005.0_SingleGammaPt35+SingleGammaPt35_pythia8_2023D41_GenSimHLBeamSpotFull+DigiFullTrigger_2023D41+RecoFullGlobal_2023D41+HARVESTFullGlobal_2023D41/cmdLog```
```bash
cmsDriver.py SingleGammaPt35_pythia8_cfi  --conditions auto:phase2_realistic -n 100 \
    --era Phase2C8_timing_layer_bar --eventcontent FEVTDEBUG --relval 9000,50 -s GEN,SIM \
    --datatier GEN-SIM --beamspot HLLHC --geometry Extended2023D41 --fileout file:step1.root  > \
    step1_SingleGammaPt35+SingleGammaPt35_pythia8_2023D41_GenSimHLBeamSpotFull+DigiFullTrigger_2023D41+RecoFullGlobal_2023D41+HARVESTFullGlobal_2023D41.log  2>&1

cmsDriver.py step2  --conditions auto:phase2_realistic \
    -s DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,HLT:@fake2 \
    --datatier GEN-SIM-DIGI-RAW -n 100 --geometry Extended2023D41 --era Phase2C8_timing_layer_bar \
    --eventcontent FEVTDEBUGHLT --filein  file:step1.root  --fileout file:step2.root  > \
    step2_SingleGammaPt35+SingleGammaPt35_pythia8_2023D41_GenSimHLBeamSpotFull+DigiFullTrigger_2023D41+RecoFullGlobal_2023D41+HARVESTFullGlobal_2023D41.log  2>&1

cmsDriver.py step3  --conditions auto:phase2_realistic -n 100 --era Phase2C8_timing_layer_bar \
    --eventcontent FEVTDEBUGHLT --runUnscheduled  \
    -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM \
    --datatier GEN-SIM-RECO --geometry Extended2023D41 --filein  file:step2.root  --fileout file:step3.root  > \
    step3_SingleGammaPt35+SingleGammaPt35_pythia8_2023D41_GenSimHLBeamSpotFull+DigiFullTrigger_2023D41+RecoFullGlobal_2023D41+HARVESTFullGlobal_2023D41.log  2>&1
```
This will create 100 events under the file ```step3.root```. Now, you can go to ```HGCalTreeMaker/test``` to run the ntuplizer on the samples you created.


## Instructions from BaylorU

```
git clone git@github.com:BaylorCMS/HGCalAnalysis.git
scramv1 b -j 8
cd HGCalAnalysis/HGCalTreeMaker/test
```

- - - -

## Sample locations

### Ntuples
Stored in /cms/data/store/user/hatake/HGCAL/

```
Parent datasets
/RelValSinglePiPt25Eta1p7_2p7/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/GEN-SIM-RECO
/SinglePiPt*Eta1p6_2p8/PhaseIITDRFall17*93X_upgrade2023_realistic_v2*/GEN-SIM-RECO

Available ntuples (examples)
root://kodiak-se.baylor.edu//store/user/hatake/HGCAL/ntuples/93x/results_pt5_v03.root
root://kodiak-se.baylor.edu//store/user/hatake/HGCAL/ntuples/93x/results_pt10_v03.root
root://kodiak-se.baylor.edu//store/user/hatake/HGCAL/ntuples/93x/results_pt15_v03.root
root://kodiak-se.baylor.edu//store/user/hatake/HGCAL/ntuples/93x/results_pt25_v03.root
root://kodiak-se.baylor.edu//store/user/hatake/HGCAL/ntuples/93x/results_pt50_v03.root
root://kodiak-se.baylor.edu//store/user/hatake/HGCAL/ntuples/93x/results_pt100_v03.root
```

### Special reco & ntuples files with HGCDigis and uncalibrated rechits included

```
Reco files:
/cms/data/store/user/hatake/HGCAL/fulreco
e.g.
/cms/data/store/user/hatake/HGCAL/fulreco/step3_pt25.root
/cms/data/store/user/hatake/HGCAL/fulreco/step3_pt25001.root
/cms/data/store/user/hatake/HGCAL/fulreco/step3_pt25002.root

Ntuple files:
/cms/data/store/user/hatake/HGCAL/ntuples/93x
e.g.
/cms/data/store/user/hatake/HGCAL/ntuples/93x/ntuples_digi_pt25_numEvent10.root
```
