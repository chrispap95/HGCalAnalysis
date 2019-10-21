
# HGCAL tuple maker
## Create Ntuples (UMD)
To create ntuples for the samples you created, do:
```bash
cmsRun run_HGCalTupleMaker_2023.py inputFiles='file:../../../29005.0_SingleGammaPt35+SingleGammaPt35_pythia8_2023D41_GenSimHLBeamSpotFull+DigiFullTrigger_2023D41+RecoFullGlobal_2023D41+HARVESTFullGlobal_2023D41/step3.root' outputFile='out.root'
```
Be sure to use the correct directory for the step3.root file.
The you can play with the ```rechitSum.cpp``` script to plot the HGCRecHit sum.

## Normal mode (BaylorU)
```
cmsRun run_HGCalTupleMaker_2023.py
cmsRun run_HGCalTupleMaker_2023.py inputFiles='file:../../../24034.0_TTbar_14TeV+TTbar_14TeV_TuneCUETP8M1_2023D28_GenSimHLBeamSpotFull14+DigiFullTrigger_2023D28+RecoFullGlobal_2023D28+HARVESTFullGlobal_2023D28/step3.root' outputFile='ntuples_TTbar.root'
```

### Full(?) content mode (including uncalibrated rechits, digis)
```
# run step3 with cmsDriver.py xxxx ... --outputCommands 'keep *_HGCalUncalibRecHit_*_*','keep *_HGCDigis*_*_*','keep *_hgcalDigis*_*_*'
# and then
cmsRun run_HGCalTupleMaker_2023_full.py
cmsRun run_HGCalTupleMaker_2023_full.py inputFiles='file:../../../24034.0_run1/step3.root' outputFile='ntuples_full_TTbar.root' >& log_full.log &
```

- - - -

## HGCAL sequence:

Kevin's March 1, 2018 talk:
Link: https://indico.cern.ch/event/709595/contributions/2916850/

### SimHits:
```
HGCHitsEE
HGCHitsHEfront
HcalHits
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsEE"       "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsHEfront"   "SIM"     
vector<PCaloHit>                      "g4SimHits"                 "HGCHitsHEback"   "SIM"     
#vector<PCaloHit>                      "g4SimHits"                 "HcalHits"        "SIM"
```

### Digis:
* https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCDigitizer.cc
* https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/python/hgcalDigitizer_cfi.py
* Additional:
* HGCEE: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCEEDigitizer.cc
* HGCHEFront: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCHEfrontDigitizer.cc
* HGC Si digitizer base: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCDigitizerBase.cc
* HGC FE electronics for Si: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCFEElectronics.cc
* HGCHEBack: https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HGCalSimProducers/src/HGCHEbackDigitizer.cc
```
HGCDigisEE
HGCDigisHEfront
HGCDigisHEback
edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >    "mix"                       "HGCDigisEE"      "HLT"     
edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >    "mix"                       "HGCDigisHEback"   "HLT"     
edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >    "mix"                       "HGCDigisHEfront"   "HLT"
and
edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >    "hgcalDigis"                "EE"              "RECO"    
edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >    "hgcalDigis"                "HEback"          "RECO"    
edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >    "hgcalDigis"                "HEfront"         "RECO"    
```
```
Silicon part utilizes HGCFEElectronics::runShaperWithToT
Scinti part utilizes HGCFEElectronics::runTrivialShaper
HEFront:
(in progress...)
HEBack:  
keV->MIP->#p.e.(poisson smeared)->(xtalk&saturation?)->MIP->(noise added)->chargeColl[i](MIPs)->ADC (10bit)
```
<!---
ADC (16 bits? 12 bits? 10 bits?)
--->

### RecHits (uncalib):
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitProducer.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/python/HGCalUncalibRecHit_cfi.py (fCPerMIP)
* Additional:
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitWorkerWeights.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecAlgos/interface/HGCalUncalibRecHitRecWeightsAlgo.h
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecAlgos/interface/HGCalRecHitSimpleAlgo.h
```
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCEEUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEFUncalibRecHits"   "RECO"    
edm::SortedCollection<HGCUncalibratedRecHit,edm::StrictWeakOrdering<HGCUncalibratedRecHit> >    "HGCalUncalibRecHit"        "HGCHEBUncalibRecHits"   "RECO"    
```
```
CEE & HEF: ADC * adcLSB (fC/ADC) / fCPerMIP_[thickness-1] -> amplitude (# of MIP equivalent)  
BH:        ADC * adcLSB (MIP/ADC)                         -> amplitude (# of MIP equivalent)  
adcLSB=0.10, fC/MIP=1.25, 2.57, 3.88 (120, 200, 300 um Si sensors)  for HGCEE & HGCHEF
adcLSB=0.25 for HGCHEB  
```

### RedHits (calib):
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalRecHitProducer.cc (weights_[layer],rcorr_[thickness] etc)
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py
* Additional:
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecProducers/plugins/HGCalRecHitWorkerSimple.cc
* https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/HGCalRecAlgos/interface/HGCalRecHitSimpleAlgo.h
```
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCEERecHits"    "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEFRecHits"   "RECO"    
edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHEBRecHits"   "RECO"    
```
```
uncalibRH.amplitude() * weights_[layer](dEdx, MeV, corrected for energy loss in passive materials) * 0.001 (MeV->GeV)  
 rcorr_[thickness] (0.88, 0.92, or 1) * cce_correction (GeV)
```
