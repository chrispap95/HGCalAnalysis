#Instructions for initial setup
First set up a CMSSW release. In this study we need to use geometry D41, so you need at least a version ```CMSSW_10_6_X``` and up. To run that you need a ```slc7_amd64_gcc820``` architecture. You can either find that by logging into ```lxplus7.cern.ch``` or ```siab-1.umd.edu```.

In your work area (say ```/data/users/$USER``` or ```/afs/cern.ch/work/${your_inital}/$USER```) do
```bash
cmsrel CMSSW_10_6_3_patch1
cd CMSSW_10_6_3_patch1/src
cmsenv
git clone https://github.com/chrispap95/HGCalAnalysis.git
scramv1 b -j 8
```
# HGCalAnalysis

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
