#ifndef fReader_h
#define fReader_h

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TInterpreter.h"
#include "TStyle.h"
#include "TLorentzVector.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

namespace globalTChain
{
 
  TTreeReader  fReader;  //!the tree reader
 
  //
  // Set up TTreeReader's
  // -- use MakeSelector of root
  //
  // Readers to access the data (delete the ones you do not need).
  
  TTreeReaderArray<double> GenParEta = {fReader, "GenParEta"};
  TTreeReaderArray<double> GenParM = {fReader, "GenParM"};
  TTreeReaderArray<double> GenParPhi = {fReader, "GenParPhi"};
  TTreeReaderArray<double> GenParPt = {fReader, "GenParPt"};
  TTreeReaderArray<double> GeneralTracksD0 = {fReader, "GeneralTracksD0"};
  TTreeReaderArray<double> GeneralTracksDZ = {fReader, "GeneralTracksDZ"};
  TTreeReaderArray<double> GeneralTracksEta = {fReader, "GeneralTracksEta"};
  TTreeReaderArray<double> GeneralTracksPhi = {fReader, "GeneralTracksPhi"};
  TTreeReaderArray<double> GeneralTracksPt = {fReader, "GeneralTracksPt"};
  TTreeReaderArray<double> PFParEta = {fReader, "PFParEta"};
  TTreeReaderArray<double> PFParM = {fReader, "PFParM"};
  TTreeReaderArray<double> PFParPhi = {fReader, "PFParPhi"};
  TTreeReaderArray<double> PFParPt = {fReader, "PFParPt"};
  
  //TTreeReaderArray<float> HBHERecHitEnergy = {fReader, "HBHERecHitEnergy"};
  //TTreeReaderArray<float> HBHERecHitEta = {fReader, "HBHERecHitEta"};
  //TTreeReaderArray<float> HBHERecHitPhi = {fReader, "HBHERecHitPhi"};
  //TTreeReaderArray<float> HBHERecHitTime = {fReader, "HBHERecHitTime"};
  TTreeReaderArray<float> HGCRecHitEnergy = {fReader, "HGCRecHitEnergy"};
  TTreeReaderArray<float> HGCRecHitEta = {fReader, "HGCRecHitEta"};
  TTreeReaderArray<float> HGCRecHitPhi = {fReader, "HGCRecHitPhi"};
  TTreeReaderArray<float> HGCRecHitPosx = {fReader, "HGCRecHitPosx"};
  TTreeReaderArray<float> HGCRecHitPosy = {fReader, "HGCRecHitPosy"};
  TTreeReaderArray<float> HGCRecHitPosz = {fReader, "HGCRecHitPosz"};
  TTreeReaderArray<float> HGCSimHitsEnergy = {fReader, "HGCSimHitsEnergy"};
  TTreeReaderArray<float> HGCSimHitsEta = {fReader, "HGCSimHitsEta"};
  TTreeReaderArray<float> HGCSimHitsPhi = {fReader, "HGCSimHitsPhi"};
  TTreeReaderArray<float> HGCSimHitsPosx = {fReader, "HGCSimHitsPosx"};
  TTreeReaderArray<float> HGCSimHitsPosy = {fReader, "HGCSimHitsPosy"};
  TTreeReaderArray<float> HGCSimHitsPosz = {fReader, "HGCSimHitsPosz"};
  //TTreeReaderArray<float> HGCSimHitsTime = {fReader, "HGCSimHitsTime"};
  //TTreeReaderArray<float> SimTracksEta = {fReader, "SimTracksEta"};
  //TTreeReaderArray<float> SimTracksPhi = {fReader, "SimTracksPhi"};
  //TTreeReaderArray<float> SimTracksPt = {fReader, "SimTracksPt"};
  //TTreeReaderArray<float> SimTracksR = {fReader, "SimTracksR"};
  //TTreeReaderArray<float> SimTracksZ = {fReader, "SimTracksZ"};
  
  TTreeReaderArray<float> PFParEcalEnergyFrac = {fReader, "PFParEcalEnergyFrac"};
  TTreeReaderArray<float> PFParHOEnergyFrac = {fReader, "PFParHOEnergyFrac"};
  TTreeReaderArray<float> PFParHcalEnergyFrac = {fReader, "PFParHcalEnergyFrac"};
  TTreeReaderArray<float> PFParHcalFrac1 = {fReader, "PFParHcalFrac1"};
  TTreeReaderArray<float> PFParHcalFrac2 = {fReader, "PFParHcalFrac2"};
  TTreeReaderArray<float> PFParHcalFrac3 = {fReader, "PFParHcalFrac3"};
  TTreeReaderArray<float> PFParHcalFrac4 = {fReader, "PFParHcalFrac4"};
  TTreeReaderArray<float> PFParHcalFrac5 = {fReader, "PFParHcalFrac5"};
  TTreeReaderArray<float> PFParHcalFrac6 = {fReader, "PFParHcalFrac6"};
  TTreeReaderArray<float> PFParHcalFrac7 = {fReader, "PFParHcalFrac7"};
  TTreeReaderArray<float> PFParTrackPt = {fReader, "PFParTrackPt"};
  TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
  TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
  TTreeReaderArray<int> GeneralTracksNValidHits = {fReader, "GeneralTracksNValidHits"};
  
  TTreeReaderArray<double> GenJetsPt = {fReader, "GenJetsPt"};
  TTreeReaderArray<double> GenJetsEta = {fReader, "GenJetsEta"};
  TTreeReaderArray<double> GenJetsPhi = {fReader, "GenJetsPhi"};
  TTreeReaderArray<double> GenJetsEnergy = {fReader, "GenJetsEnergy"};
  TTreeReaderArray<double> PFJetsPt = {fReader, "PFJetsPt"};
  TTreeReaderArray<double> PFJetsEta = {fReader, "PFJetsEta"};
  TTreeReaderArray<double> PFJetsPhi = {fReader, "PFJetsPhi"};
  TTreeReaderArray<double> PFJetsEnergy = {fReader, "PFJetsEnergy"};
  
  TTreeReaderArray<float> PFJetsChargedHadronMultiplicity = {fReader, "PFJetsChargedHadronMultiplicity"};
  TTreeReaderArray<float> PFJetsElectronMultiplicity = {fReader, "PFJetsElectronMultiplicity"};
  TTreeReaderArray<float> PFJetsMuonMultiplicity = {fReader, "PFJetsMuonMultiplicity"};
  TTreeReaderArray<float> PFJetsNeutralHadronMultiplicity = {fReader, "PFJetsNeutralHadronMultiplicity"};
  TTreeReaderArray<float> PFJetsPhotonEnergyFraction = {fReader, "PFJetsPhotonEnergyFraction"};
  TTreeReaderArray<float> PFJetsPhotonMultiplicity = {fReader, "PFJetsPhotonMultiplicity"};
  
  TTreeReaderArray<float> PFJetsElectronEnergyFraction = {fReader, "PFJetsElectronEnergyFraction"};
  
  TTreeReaderArray<float> PFJetsrecoJetsHFEMEnergyFraction = {fReader, "PFJetsrecoJetsHFEMEnergyFraction"};
  TTreeReaderArray<float> PFJetsrecoJetsHFHadronEnergyFraction = {fReader, "PFJetsrecoJetsHFHadronEnergyFraction"};
  TTreeReaderArray<float> PFJetsrecoJetschargedEmEnergyfraction = {fReader, "PFJetsrecoJetschargedEmEnergyFraction"};
  TTreeReaderArray<float> PFJetsrecoJetschargedHadronEnergyFraction = {fReader, "PFJetsrecoJetschargedHadronEnergyFraction"};
  TTreeReaderArray<float> PFJetsrecoJetsmuonEnergyFraction = {fReader, "PFJetsrecoJetsmuonEnergyFraction"};
  TTreeReaderArray<float> PFJetsrecoJetsneutralEmEnergyFraction = {fReader, "PFJetsrecoJetsneutralEmEnergyFraction"};
  TTreeReaderArray<float> PFJetsrecoJetsneutralEnergyFraction = {fReader, "PFJetsrecoJetsneutralEnergyFraction"};
  
  TTreeReaderArray<double> PFClusterECALPt = {fReader, "PFClusterECALPt"};
  TTreeReaderArray<double> PFClusterECALEnergy = {fReader, "PFClusterECALEnergy"};
  TTreeReaderArray<double> PFClusterECALEta = {fReader, "PFClusterECALEta"};
  TTreeReaderArray<double> PFClusterECALPhi = {fReader, "PFClusterECALPhi"};
  TTreeReaderArray<double> PFClusterECALPosX = {fReader, "PFClusterECALPosX"};
  TTreeReaderArray<double> PFClusterECALPosY = {fReader, "PFClusterECALPosY"};
  TTreeReaderArray<double> PFClusterECALPosZ = {fReader, "PFClusterECALPosZ"};
  TTreeReaderArray<double> PFClusterECALSize = {fReader, "PFClusterECALSize"};
  TTreeReaderArray<float> PFClusterECALTime = {fReader, "PFClusterECALTime"};
  TTreeReaderArray<float> PFClusterECALTimeE = {fReader, "PFClusterECALTimeE"};
  TTreeReaderArray<double> PFClusterECALDepth = {fReader, "PFClusterECALDepth"};
  
  TTreeReaderArray<double> PFClusterHCALPt = {fReader, "PFClusterHCALPt"};
  TTreeReaderArray<double> PFClusterHCALEnergy = {fReader, "PFClusterHCALEnergy"};
  TTreeReaderArray<double> PFClusterHCALEta = {fReader, "PFClusterHCALEta"};
  TTreeReaderArray<double> PFClusterHCALPhi = {fReader, "PFClusterHCALPhi"};
  TTreeReaderArray<double> PFClusterHCALPosX = {fReader, "PFClusterHCALPosX"};
  TTreeReaderArray<double> PFClusterHCALPosY = {fReader, "PFClusterHCALPosY"};
  TTreeReaderArray<double> PFClusterHCALPosZ = {fReader, "PFClusterHCALPosZ"};
  TTreeReaderArray<double> PFClusterHCALSize = {fReader, "PFClusterHCALSize"};
  TTreeReaderArray<float> PFClusterHCALTime = {fReader, "PFClusterHCALTime"};
  TTreeReaderArray<float> PFClusterHCALTimeE = {fReader, "PFClusterHCALTimeE"};
  TTreeReaderArray<double> PFClusterHCALDepth = {fReader, "PFClusterHCALDepth"};
  
  TTreeReaderArray<double> PFClusterHGCalPt = {fReader, "PFClusterHGCalPt"};
  TTreeReaderArray<double> PFClusterHGCalEnergy = {fReader, "PFClusterHGCalEnergy"};
  TTreeReaderArray<double> PFClusterHGCalEta = {fReader, "PFClusterHGCalEta"};
  TTreeReaderArray<double> PFClusterHGCalPhi = {fReader, "PFClusterHGCalPhi"};
  TTreeReaderArray<double> PFClusterHGCalPosX = {fReader, "PFClusterHGCalPosX"};
  TTreeReaderArray<double> PFClusterHGCalPosY = {fReader, "PFClusterHGCalPosY"};
  TTreeReaderArray<double> PFClusterHGCalPosZ = {fReader, "PFClusterHGCalPosZ"};
  TTreeReaderArray<double> PFClusterHGCalSize = {fReader, "PFClusterHGCalSize"};
  TTreeReaderArray<float> PFClusterHGCalTime = {fReader, "PFClusterHGCalTime"};
  TTreeReaderArray<float> PFClusterHGCalTimeE = {fReader, "PFClusterHGCalTimeE"};
  TTreeReaderArray<double> PFClusterHGCalDepth = {fReader, "PFClusterHGCalDepth"};
  
  TTreeReaderArray<double> PFClusterHFPt = {fReader, "PFClusterHFPt"};
  TTreeReaderArray<double> PFClusterHFEnergy = {fReader, "PFClusterHFEnergy"};
  TTreeReaderArray<double> PFClusterHFEta = {fReader, "PFClusterHFEta"};
  TTreeReaderArray<double> PFClusterHFPhi = {fReader, "PFClusterHFPhi"};
  TTreeReaderArray<double> PFClusterHFPosX = {fReader, "PFClusterHFPosX"};
  TTreeReaderArray<double> PFClusterHFPosY = {fReader, "PFClusterHFPosY"};
  TTreeReaderArray<double> PFClusterHFPosZ = {fReader, "PFClusterHFPosZ"};
  TTreeReaderArray<double> PFClusterHFSize = {fReader, "PFClusterHFSize"};
  TTreeReaderArray<float> PFClusterHFTime = {fReader, "PFClusterHFTime"};
  TTreeReaderArray<float> PFClusterHFTimeE = {fReader, "PFClusterHFTimeE"};
  TTreeReaderArray<double> PFClusterHFDepth = {fReader, "PFClusterHFDepth"};
  
  TTreeReaderArray<double> PFClusterHOPt = {fReader, "PFClusterHOPt"};
  TTreeReaderArray<double> PFClusterHOEnergy = {fReader, "PFClusterHOEnergy"};
  TTreeReaderArray<double> PFClusterHOEta = {fReader, "PFClusterHOEta"};
  TTreeReaderArray<double> PFClusterHOPhi = {fReader, "PFClusterHOPhi"};
  TTreeReaderArray<double> PFClusterHOPosX = {fReader, "PFClusterHOPosX"};
  TTreeReaderArray<double> PFClusterHOPosY = {fReader, "PFClusterHOPosY"};
  TTreeReaderArray<double> PFClusterHOPosZ = {fReader, "PFClusterHOPosZ"};
  TTreeReaderArray<double> PFClusterHOSize = {fReader, "PFClusterHOSize"};
  TTreeReaderArray<float> PFClusterHOTime = {fReader, "PFClusterHOTime"};
  TTreeReaderArray<float> PFClusterHOTimeE = {fReader, "PFClusterHOTimeE"};
  TTreeReaderArray<double> PFClusterHODepth = {fReader, "PFClusterHODepth"};

  TTreeReaderArray<double> PFClusterHGCalFromMultiCLPt = {fReader, "PFClusterHGCalFromMultiCLPt"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLEnergy = {fReader, "PFClusterHGCalFromMultiCLEnergy"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLEta = {fReader, "PFClusterHGCalFromMultiCLEta"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLPhi = {fReader, "PFClusterHGCalFromMultiCLPhi"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLPosX = {fReader, "PFClusterHGCalFromMultiCLPosX"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLPosY = {fReader, "PFClusterHGCalFromMultiCLPosY"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLPosZ = {fReader, "PFClusterHGCalFromMultiCLPosZ"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLSize = {fReader, "PFClusterHGCalFromMultiCLSize"};
  TTreeReaderArray<float> PFClusterHGCalFromMultiCLTime = {fReader, "PFClusterHGCalFromMultiCLTime"};
  TTreeReaderArray<float> PFClusterHGCalFromMultiCLTimeE = {fReader, "PFClusterHGCalFromMultiCLTimeE"};
  TTreeReaderArray<double> PFClusterHGCalFromMultiCLDepth = {fReader, "PFClusterHGCalFromMultiCLDepth"};
  
  //TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
  //TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
  //TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
  //TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
  //TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
  //TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
  //TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
  TTreeReaderArray<int> HGCRecHitIndex = {fReader, "HGCRecHitIndex"};
  TTreeReaderArray<int> HGCRecHitLayer = {fReader, "HGCRecHitLayer"};
  TTreeReaderArray<int> HGCSimHitsCellU = {fReader, "HGCSimHitsCellU"};
  TTreeReaderArray<int> HGCSimHitsCellV = {fReader, "HGCSimHitsCellV"};
  TTreeReaderArray<int> HGCSimHitsIEta = {fReader, "HGCSimHitsIEta"};
  TTreeReaderArray<int> HGCSimHitsIPhi = {fReader, "HGCSimHitsIPhi"};
  TTreeReaderArray<int> HGCSimHitsIndex = {fReader, "HGCSimHitsIndex"};
  TTreeReaderArray<int> HGCSimHitsLayer = {fReader, "HGCSimHitsLayer"};
  TTreeReaderArray<int> HGCSimHitsSubdet = {fReader, "HGCSimHitsSubdet"};
  TTreeReaderArray<int> HGCSimHitsWaferU = {fReader, "HGCSimHitsWaferU"};
  TTreeReaderArray<int> HGCSimHitsWaferV = {fReader, "HGCSimHitsWaferV"};
  //TTreeReaderArray<int> SimTracksCharge = {fReader, "SimTracksCharge"};
  //TTreeReaderArray<int> SimTracksPID = {fReader, "SimTracksPID"};
  
  TTreeReaderArray<int> PFParPdgId = {fReader, "PFParPdgId"};
  TTreeReaderArray<int> PFParStatus = {fReader, "PFParStatus"};
  TTreeReaderArray<double> PFParMET = {fReader, "PFParMET"};
  TTreeReaderArray<double> PFMET = {fReader, "PFMET"};
  TTreeReaderArray<double> GenMET = {fReader, "GenMET"};
  TTreeReaderValue<UInt_t> bx = {fReader, "bx"};
  TTreeReaderValue<UInt_t> event = {fReader, "event"};
  TTreeReaderValue<UInt_t> ls = {fReader, "ls"};
  TTreeReaderValue<UInt_t> orbit = {fReader, "orbit"};
  TTreeReaderValue<UInt_t> run = {fReader, "run"};

}

#endif
