// ------------------------------------------------------------------------------------
//  ROOT macro that produces average RecHit energy from PFG ntuples
//
//  Author : Ken H
//  Written on May 24, 2018
// ------------------------------------------------------------------------------------
//  
// Pre-requisite :
//
//   You should have the PFG ntuple for the Run from which you want to do a measurement. 
//   Instruction on how to make PFG ntuples can be found here : FIXME link here 
//
//   You should have "Fig" directory for plots 
//
// Usage : 
//
//   $ root -b  
//   root> .L ana_PFStudy.C+
//   root> ana_PFStudy("/cms/data/store/user/hatake/HCAL/ntuples/10_2_x/pi50_trees_MCfull_CMSSW_10_2_0_pre3_*.root","hcal_timestudy_pi50_histograms.root")
//   or
//   root> ana_PFStudy("list_trees_pi50_MCfull_CMSSW_10_2_0_pre3.txt","hcal_timestudy_pi50_histograms.root")
//   or
//   from command line:
/*
     root.exe -b -q 'ana_PFStudy.C++("trees_relval_ttbar_phase2_age_new2_4500ultimate.root","hcal_noisestudy_histograms_age_new2_4500ultimate.root")'
     root.exe -b -q 'ana_PFStudy.C++("trees_relval_ttbar_phase2_age_org.root","hcal_noisestudy_histograms_age_org.root")'
     root.exe -b -q 'ana_PFStudy.C++("trees_relval_ttbar_phase2_noage.root","hcal_noisestudy_histograms_noage.root")'
 */
//    
// -----------------------------------------------------------------------------------
// 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()
#include <algorithm> 

#include "../include/fReader.h" // Read from TTRee

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

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;
using namespace globalTChain;

bool DRAWPLOTS  = false;  // draw plots or not (make "Fig" directory first before turning this on)
bool VERBOSE    = false;  // print out mean +/- sigma for each channel or not

// Assemble a list of inputfiles

std::vector<std::string> GetInputFiles(TString geoConfig)
{
  int numFiles = 20;
  std::string path = "/cms/data/store/user/bcaraway/condor/outputs/";
  std::string startName = "ttbar_10_4_";
  std::string midName = "_pt25_numEvent1000_CMSSW_10_4_0_pre2_";
  std::string endName = "_v01.root";
  std::vector<std::string> inputFiles;
  
  for( int iFile = 0 ; iFile<numFiles ; iFile++ )
    {
      std::ostringstream fileName;
      fileName << path << startName << geoConfig << midName << iFile << endName;
      inputFiles.push_back(fileName.str());
     
    }

  return inputFiles;

}

//
// Book histograms
//

void book1D(TList *v_hist, std::string name, int n, double min, double max);
void book1DProf(TList *v_hist, std::string name, int n, double min, double max, double ymin, double ymax, Option_t *option);

void book2D(TList *v_hist, std::string name, int xn, double xmin, double xmax, int yn, double ymin, double ymax); // added by Bryan

void bookHistograms(TList *v_hist);

//
// Fill histograms
//
void fill1D(TList *v_hist, std::string name, double value);
void fill1D(TList *v_hist, std::string name, double value, double valuey);
void fill1DProf(TList *v_hist, std::string name, double value, double valuey);

void fill2D(TList *v_hist, std::string name, double valuex, double valuey);  // added by Bryan 

//
// Aux
//
void relabelProfA(TList *v_hist, std::string name);
void relabel2D   (TList *v_hist, std::string name); // added by Bryan
//
// Main analyzer
//
void PFCheckRun(std::vector<std::string> inputFiles, TString outfile, int maxevents, int option=2) 
{ 

   cout << "[PF analyzer] Running option " << option << " for " << endl; 

   // fit pannel display option
   gStyle->SetOptFit(1011);

   //
   // Get the tree from the PFG ntuple 
   //
   TChain *ch = new TChain("hgcalTupleTree/tree");

   //std::string filename(rootfile);
   //std::string::size_type idx;
   //idx = filename.rfind('.');
   //std::string extension = filename.substr(idx+1);
   //std::string line;
   
   //if(idx != std::string::npos && extension=="txt")
   //  {
   //    std::cout << rootfile << " " << extension << std::endl;
   //    std::ifstream in(rootfile);
   //    while (std::getline(in, line)) {     // Process line
   //	 if (line.size()>0) ch->Add(line.c_str());
   //    }
   //  }
   //else
   //  {
   //    // No extension found
   //    ch->Add(rootfile);
   //  }


   for (unsigned int iFile=0; iFile<inputFiles.size(); ++iFile) {
    ch->Add(inputFiles[iFile].c_str());
    std::cout<<inputFiles[iFile]<<std::endl;
   }

   printf("%d;\n",ch->GetNtrees());
   printf("%lld;\n",ch->GetEntries());
   

   fReader.SetTree(ch);  //the tree reader (Defined in fReader.h)

   //
   // Define histograms to fill
   //
   TList *v_hist = new TList();
   
   bookHistograms(v_hist); // most of histograms booked here
   
   //
   // Loop over entries
   //
   unsigned int nentries = (Int_t)ch->GetEntries();
   cout << "[HGCal analyzer] The number of entries is: " << nentries << endl;

   bool debug=false;
   bool debug_met =false;
   bool print_prev = false;

   //---------------------------------------------------------------------------------------------------------
   // main event loop
   //---------------------------------------------------------------------------------------------------------
   
   TLorentzVector tlzv;      // for MET
   TLorentzVector tlzv_temp; // for MET
   
   std::map<int, TLorentzVector> pfpar_tlzv;
   std::map<int, TString>        pfpar_type;
   pfpar_type[1] = "chargedHadron"; pfpar_type[2] = "electron"; pfpar_type[3] = "muon";
   pfpar_type[4] = "photon"; pfpar_type[5] = "neutralHadron"; pfpar_type[6] = "HFHadron";
   pfpar_type[7] = "HFPhoton";
   								  
   
   int ievent=0;
   while (fReader.Next()) {
  
     // Progress indicator 
     ievent++;
     if(ievent%1000==0) cout << "[HCAL analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
     if (maxevents>0 && ievent>maxevents) Break;
     
     //--------------------
     // Loop over PF candidates
     //--------------------
     
     tlzv.SetPtEtaPhiM(0.0,0.0,0.0,0.0);       // for MET
     tlzv_temp.SetPtEtaPhiM(0.0,0.0,0.0,0.0);  // for MET
     for (int i = 1; i <= 7 ; i++){            // for summing pt of respective pfparticle
       pfpar_tlzv[i].SetPtEtaPhiM(0.0,0.0,0.0,0.0); 
     }
     double delta_phi = 2*TMath::Pi();
     double delta_eta = 0.1;

     double PFMass = 0; // for hardcoded mass in GeV, PFM suffers from a rounding error.     
     

     // Begin Loop
     for (int ipfcand = 0, npfcand =  PFParPt.GetSize(); ipfcand < npfcand; ++ipfcand) {
       
       //if (PFParPdgId[ipfcand] == 7) std::cout << "PdgID: " << PFParPdgId[ipfcand] <<", with Mass: "<<PFParM[ipfcand]<<std::endl; // to debug PFParM
       
       // set the mass of PF cand

       
       if (PFParPdgId[ipfcand] == 1) PFMass = .13957 ; // charged hadron --> Pion
       if (PFParPdgId[ipfcand] == 2) PFMass = .000511 ; // electron
       if (PFParPdgId[ipfcand] == 3) PFMass = .1057 ; // muon
       if (PFParPdgId[ipfcand] == 4) PFMass = 0. ; // photon
       if (PFParPdgId[ipfcand] == 5) PFMass = .49761 ; // neutral hadron --> K0L
       if (PFParPdgId[ipfcand] == 6) PFMass = .13957 ; // HF Hadron --> Pion
       if (PFParPdgId[ipfcand] == 7) PFMass = 0. ; // HF EM Particle --> Photon
	 
       tlzv_temp.SetPtEtaPhiM(PFParPt[ipfcand],PFParEta[ipfcand],PFParPhi[ipfcand],PFMass);
       std::string strtmp;

       pfpar_tlzv[PFParPdgId[ipfcand]] += tlzv_temp;
       ///=================================PT vs Eta======================
       
       fill1D(v_hist, static_cast<std::string>("PFPt_Eta_"+pfpar_type[PFParPdgId[ipfcand]]), tlzv_temp.Eta(), tlzv_temp.Pt()/(nentries*delta_phi*delta_eta));
       
       ////
       //// Bryan's studies
       ////
       strtmp = "PFTask_PdgId_vs_pt";  
       fill2D(v_hist, strtmp, PFParPdgId[ipfcand] , PFParPt[ipfcand]);       
       // Plotting PF eta
       strtmp = "PFTask_PFEta_All";
       fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       if (PFParPdgId[ipfcand] == 1){  // for Charged Hadron
	 strtmp = "PFTask_PFEta_ChargedHadron";
	 fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       }
       if (PFParPdgId[ipfcand] == 5){  // for Neutral Hadron
	 strtmp = "PFTask_PFEta_NeutralHadron";
	 fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       }
       if (PFParPdgId[ipfcand] == 2){  // for Electron
	 strtmp = "PFTask_PFEta_Electron";
	 fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       }
       if (PFParPdgId[ipfcand] == 3){  // for Muon
	 strtmp = "PFTask_PFEta_Muon";
	 fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       }
       if (PFParPdgId[ipfcand] == 4){  // for Photon
	 strtmp = "PFTask_PFEta_Photon";
	 fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       }
       if (PFParPdgId[ipfcand] == 6) { // for HF Hadron
	 strtmp = "PFTask_PFEta_HF_Hadron";
	 fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       }
       if (PFParPdgId[ipfcand] == 7) { // for HF Photon
	 strtmp = "PFTask_PFEta_HF_Photon";
	 fill1D(v_hist, strtmp, PFParEta[ipfcand]);
       }
       ///
       /// Adding up for MET
       ///
       
       
       if (debug_met){
	 std::cout << "Pt: " << PFParPt[ipfcand] << ", Eta: " << PFParEta[ipfcand] << ", Phi: " << PFParPhi[ipfcand] << ", Mass: " << PFMass << std::endl;	 
       }
       tlzv_temp.SetPtEtaPhiM(PFParPt[ipfcand],PFParEta[ipfcand],PFParPhi[ipfcand],PFMass);
       tlzv += tlzv_temp;
       

     } // PF candidiate loop
     for (int i = 1; i <= 7 ; i++){
       fill1D(v_hist, static_cast<std::string>("PFPtvsEta_"+pfpar_type[i]),  pfpar_tlzv[i].Eta(),  pfpar_tlzv[i].Pt()/(nentries*delta_phi*delta_eta));
     }
     // catagorize the event by gen particle makeup
     int n_neutral_gen = 0;
     int ngen = GenParPdgId.GetSize();
     for (int igen = 0; igen<ngen; igen++){
       if (GenParPdgId[igen] == 130 || GenParPdgId[igen] == 111 || GenParPdgId[igen] == 310) n_neutral_gen++; // counts neutral particles 
     }
     double gen_ratio = static_cast<double>(n_neutral_gen)/static_cast<double>(ngen);// std::cout<<gen_ratio<<"\t"<<n_neutral_gen<<"\t"<<ngen<<std::endl;
     fill1D(v_hist, "Neutral_Gen_Particles",n_neutral_gen);
     fill1D(v_hist, "Neutral_Gen_Ratio",gen_ratio);

     bool debug_met_V2    = false; 
     //=========ONLY SET ONE TO "TRUE"================
     bool low_GenMET        = false; // set to true if doing low Gen study
     bool jet_at_boundary   = true; // set to true if doing at least one jet at HGCal boundary study
     bool norm_jet          = true; // no constaints 
     bool large_neutral     = false; // jet makeup consists of a large amount of neutral hadrons/Em
     bool gen_neutral_ratio = false;  // gen particle ratio consists of a noticable amount of neutral particls
     //===============================================

     //if  (tlzv.Pt() > 200 && debug_met_V2){
     //  std::cout << "MET: "<<tlzv.Pt()<<std::endl;
     //}
     fill1D(v_hist, "PFTask_MET", tlzv.Pt());  // for MET
     
     fill1D(v_hist, "PFParMET", PFParMET[0]);     // MET calculated before ntuple step
     //fill1D(v_hist, "PFMET"   , PFMET[0]   );     //  MET
     fill1D(v_hist, "GenMET"  , GenMET[0]  );     // Gen MET
     
     if (GenMET[0] < 10 && low_GenMET) {
       fill1D(v_hist, "PFMET", PFMET[0]   ); // exclude events where met is largely due to neutrinos 
     }

     if (tlzv.Pt() - PFMET[0] > abs(5.)){ 
       std::cout << "\nCalculated MET:\t\t " << tlzv.Pt() <<"\nCalculated before ntuple MET:\t " <<PFParMET[0] << "\nTrue MET:\t\t "<<PFMET[0]<<std::endl;
     }
     ///
     /// DEBUG MET
     ///
     
     if ((GenMET[0]< 10 && PFMET[0] > 200) && debug_met_V2){
       
       tlzv.SetPtEtaPhiM(0.0,0.0,0.0,0.0);       // for MET
       tlzv_temp.SetPtEtaPhiM(0.0,0.0,0.0,0.0);  // for MET

       for (int ipfcand = 0, npfcand =  PFParPt.GetSize(); ipfcand < npfcand; ++ipfcand) {
	 
	 if (PFParPdgId[ipfcand] == 1) PFMass = .13957 ; // charged hadron --> Pion
	 if (PFParPdgId[ipfcand] == 2) PFMass = .000511 ; // electron
	 if (PFParPdgId[ipfcand] == 3) PFMass = .1057 ; // muon
	 if (PFParPdgId[ipfcand] == 4) PFMass = 0. ; // photon
	 if (PFParPdgId[ipfcand] == 5) PFMass = .49761 ; // neutral hadron --> K0L
	 if (PFParPdgId[ipfcand] == 6) PFMass = .13957 ; // HF Hadron --> Pion
	 if (PFParPdgId[ipfcand] == 7) PFMass = 0. ; // HF EM Particle --> Photon

	 tlzv_temp.SetPtEtaPhiM(PFParPt[ipfcand],PFParEta[ipfcand],PFParPhi[ipfcand],PFMass);
	 tlzv += tlzv_temp;
	 
	 if (PFParPt[ipfcand] > 10){
	   std::cout << "PFID: "<< PFParPdgId[ipfcand] << ", Pt: " << PFParPt[ipfcand] << ", Eta: " << PFParEta[ipfcand] << ", Phi: " << PFParPhi[ipfcand] << ", Mass: " << PFMass << std::endl;
	   //std::cout<< "Summed PT: " << tlzv.Pt() << std::endl;
	 }
       }
     }

     ///
     /// END DEBUG MET
     ///

     ///
     /// Gen Jet and PF Jet analysis
     ///
     bool boundary_eta = false;
     int gen_boundary_jets = 0;
     TLorentzVector genJet_tlzv; genJet_tlzv.SetPtEtaPhiE(0.,0.,0.,0.);
     std::vector<TLorentzVector> tlzv_vec;
     TLorentzVector genJet_tlzv_hcal; genJet_tlzv.SetPtEtaPhiE(0.,0.,0.,0.);
     std::vector<TLorentzVector> tlzv_vec_hcal;
     TLorentzVector genJet_tlzv_base; genJet_tlzv.SetPtEtaPhiE(0.,0.,0.,0.);
     std::vector<TLorentzVector> tlzv_vec_base;
     // Gen
     tlzv.SetPtEtaPhiE(0.0,0.0,0.0,0.0);       // for MET                                                                                                 
     tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);  // for MET   

     for (int ijet = 0, njet = GenJetsPt.GetSize(); ijet < njet; ++ijet) {      
       tlzv_temp.SetPtEtaPhiE(GenJetsPt[ijet],GenJetsEta[ijet],GenJetsPhi[ijet],GenJetsEnergy[ijet]);
       tlzv += tlzv_temp;
       fill1D(v_hist,"GenJetEta",GenJetsEta[ijet]);
       double absJetEta = abs(GenJetsEta[ijet]);
       if ((absJetEta <= 3.10 && absJetEta > 2.80) && GenJetsPt[ijet] >= 20){// || (absJetEta < 1.6 && absJetEta >= 1.4)){
	 gen_boundary_jets++;
	 boundary_eta = true;
	 genJet_tlzv = tlzv_temp;
	 tlzv_vec.push_back(genJet_tlzv);
       }
       if ((absJetEta <= 1.60 && absJetEta > 1.40) && GenJetsPt[ijet] >= 20){// || (absJetEta < 1.6 && absJetEta >= 1.4)){
	 boundary_eta = true;
	 genJet_tlzv_hcal = tlzv_temp;
	 tlzv_vec_hcal.push_back(genJet_tlzv_hcal);
       }
       if ((absJetEta <= 2.80 && absJetEta > 2.40) && GenJetsPt[ijet] >= 20){// || (absJetEta < 1.6 && absJetEta >= 1.4)){
	 boundary_eta = true;
	 genJet_tlzv_base = tlzv_temp;
	 tlzv_vec_base.push_back(genJet_tlzv_base);
       }
     }
     //if (boundary_eta == true) std::cout<<"Gen jets at boundary:\t\t"<<gen_boundary_jets<<std::endl;
     fill1D(v_hist,"GenJetMET",tlzv.Pt());

     // PF

    
     bool neutral_event = false; 
     int n_jets_eta = 0;
     int n_neutral_jets = 0;

     tlzv.SetPtEtaPhiE(0.0,0.0,0.0,0.0);       // for MET                                                                                                 
     tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);  // for MET   
     
     for (int ijet = 0, njet = PFJetsPt.GetSize(); ijet < njet; ++ijet) {
      
       tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet],PFJetsEta[ijet],PFJetsPhi[ijet],PFJetsEnergy[ijet]);
       tlzv += tlzv_temp;
       
       if (norm_jet) fill1D(v_hist,"PFJetEta",PFJetsEta[ijet]);
       
       fill1D(v_hist,"PFJetHFEMEFraction",PFJetsrecoJetsHFEMEnergyFraction[ijet]);
       fill1D(v_hist,"PFJetchargedHFHadronEFraction" ,PFJetsrecoJetsHFHadronEnergyFraction[ijet]);
       fill1D(v_hist,"PFJetchargedEMEFraction" ,PFJetsrecoJetschargedEmEnergyfraction[ijet]);
       fill1D(v_hist,"PFJetchargedHadronEFraction" ,PFJetsrecoJetschargedHadronEnergyFraction[ijet]);
       fill1D(v_hist,"PFJetMuonEFraction" ,PFJetsrecoJetsmuonEnergyFraction[ijet]);
       fill1D(v_hist,"PFJetneutralEMEFraction" ,PFJetsrecoJetsneutralEmEnergyFraction[ijet]);
       fill1D(v_hist,"PFJetneutralEFraction" ,PFJetsrecoJetsneutralEnergyFraction[ijet]);
       
       double absJetEta = abs(PFJetsEta[ijet]);

       if ((absJetEta <= 3.0 && absJetEta > 2.8) || (absJetEta < 1.6 && absJetEta >= 1.4)){
	 //boundary_eta = true;
	 n_jets_eta++ ;
       }
     
       
       if ((PFJetsrecoJetsneutralEnergyFraction[ijet] + PFJetsrecoJetsneutralEmEnergyFraction[ijet]) > .30) n_neutral_jets++;
       //if (gen_ratio > .09) fill1D(v_hist,"PFtotalneutralEFraction", PFJetsrecoJetsneutralEnergyFraction[ijet] + PFJetsrecoJetsneutralEmEnergyFraction[ijet]);
     }	 
     
     fill1D(v_hist, "n_neutral_jets", n_neutral_jets);
     if (n_neutral_jets > 9) neutral_event = true; // event has a large amount of neutral jets

     //std::cout<<"# of Jets at the boundary:\t"<<n_jets_eta<<std::endl;
     if (norm_jet) fill1D(v_hist, "PFJetMET", tlzv.Pt());
     //  jet at boundary ===================================================================================
     if (boundary_eta && jet_at_boundary){ // at least one jet at one of the hgcal boundaries
       // for the beampipe, eta = 3.0 resolution
       for (int igen = 0; igen != tlzv_vec.size(); igen++){
	 tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);         
	 int ijet_store = 0;
	 double min_dR = .21;
	 bool matched_jet = false;
	 for (int ijet = 0, njet = PFJetsPt.GetSize(); ijet < njet; ++ijet) {
	   tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet],PFJetsEta[ijet],PFJetsPhi[ijet],PFJetsEnergy[ijet]);
	   if (tlzv_vec[igen].DeltaR(tlzv_temp) <= 0.2){
	     matched_jet = true;
	     if ( tlzv_vec[igen].DeltaR(tlzv_temp) < min_dR ){
	       min_dR = tlzv_vec[igen].DeltaR(tlzv_temp);
	       ijet_store = ijet;
	     }
	   }
	 }
	 if (matched_jet) {
	   tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet_store],PFJetsEta[ijet_store],PFJetsPhi[ijet_store],PFJetsEnergy[ijet_store]);
	   double jet_reso = (tlzv_temp.Pt() - tlzv_vec[igen].Pt())/tlzv_vec[igen].Pt();
	   fill1D(v_hist, "PFJet_Reso_beam", jet_reso);
	 }
       }
       // for outer boundary, eta = 1.4 resolution
       for (int igen = 0; igen != tlzv_vec_hcal.size(); igen++){
	 tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);         
	 int ijet_store = 0;
	 double min_dR = .21;
	 bool matched_jet = false;
	 for (int ijet = 0, njet = PFJetsPt.GetSize(); ijet < njet; ++ijet) {
	   tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet],PFJetsEta[ijet],PFJetsPhi[ijet],PFJetsEnergy[ijet]);
	   if (tlzv_vec_hcal[igen].DeltaR(tlzv_temp) <= 0.2){
	     matched_jet = true;
	     if ( tlzv_vec_hcal[igen].DeltaR(tlzv_temp) < min_dR ){
	       min_dR = tlzv_vec_hcal[igen].DeltaR(tlzv_temp);
	       ijet_store = ijet;
	     }
	   }
	 }
	 if (matched_jet) {
	   tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet_store],PFJetsEta[ijet_store],PFJetsPhi[ijet_store],PFJetsEnergy[ijet_store]);
	   double jet_reso_hcal = (tlzv_temp.Pt() - tlzv_vec_hcal[igen].Pt())/tlzv_vec_hcal[igen].Pt();
	   fill1D(v_hist, "PFJet_Reso_hcal", jet_reso_hcal);
	 }
       }
       // for the baseline resolution
       for (int igen = 0; igen != tlzv_vec_base.size(); igen++){
	 tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);         
	 int ijet_store = 0;
	 double min_dR = .21;
	 bool matched_jet = false;
	 for (int ijet = 0, njet = PFJetsPt.GetSize(); ijet < njet; ++ijet) {
	   tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet],PFJetsEta[ijet],PFJetsPhi[ijet],PFJetsEnergy[ijet]);
	   if (tlzv_vec_base[igen].DeltaR(tlzv_temp) <= 0.2){
	     matched_jet = true;
	     if ( tlzv_vec_base[igen].DeltaR(tlzv_temp) < min_dR ){
	       min_dR = tlzv_vec_base[igen].DeltaR(tlzv_temp);
	       ijet_store = ijet;
	     }
	   }
	 }
	 if (matched_jet) {
	   tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet_store],PFJetsEta[ijet_store],PFJetsPhi[ijet_store],PFJetsEnergy[ijet_store]);
	   double jet_reso_base = (tlzv_temp.Pt() - tlzv_vec_base[igen].Pt())/tlzv_vec_base[igen].Pt();
	   fill1D(v_hist, "PFJet_Reso_base", jet_reso_base);
	 }
       }
       //============
     }// endif ======
     // low GenMET =========================================================================================
     if (GenMET[0] < 10 && low_GenMET) { // requiring gen met to be low elliminates high neutrino events 
       tlzv.SetPtEtaPhiE(0.0,0.0,0.0,0.0);       // for MET                                                                                                 
       tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);  // for MET   
     
       for (int ijet = 0, njet = PFJetsPt.GetSize(); ijet < njet; ++ijet) {
       
	 tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet],PFJetsEta[ijet],PFJetsPhi[ijet],PFJetsEnergy[ijet]);
	 tlzv += tlzv_temp;
       
	 fill1D(v_hist,"PFJetEta",PFJetsEta[ijet]);
       }
       fill1D(v_hist, "PFJetMET", tlzv.Pt());
     }
     // large amount of neutral jet events ================================================================
     if (neutral_event && large_neutral){
       tlzv.SetPtEtaPhiE(0.0,0.0,0.0,0.0);       // for MET                                                                                                 
       tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);  // for MET   
     
       for (int ijet = 0, njet = PFJetsPt.GetSize(); ijet < njet; ++ijet) {
       
	 tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet],PFJetsEta[ijet],PFJetsPhi[ijet],PFJetsEnergy[ijet]);
	 tlzv += tlzv_temp;
       
	 fill1D(v_hist,"PFJetEta",PFJetsEta[ijet]);
       }
       fill1D(v_hist, "PFJetMET", tlzv.Pt());
       fill1D(v_hist, "PFMET", PFMET[0]   );
     } 
     // decent amount of neutral gen particles ===========================================================
     if (gen_ratio > .09 && gen_neutral_ratio){
       tlzv.SetPtEtaPhiE(0.0,0.0,0.0,0.0);       // for MET                                                                                                 
       tlzv_temp.SetPtEtaPhiE(0.0,0.0,0.0,0.0);  // for MET   
       
       for (int ijet = 0, njet = PFJetsPt.GetSize(); ijet < njet; ++ijet) {
       
	 tlzv_temp.SetPtEtaPhiE(PFJetsPt[ijet],PFJetsEta[ijet],PFJetsPhi[ijet],PFJetsEnergy[ijet]);
	 tlzv += tlzv_temp;
       
	 fill1D(v_hist,"PFJetEta",PFJetsEta[ijet]);
       }
       fill1D(v_hist, "PFJetMET", tlzv.Pt());
       fill1D(v_hist, "PFMET", PFMET[0]   );
     } 
      
   }   // Event loop ends
   //---------------------------------------------------------------------------------------------------------
   // main event loop ends
   //---------------------------------------------------------------------------------------------------------

   // output file for histograms
   std::cout<<"Output file: \t"<<outfile<<std::endl;
   TFile file_out(outfile,"RECREATE");

   v_hist->Write();
   
   file_out.ls();
   file_out.Close();

}

//
// Main function
//
//void ana_PFStudy(TString rootfile="../../HGCalTreeMaker/test/ttbar_10_4_D30_pt25.root",TString outfile="pfstudy_histograms.root",int maxevents=-1)
// "D30" for D30 geo, "D28" for D28
void ana_main(TString geoType)
{
  int maxevents=20000;
  // edit 
  bool test_file = false; // if testing setup with single file (will have to edit below for file choice)
  TString rootfile = "./ttbar_10_4_"+geoType+"_pt25.root";
  TString outfile  = "PF"+geoType+"_histos.root";

  std::vector<std::string> inputFiles;
  if (!test_file){
      inputFiles = GetInputFiles(geoType);
    }
  else inputFiles.push_back(static_cast<std::string>(rootfile));

  PFCheckRun(inputFiles, outfile, maxevents, 0);
}

void ana_PFStudy(std::vector<TString> geoTypes = {"D28","D30"}){ // loop over HGCal geometry configurations
  for (int i = 0; i != geoTypes.size(); i++){
    ana_main(geoTypes[i]);
  }
}
 

//
// --- Aux ---
//

//
// Book 1D histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max)
{
  TH1D *htemp = new TH1D(name.c_str(), name.c_str(), n, min, max);
  v_hist->Add(htemp);
}
//
// Book 1D profile histograms
//
void book1DProf(TList *v_hist, std::string name, int n, double min, double max, double ymin, double ymax, Option_t *option="")
{
  TProfile *htemp = new TProfile(name.c_str(), name.c_str(), n, min, max, ymin, ymax, option);
  v_hist->Add(htemp);
}
//
// Book 2D profile histograms 
//
void book2D(TList *v_hist, std::string name, int xn, double xmin, double xmax, int yn, double ymin, double ymax) 
{
  TH2D *htemp = new TH2D(name.c_str(), name.c_str(), xn, xmin, xmax, yn, ymin, ymax);
  v_hist->Add(htemp);
}
//
// Book histograms
//
void bookHistograms(TList *v_hist)
{

  Char_t histo[100];
  std::string strtmp;
  
  //
  // Booking histograms
  // 

  sprintf(histo, "PFTask_PdgId_vs_pt");          
  book2D(    v_hist, histo, 10., 0.5, 10.5, 100., 0., 10.); 
  relabel2D( v_hist, histo);

  std::map<int, TString>        pfpar_type;
  pfpar_type[1] = "chargedHadron"; pfpar_type[2] = "electron"; pfpar_type[3] = "muon";
  pfpar_type[4] = "photon"; pfpar_type[5] = "neutralHadron"; pfpar_type[6] = "HFHadron";
  pfpar_type[7] = "HFPhoton";  
  TString pfpar_pre = "PFPt_Eta_";

  for (int i = 1; i <= 7; i++){
    book1D( v_hist, static_cast<std::string>("PFPt_Eta_" +pfpar_type[i]), 100, -5, 5);
  }
  for (int i = 1; i <= 7; i++){
    book1D( v_hist, static_cast<std::string>("PFPtvsEta_"+pfpar_type[i]), 100, -5, 5);
  }
  ////
  //// Bryan's Studies 
  ////

  sprintf(histo,"Neutral_Gen_Particles");
  book1D( v_hist, histo, 100, 0. , 500);
  sprintf(histo,"Neutral_Gen_Ratio");
  book1D( v_hist, histo, 110, 0. , 1.1);
  
  sprintf(histo, "PFTask_PFEta_All");
  book1D(v_hist, histo, 100, -5., 5.);
  
  sprintf(histo, "PFTask_PFEta_ChargedHadron");
  book1D(v_hist, histo, 100, -5., 5.);
  sprintf(histo, "PFTask_PFEta_NeutralHadron");
  book1D(v_hist, histo, 100, -5., 5.);
  sprintf(histo, "PFTask_PFEta_Electron");
  book1D(v_hist, histo, 100, -5., 5.);
  sprintf(histo, "PFTask_PFEta_Muon");
  book1D(v_hist, histo, 100, -5., 5.);
  sprintf(histo, "PFTask_PFEta_Photon");
  book1D(v_hist, histo, 100, -5., 5.);
  sprintf(histo, "PFTask_PFEta_HF_Hadron");
  book1D(v_hist, histo, 100, -5., 5.);
  sprintf(histo, "PFTask_PFEta_HF_Photon");
  book1D(v_hist, histo, 100, -5., 5.);

  ///
  /// MET
  ///
  
  sprintf(histo, "PFTask_MET");
  book1D(v_hist, histo, 300, 0., 300.);
  sprintf(histo, "PFParMET");
  book1D(v_hist, histo, 300, 0., 300.);
  sprintf(histo, "PFMET");
  book1D(v_hist, histo, 300, 0., 300.);

  sprintf(histo, "GenMET");
  book1D(v_hist, histo, 300, 0., 300.);

  ///
  /// Jets
  ///

  sprintf(histo, "PFJet_Reso_beam");
  book1D(v_hist, histo, 150, -1.5, 1.5 );
  sprintf(histo, "PFJet_Reso_hcal");
  book1D(v_hist, histo, 150, -1.5, 1.5 );
  sprintf(histo, "PFJet_Reso_base");
  book1D(v_hist, histo, 150, -1.5, 1.5 );
  
  sprintf(histo, "GenJetMET");
  book1D(v_hist, histo, 300, 0., 300.);
  sprintf(histo, "PFJetMET");
  book1D(v_hist, histo, 300, 0., 300.);

  sprintf(histo, "GenJetEta");
  book1D(v_hist, histo, 60, -6., 6.  );
  sprintf(histo, "PFJetEta");
  book1D(v_hist, histo, 60, -6., 6.  );

  sprintf(histo, "n_neutral_jets");
  book1D(v_hist, histo, 100, 0., 100.  );


  sprintf(histo, "PFJetchargedHadronEFraction");
  book1D(v_hist, histo, 110, 0., 1.1 );
  sprintf(histo, "PFJetneutralEFraction");
  book1D(v_hist, histo, 110, 0., 1.1 );
  sprintf(histo, "PFJetchargedEMEFraction");
  book1D(v_hist, histo, 110, 0., 1.1 );
  sprintf(histo, "PFJetMuonEFraction");
  book1D(v_hist, histo, 110, 0., 1.1 );
  sprintf(histo, "PFJetneutralEMEFraction");
  book1D(v_hist, histo, 110, 0., 1.1 );
  sprintf(histo, "PFJetHFEMEFraction");
  book1D(v_hist, histo, 110, 0., 1.1 );
  sprintf(histo, "PFJetchargedHFHadronEFraction");
  book1D(v_hist, histo, 110, 0., 1.1 );
  sprintf(histo, "PFtotalneutralEFraction");
  book1D(v_hist, histo, 110, 0., 1.1  );
  
}
//
// relabel 1DProf histograms
//
void relabelProfA(TList *v_hist, std::string name)
{
  TProfile* htemp = (TProfile*) v_hist->FindObject(name.c_str());
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
}

void relabel2D(TList *v_hist, std::string name) // added by Bryan 
{
  TH2D* htemp = (TH2D*) v_hist->FindObject(name.c_str());
  htemp->GetXaxis()->SetBinLabel(1,"Charged Had");
  htemp->GetXaxis()->SetBinLabel(2,"Electron");
  htemp->GetXaxis()->SetBinLabel(4,"Photon");
  htemp->GetXaxis()->SetBinLabel(5,"Nuetral Had");
}
//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value)
{
  TH1F* htemp = (TH1F*) v_hist->FindObject(name.c_str());
  htemp->Fill(value);
}
void fill1D(TList *v_hist, std::string name, double value, double valuey)
{
  TH1F* htemp = (TH1F*) v_hist->FindObject(name.c_str());
  htemp->Fill(value, valuey);
}
//
// Fill 1D Profile histograms
//
void fill1DProf(TList *v_hist, std::string name, double value, double valuey)
{
  TProfile* htemp = (TProfile*) v_hist->FindObject(name.c_str());
  htemp->Fill(value,valuey);
}
// 
// Fill 2D histograms
//
void fill2D(TList *v_hist, std::string name, double valuex, double valuey)
{
  TH2D* h_temp = (TH2D*) v_hist->FindObject(name.c_str());
  h_temp->Fill(valuex,valuey);
}
