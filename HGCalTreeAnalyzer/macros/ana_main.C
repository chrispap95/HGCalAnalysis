#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()
#include <algorithm> 

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
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

// HOW TO USE:
// Change "numFiles" to desired amount of files to run over (max: 20)
// Change "geoType" to the type of geometry 
// outFile names are at the bottom
// change bool statement if you want to run debug file
 
std::vector<std::string> GetInputFiles(TString geoConfig)
{
  int numFiles = 1;
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

void WriteToOutput(int n_f,  TH2F* h[n_f], TString outputFile, TString fileOption)
{
  TFile *f = new TFile(outputFile,fileOption); // output file 
  for ( int i = 1; i <= n_f; i++)
    {
      h[i] -> Write();
    }
  f->Close();
}


void ana_main()
{
  bool test_file = false; // if testing setup with single file (will have to edit below for file choice)

  TString geoType = "D30" ; // "D30" for D30 geo, "D28" for D28
  std::vector<std::string> inputFiles;
  if (!test_file){
      inputFiles = GetInputFiles(geoType); 
  }
  else inputFiles.push_back("ttbar_10_4_D30_pt25.root"); // edit here if test file!!!!!!!!!!!!!!===================

  TChain *ch = new TChain("hgcalTupleTree/tree");
  
  for (unsigned int iFile=0; iFile<inputFiles.size(); ++iFile) {
    ch->Add(inputFiles[iFile].c_str());
    std::cout<<inputFiles[iFile]<<std::endl;
  }
  //ch->Print();
  //   printf("%d;\n",ch->GetNtrees());
  printf("%lld;\n",ch->GetEntries());
  
  TTreeReader     fReader(ch);  //!the tree reader
   
  //
  // Set up TTreeReader's
  // -- use MakeSelector of root
  //
  
  // Readers to access the data (delete the ones you do not need).

  //TTreeReaderArray<double> GenParEta = {fReader, "GenParEta"};
  //TTreeReaderArray<double> GenParM = {fReader, "GenParM"};
  //TTreeReaderArray<double> GenParPhi = {fReader, "GenParPhi"};
  //TTreeReaderArray<double> GenParPt = {fReader, "GenParPt"};
  //TTreeReaderArray<double> GeneralTracksD0 = {fReader, "GeneralTracksD0"};
  //TTreeReaderArray<double> GeneralTracksDZ = {fReader, "GeneralTracksDZ"};
  //TTreeReaderArray<double> GeneralTracksEta = {fReader, "GeneralTracksEta"};
  //TTreeReaderArray<double> GeneralTracksPhi = {fReader, "GeneralTracksPhi"};
  //TTreeReaderArray<double> GeneralTracksPt = {fReader, "GeneralTracksPt"};
  //TTreeReaderArray<float> HBHERecHitEnergy = {fReader, "HBHERecHitEnergy"};
  //TTreeReaderArray<float> HBHERecHitEta = {fReader, "HBHERecHitEta"};
  //TTreeReaderArray<float> HBHERecHitPhi = {fReader, "HBHERecHitPhi"};
  //TTreeReaderArray<float> HBHERecHitTime = {fReader, "HBHERecHitTime"};
  //TTreeReaderArray<float> HGCDigiCharge = {fReader, "HGCDigiCharge"};
  TTreeReaderArray<float> HGCDigiEta = {fReader, "HGCDigiEta"};
  TTreeReaderArray<float> HGCDigiPhi = {fReader, "HGCDigiPhi"};
  TTreeReaderArray<float> HGCDigiPosx = {fReader, "HGCDigiPosx"};
  TTreeReaderArray<float> HGCDigiPosy = {fReader, "HGCDigiPosy"};
  TTreeReaderArray<float> HGCDigiPosz = {fReader, "HGCDigiPosz"};
  //TTreeReaderArray<float> HGCRecHitEnergy = {fReader, "HGCRecHitEnergy"};
  //TTreeReaderArray<float> HGCRecHitEta = {fReader, "HGCRecHitEta"};
  //TTreeReaderArray<float> HGCRecHitPhi = {fReader, "HGCRecHitPhi"};
  TTreeReaderArray<float> HGCRecHitPosx = {fReader, "HGCRecHitPosx"};
  TTreeReaderArray<float> HGCRecHitPosy = {fReader, "HGCRecHitPosy"};
  //TTreeReaderArray<float> HGCRecHitPosz = {fReader, "HGCRecHitPosz"};
  //TTreeReaderArray<float> HGCSimHitsEnergy = {fReader, "HGCSimHitsEnergy"};
  TTreeReaderArray<float> HGCSimHitsEta = {fReader, "HGCSimHitsEta"};
  //TTreeReaderArray<float> HGCSimHitsPhi = {fReader, "HGCSimHitsPhi"};
  TTreeReaderArray<float> HGCSimHitsPosx = {fReader, "HGCSimHitsPosx"};
  TTreeReaderArray<float> HGCSimHitsPosy = {fReader, "HGCSimHitsPosy"};
  //TTreeReaderArray<float> HGCSimHitsPosz = {fReader, "HGCSimHitsPosz"};
  //TTreeReaderArray<float> HGCSimHitsTime = {fReader, "HGCSimHitsTime"};
  //TTreeReaderArray<float> HGCUncalibratedRecHitAmplitude = {fReader, "HGCUncalibratedRecHitAmplitude"};
  //TTreeReaderArray<float> HGCUncalibratedRecHitEta = {fReader, "HGCUncalibratedRecHitEta"};
  //TTreeReaderArray<float> HGCUncalibratedRecHitPhi = {fReader, "HGCUncalibratedRecHitPhi"};
  //TTreeReaderArray<float> HGCUncalibratedRecHitPosx = {fReader, "HGCUncalibratedRecHitPosx"};
  //TTreeReaderArray<float> HGCUncalibratedRecHitPosy = {fReader, "HGCUncalibratedRecHitPosy"};
  //TTreeReaderArray<float> HGCUncalibratedRecHitPosz = {fReader, "HGCUncalibratedRecHitPosz"};
  //TTreeReaderArray<float> SimTracksEta = {fReader, "SimTracksEta"};
  //TTreeReaderArray<float> SimTracksPhi = {fReader, "SimTracksPhi"};
  //TTreeReaderArray<float> SimTracksPt = {fReader, "SimTracksPt"};
  //TTreeReaderArray<float> SimTracksR = {fReader, "SimTracksR"};
  //TTreeReaderArray<float> SimTracksZ = {fReader, "SimTracksZ"};
  //TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
  //TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
  //TTreeReaderArray<int> GeneralTracksNValidHits = {fReader, "GeneralTracksNValidHits"};
  //TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
  //TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
  //TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
  //TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
  //TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
  //TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
  //TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
  TTreeReaderArray<int> HGCDigiCellU = {fReader, "HGCDigiCellU"};
  TTreeReaderArray<int> HGCDigiCellV = {fReader, "HGCDigiCellV"};
  TTreeReaderArray<int> HGCDigiIEta = {fReader, "HGCDigiIEta"};
  TTreeReaderArray<int> HGCDigiIPhi = {fReader, "HGCDigiIPhi"};
  TTreeReaderArray<int> HGCDigiIndex = {fReader, "HGCDigiIndex"};
  TTreeReaderArray<int> HGCDigiLayer = {fReader, "HGCDigiLayer"};
  TTreeReaderArray<int> HGCDigiWaferU = {fReader, "HGCDigiWaferU"};
  TTreeReaderArray<int> HGCDigiWaferV = {fReader, "HGCDigiWaferV"};
  TTreeReaderArray<int> HGCRecHitIndex = {fReader, "HGCRecHitIndex"};
  TTreeReaderArray<int> HGCRecHitLayer = {fReader, "HGCRecHitLayer"};
  TTreeReaderArray<int> HGCSimHitsCellU = {fReader, "HGCSimHitsCellU"};
  TTreeReaderArray<int> HGCSimHitsCellV = {fReader, "HGCSimHitsCellV"};
  //TTreeReaderArray<int> HGCSimHitsIEta = {fReader, "HGCSimHitsIEta"};
  //TTreeReaderArray<int> HGCSimHitsIPhi = {fReader, "HGCSimHitsIPhi"};
  TTreeReaderArray<int> HGCSimHitsIndex = {fReader, "HGCSimHitsIndex"};
  TTreeReaderArray<int> HGCSimHitsLayer = {fReader, "HGCSimHitsLayer"};
  //TTreeReaderArray<int> HGCSimHitsSubdet = {fReader, "HGCSimHitsSubdet"};
  TTreeReaderArray<int> HGCSimHitsWaferU = {fReader, "HGCSimHitsWaferU"};
  TTreeReaderArray<int> HGCSimHitsWaferV = {fReader, "HGCSimHitsWaferV"};
  //TTreeReaderArray<int> HGCUncalibratedRecHitIndex = {fReader, "HGCUncalibratedRecHitIndex"};
  //TTreeReaderArray<int> HGCUncalibratedRecHitLayer = {fReader, "HGCUncalibratedRecHitLayer"};
  //TTreeReaderArray<int> SimTracksCharge = {fReader, "SimTracksCharge"};
  //TTreeReaderArray<int> SimTracksPID = {fReader, "SimTracksPID"};
  TTreeReaderValue<UInt_t> bx = {fReader, "bx"};
  TTreeReaderValue<UInt_t> event = {fReader, "event"};
  TTreeReaderValue<UInt_t> ls = {fReader, "ls"};
  TTreeReaderValue<UInt_t> orbit = {fReader, "orbit"};
  TTreeReaderValue<UInt_t> run = {fReader, "run"};
  TTreeReaderArray<unsigned short> HGCDigiADC = {fReader, "HGCDigiADC"};

  // Histos defined here

  TH2F *h_digi_cee_waferuv[28];
  for(int i=1; i<=28; i++) 
    {
      std::ostringstream name;
      name <<"h_digi_cee_waferuv_"<<i;
      h_digi_cee_waferuv[i] = new TH2F(name.str().c_str(),name.str().c_str(),31,-15.5,15.5,31,-15.5,15.5);
    }

  TH2F *h_digi_ceh_waferuv[24];
  for(int i=1; i<=24; i++) 
    {
      std::ostringstream name;
      name <<"h_digi_ceh_waferuv_"<<i;
      h_digi_ceh_waferuv[i] = new TH2F(name.str().c_str(),name.str().c_str(),31,-15.5,15.5,31,-15.5,15.5);
    }

    TH2F *h_digi_cee_celluv[28];
  for(int i=1; i<=28; i++) 
    {
      std::ostringstream name;
      name <<"h_digi_cee_celluv_"<<i;
      h_digi_cee_celluv[i] = new TH2F(name.str().c_str(),name.str().c_str(),36,-5.5,30.5,36,-5.5,30.5);
    }

  TH2F *h_digi_ceh_celluv[24];
  for(int i=1; i<=24; i++) 
    {
      std::ostringstream name;
      name <<"h_digi_ceh_celluv_"<<i;
      h_digi_ceh_celluv[i] = new TH2F(name.str().c_str(),name.str().c_str(),36,-5.5,30.5,36,-5.5,30.5);
    }
  
  TH2F *h_digi_cee_cellxy[28];
  for(int i=1; i<=28; i++) 
    {
      std::ostringstream name;
      name <<"h_digi_cee_cellxy_"<<i;
      h_digi_cee_cellxy[i] = new TH2F(name.str().c_str(),name.str().c_str(),1000,-750,750,750,-750,750);
    }
  TH2F *h_digi_cee_xy[28];
  for(int i=1; i<=28; i++) 
    {
      std::ostringstream name;
      name <<"h_digi_cee_xy_"<<i;
      h_digi_cee_xy[i] = new TH2F(name.str().c_str(),name.str().c_str(),300,-300,300,300,-300,300);
    }

  TH2F *h_digi_ceh_cellxy[24];
  for(int i=1; i<=24; i++) 
    {
      std::ostringstream name;
      name <<"h_digi_ceh_cellxy_"<<i;
      h_digi_ceh_cellxy[i] = new TH2F(name.str().c_str(),name.str().c_str(),1000,-750,750,750,-750,750);
    }

  TH2F *h_sim_cee_cellxy[28];
  for(int i=1; i<=28; i++) 
    {
      std::ostringstream name;
      name <<"h_sim_cee_cellxy_"<<i;
      h_sim_cee_cellxy[i] = new TH2F(name.str().c_str(),name.str().c_str(),1000,-750,750,750,-750,750);
    }

  TH2F *h_sim_ceh_cellxy[24];
  for(int i=1; i<=24; i++) 
    {
      std::ostringstream name;
      name <<"h_sim_ceh_cellxy_"<<i;
      h_sim_ceh_cellxy[i] = new TH2F(name.str().c_str(),name.str().c_str(),1000,-750,750,750,-750,750);
    }
  TH2F *h_sim_cee_xy[28];
  for(int i=1; i<=28; i++) 
    {
      std::ostringstream name;
      name <<"h_sim_cee_xy_"<<i;
      h_sim_cee_xy[i] = new TH2F(name.str().c_str(),name.str().c_str(),300,-300,300,300,-300,300);
    }
  TH2F *h_reco_cee_xy[28];
  for(int i=1; i<=28; i++) 
    {
      std::ostringstream name;
      name <<"h_reco_cee_xy_"<<i;
      h_reco_cee_xy[i] = new TH2F(name.str().c_str(),name.str().c_str(),300,-300,300,300,-300,300);
    }
  // to figure out fine from coarse wafers 
  ofstream csv_file;
  bool openFile = false; // toggle off if you dont want to slow down code 

  if (openFile)
    {
      csv_file.open ("cell_uv_map.csv");
      csv_file << "Layer #, Wafer U, Wafer V, Cell U, Cell V\n";
    }

  // Main Loop Start
  unsigned int nentries = (Int_t)ch->GetEntries();
  int ievent=0; int maxevents = -1; int skipevents = 0;
  while (fReader.Next()) 
    {
  
      // Progress indicator 
      ievent++;
      if(ievent%50==0) cout << "[HGCAL Response analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
      if (maxevents>0 && ievent>maxevents) break;
      if (ievent<=skipevents) continue;
      
      //std::cout << "Simhit size: " << HGCSimHitsEnergy.GetSize() << std::endl;
      //std::cout << "Digi size:   " << HGCDigiADC.GetSize() << std::endl;
      //std::cout << "Rechit size: " << HGCRecHitEnergy.GetSize() << std::endl;

      // Analyze over Simhits, Rechits, and Digis here
      // Index definition: 0 -- CEE, 1 -- CEH Si, 2 -- CEH Scint
      // CEE -- 1-28 layers, CEH Si -- 1-24 layers, CEH Scint -- 9-24 layers

      // ================RECO==========================
      for (int irc = 0, nrc = HGCRecHitPosx.GetSize(); irc <nrc; ++irc){
	  int l_hit = HGCRecHitLayer[irc];
	  if (HGCRecHitIndex[irc] == 0){
	    h_reco_cee_xy[l_hit]->Fill(HGCRecHitPosx[irc],HGCRecHitPosy[irc]);
	  }
      }
      // ================DIGIS=========================
      for (int irc = 0, nrc = HGCDigiEta.GetSize(); irc <nrc; ++irc)
	{
	  int l_hit = HGCDigiLayer[irc];
	  int u_wafer = HGCDigiWaferU[irc];
	  int v_wafer = HGCDigiWaferV[irc];
	  int u_cell = HGCDigiCellU[irc];
	  int v_cell = HGCDigiCellV[irc];	 

	  int num_cell = 12; // 8 for coarse, 12 for fine ================ NEED TO FIGURE OUT WHICH CELLS ARE COARSE AND FINE
	  double x_wafer= -2*u_wafer+v_wafer;
	  double y_wafer= 2*v_wafer;
	  double x_cell= 1.5*(v_cell-num_cell)+1;
	  double y_cell= u_cell-0.5*(num_cell+v_cell);
	  double x_test = x_cell+1.65*num_cell*x_wafer;
	  double y_test = y_cell+0.9*num_cell*y_wafer;
	  
	  std::ostringstream csv_line;

	  if (HGCDigiIndex[irc] == 0)
	    {
	      h_digi_cee_waferuv[l_hit]->Fill(u_wafer,v_wafer);
	      h_digi_cee_celluv[l_hit]->Fill(u_cell,v_cell);
	      h_digi_cee_cellxy[l_hit]->Fill(x_test,y_test);
	      h_digi_cee_xy[l_hit]->Fill(HGCDigiPosx[irc],HGCDigiPosy[irc]);
	      if (openFile)
		{
		  if (u_cell > 15 || ((abs(u_wafer) <= 2 && abs(v_wafer) < 2) || (abs(u_wafer) < 2 && abs(v_wafer) <= 2)))
		    {
		      csv_line << l_hit-1 <<","<<u_wafer<<","<<v_wafer<<","<<u_cell<<","<<v_cell<<"\n";
		      csv_file << csv_line.str().c_str();
		    }
		}
	    }
	  if (HGCDigiIndex[irc] == 1)
	    {
	      h_digi_ceh_waferuv[l_hit]->Fill(u_wafer,v_wafer);
	      h_digi_ceh_celluv[l_hit]->Fill(u_cell,v_cell);
	      //if ( l_hit == 6 && ((u_cell == 0 && v_cell == 11) || (u_cell == 12 && v_cell == 0) || (u_cell == 23 &&v_cell == 11) || (u_cell == 12 &&v_cell == 23) || (u_cell == 23 &&v_cell == 23))){
	      h_digi_ceh_cellxy[l_hit]->Fill(x_test,y_test);
	      //std::cout<<"Wafer U: " << u_wafer<<", Wafer V: "<<v_wafer<<std::endl;

	      if (openFile)
		{
		  if ((u_cell > 15 || ((abs(u_wafer) <= 2 && abs(v_wafer) < 2) || (abs(u_wafer) < 2 && abs(v_wafer) <= 2))) && l_hit <= 5)
		    {
		      csv_line << l_hit+27 <<","<<u_wafer<<","<<v_wafer<<","<<u_cell<<","<<v_cell<<"\n";
		      csv_file << csv_line.str().c_str();
		    }
		}
	      
	    }
	  
	}
      // =====================SIMHITS STUDY==============
      for (int irc = 0, nrc = HGCSimHitsEta.GetSize(); irc <nrc; ++irc)
	{
	  int l_hit = HGCSimHitsLayer[irc];
	  int u_wafer = HGCSimHitsWaferU[irc];
	  int v_wafer = HGCSimHitsWaferV[irc];
	  int u_cell = HGCSimHitsCellU[irc];
	  int v_cell = HGCSimHitsCellV[irc];	 

	  int num_cell = 12; // 8 for coarse, 12 for fine ================ NEED TO FIGURE OUT WHICH CELLS ARE COARSE AND FINE
	  double x_wafer= -2*u_wafer+v_wafer;
	  double y_wafer= 2*v_wafer;
	  double x_cell= 1.5*(v_cell-num_cell)+1;
	  double y_cell= u_cell-0.5*(num_cell+v_cell);
	  double x_test = x_cell+1.65*num_cell*x_wafer;
	  double y_test = y_cell+0.9*num_cell*y_wafer;
	  
	  //std::ostringstream csv_line;

	  if (HGCSimHitsIndex[irc] == 0)
	    {
	      //h_sim_cee_waferuv[l_hit]->Fill(u_wafer,v_wafer);
	      //h_sim_cee_celluv[l_hit]->Fill(u_cell,v_cell);
	      h_sim_cee_cellxy[l_hit]->Fill(x_test,y_test);
	      h_sim_cee_xy[l_hit]->Fill(HGCSimHitsPosx[irc],HGCSimHitsPosy[irc]);
	//if (openFile)
	//	{
	//	  if (u_cell > 15 || ((abs(u_wafer) <= 2 && abs(v_wafer) < 2) || (abs(u_wafer) < 2 && abs(v_wafer) <= 2)))
	//	    {
	//	      csv_line << l_hit-1 <<","<<u_wafer<<","<<v_wafer<<","<<u_cell<<","<<v_cell<<"\n";
	//	      csv_file << csv_line.str().c_str();
	//	    }
	//	}
	    }
	  if (HGCSimHitsIndex[irc] == 1)
	    {
	      //h_sim_ceh_waferuv[l_hit]->Fill(u_wafer,v_wafer);
	      //h_sim_ceh_celluv[l_hit]->Fill(u_cell,v_cell);
	      //if ( l_hit == 6 && ((u_cell == 0 && v_cell == 11) || (u_cell == 12 && v_cell == 0) || (u_cell == 23 &&v_cell == 11) || (u_cell == 12 &&v_cell == 23) || (u_cell == 23 &&v_cell == 23))){
	      h_sim_ceh_cellxy[l_hit]->Fill(x_test,y_test);
	      //std::cout<<"Wafer U: " << u_wafer<<", Wafer V: "<<v_wafer<<std::endl;

	//if (openFile)
	//	{
	//	  if ((u_cell > 15 || ((abs(u_wafer) <= 2 && abs(v_wafer) < 2) || (abs(u_wafer) < 2 && abs(v_wafer) <= 2))) && l_hit <= 5)
	//	    {
	//	      csv_line << l_hit+27 <<","<<u_wafer<<","<<v_wafer<<","<<u_cell<<","<<v_cell<<"\n";
	//	      csv_file << csv_line.str().c_str();
	//	    }
	//	}
	      
	    }
	  
	}
    }
  
  if (openFile) csv_file.close();

  //end of main event loop
  // write histos to output file
  // Note: cee->28, ceh->24
  int cee = 28; int ceh = 24;
  //WriteToOutput(cee,h_digi_cee_waferuv,"waferuv_10_4_D30.root","RECREATE");
  //WriteToOutput(ceh,h_digi_ceh_waferuv,"waferuv_10_4_D30.root","UPDATE");
  //
  //WriteToOutput(cee,h_digi_cee_celluv,"celluv_10_4_D30.root","RECREATE");
  //WriteToOutput(ceh,h_digi_ceh_celluv,"celluv_10_4_D30.root","UPDATE");
  //
  //WriteToOutput(cee,h_digi_cee_cellxy,"cellmap_10_4_D30.root","RECREATE");
  //WriteToOutput(ceh,h_digi_ceh_cellxy,"cellmap_10_4_D30.root","UPDATE");
  //
  //WriteToOutput(cee,h_sim_cee_cellxy,"cellmap_10_4_D30_sim.root","RECREATE");
  //WriteToOutput(ceh,h_sim_ceh_cellxy,"cellmap_10_4_D30_sim.root","UPDATE");
  WriteToOutput(cee,h_digi_cee_xy,"waferxy_10_4_D30.root","RECREATE"); 
  WriteToOutput(cee,h_sim_cee_xy,"waferxy_10_4_D30.root","UPDATE");
  WriteToOutput(cee,h_reco_cee_xy,"waferxy_10_4_D30.root","UPDATE");
  
}



int mian()
{
  ana_main();

  return 0;
}
