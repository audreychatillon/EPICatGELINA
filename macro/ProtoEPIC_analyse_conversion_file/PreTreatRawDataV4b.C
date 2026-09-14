#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <unordered_map>

#include <Riostream.h>
#include <TROOT.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TTree.h>

#include "../ClassDef/EpicRawTree.h"

void run(int run_number, int P, int HV) {

  // ==========================================================================================
  // === VARIABLES

  char   name[100];
  int    Q1max;
  short  det;
  short  anode;
  double tof_raw;
  double t_cfd;
  double t_qmax;
  double qmax;
  double q1;
  double q2;
  double q3;
  double q4;
  vector<double> fQmax_Sampler;

  int    index;
  int    i_qmax;
  double q_qmax;

  TFile * fout = new TFile(Form("../../output/pretreat/pretreat%i_%imbar_%iV.root",run_number,P,HV),"recreate"); 
  TTree * tout = new TTree("EpicPreTreat",Form("pretreated raw data, cross talk removed run%i P=%imbar, HV=%iV",run_number,P,HV));
  tout->Branch("fFC_DetNbr",&det);
  tout->Branch("fFC_AnodeNbr",&anode);
  tout->Branch("fFC_TofRaw",&tof_raw);
  tout->Branch("fFC_TimeCfd",&t_cfd);
  tout->Branch("fFC_TimeQmax",&t_qmax);
  tout->Branch("fFC_Qmax",&qmax);
  tout->Branch("fFC_Q1",&q1);
  tout->Branch("fFC_Q2",&q2);
  tout->Branch("fFC_Q3",&q3);
  tout->Branch("fFC_Q4",&q4);
  tout->Branch("fQmax_Sampler",&fQmax_Sampler);


  // ===========================================================================================
  // === INPUT DATA
  TChain *ch = new TChain("EpicRawTree");
  ch->Add(Form("../../output/conversion/Test%i_V4b_%imbar_%iV.root",run_number,P,HV)); // in raw data: 1FC installed
  EpicRawTree raw(ch);


  // ===========================================================================================
  // === LOOP
  Long64_t nentries = (Long64_t)ch->GetEntries();
  cout << "nentries = " << nentries << endl;
  for (Long64_t entry = 0; entry < nentries; entry++) {

    ch->GetEntry(entry);
    if ((entry % 5000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;

    int fFC_size = (int)raw.fFC_AnodeNbr.size();
    if (fFC_size <= 0)        continue;
    if (raw.fQmax_Index < 0)  continue; // exclude fHF data 
    i_qmax = raw.fQmax_Index;
    det    = raw.fFC_DetNbr[i_qmax];
    anode   = raw.fFC_AnodeNbr[i_qmax];
    qmax    = raw.fFC_Qmax[i_qmax];
    q1      = raw.fFC_Q1[i_qmax];
    q2      = raw.fFC_Q2[i_qmax];
    q3      = raw.fFC_Q3[i_qmax];
    q4      = raw.fFC_Q4[i_qmax];
    t_cfd   = raw.fFC_TimeCfd[i_qmax];
    t_qmax  = raw.fFC_TimeQmax[i_qmax];
    tof_raw = raw.fFC_TofRaw[i_qmax];
    fQmax_Sampler = raw.fQmax_Sampler;

    if (q3 > 0) {
	tout->Fill();
    } // end of if(q3)
    
    fQmax_Sampler.clear();
  } // end of loop over the entries
  cout << endl;

  fout->cd();
  tout->Write();
}
