#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "/home/audrey/.local/nptool/default/include/EpicData.h"
//#include "/home/epic/.local/nptool/default/include/EpicData.h"

using namespace std;

void run(int runid=17)
{

  // === =========================================================
  // === variables 
  
  char   name[100];
  double mHF_time_curr = -1; // for each HF
  double thf_prev = -1;      // for each HF without reflexion
  double thf_curr = -1;      // for each HF without reflexion
  double thf_time_curr = -1; // for each HF without reflexion with DT=2.5ms
  
  struct FCdata{
  	short  mFC_anode;
  	double mFC_time;
  	double mFC_qmax;
  	double mFC_q1;
  	double mFC_q2;
  	double mFC_q3;
  	bool   mFC_alpha;
  	bool   mFC_fission;
  	double mHF_curr;
  };
  vector<FCdata> pendingFC;
  
  // === for output tree
  
  short  fFC_anode;
  double fFC_time;
  double fFC_qmax;
  double fFC_q1;
  double fFC_q2;
  double fFC_q3;
  bool   fFC_alpha;
  bool   fFC_fission;
  double fHF_time_curr;
  double fHF_time_next;
  
  // === output tree
  TFile * fout = new TFile(Form("../../output/pretreat/pretreat_run%i.root",runid),"recreate");
  TTree * tout = new TTree("EpicPretreatData","no cross talk, no reflexion, HF_curr, HF_next");
  tout->Branch("fFC_anode",&fFC_anode);
  tout->Branch("fFC_time",&fFC_time);
  tout->Branch("fFC_qmax",&fFC_qmax);
  tout->Branch("fFC_q1",&fFC_q1);
  tout->Branch("fFC_q2",&fFC_q2);
  tout->Branch("fFC_q3",&fFC_q3);
  tout->Branch("fFC_alpha",&fFC_alpha);
  tout->Branch("fFC_fission",&fFC_fission);
  tout->Branch("fHF_time_curr",&fHF_time_curr);
  tout->Branch("fHF_time_next",&fHF_time_next);

  // === =========================================================
  // === TCutG for alpha and fission
  TFile * fcutg = new TFile("2DdiscriTCutG.root","read");
  TCutG * cutalpha[10];
  TCutG * cutfission[10];
  for(short a = 0 ; a < 10 ; a++){
    cutalpha[a] = (TCutG*)fcutg->Get(Form("A%02d_alpha",a+1));
    cutfission[a] = (TCutG*)fcutg->Get(Form("A%02d_fission",a+1));
  } 
  fcutg->Close();

  // === =========================================================
  // === histograms to understand the HF 
  
  TH1D * h1_THF = new TH1D("THF","THF",86400,0,86400);
  TH1D * h1_DeltaTHF = new TH1D("DeltaTHF","DeltaTHF",5100,0,5100000);
  TH1D * h1_DeltaTHFzoom = new TH1D("DeltaTHFzoom","DeltaTHFzoom",20000,2499000,2501000);
  TH1D * h1_onlyR_THF = new TH1D("THF_onlyReflex","THF_onlyReflex",86400,0,86400);
  TH1D * h1_onlyR_DeltaTHFzoom = new TH1D("DeltaTHFzoom_onlyReflex","DeltaTHFzoom_onlyReflex",10000,850,950);
  TH1D * h1_noR_THF = new TH1D("THF_noReflex","THF_noReflex",86400,0,86400);
  TH1D * h1_noR_DeltaTHF = new TH1D("DeltaTHF_noReflex","DeltaTHF_noReflex",5100,0,5100000);
  TH1D * h1_noR_DeltaTHFzoom = new TH1D("DeltaTHFzoom_noReflex","DeltaTHFzoom_noReflex",20000,2499000,2501000);
  TH1D * h1_noR_Good_THF = new TH1D("THF_noReflex_GoodDeltaT","THF_noReflex_GoodDeltaT",86400,0,86400);
  
  h1_THF->GetXaxis()->SetTitle("Time [s]");
  h1_DeltaTHF->GetXaxis()->SetTitle("Time [ns]");
  h1_DeltaTHFzoom->GetXaxis()->SetTitle("Time [ns], 100 ps/bin");

  h1_onlyR_THF->SetLineColor(kBlack);
  h1_onlyR_DeltaTHFzoom->SetLineColor(kBlack);
  
  h1_noR_THF->SetLineColor(kRed);
  h1_noR_DeltaTHF->SetLineColor(kRed);
  h1_noR_DeltaTHFzoom->SetLineColor(kRed);
  
  h1_noR_Good_THF->SetLineColor(8);
  h1_noR_Good_THF->SetLineStyle(2);
  
  h1_THF->SetDirectory(0);
  h1_DeltaTHF->SetDirectory(0);
  h1_DeltaTHFzoom->SetDirectory(0);
  h1_onlyR_THF->SetDirectory(0);
  h1_onlyR_DeltaTHFzoom->SetDirectory(0);
  h1_noR_THF->SetDirectory(0);
  h1_noR_DeltaTHF->SetDirectory(0);
  h1_noR_DeltaTHFzoom->SetDirectory(0);
  h1_noR_Good_THF->SetDirectory(0);
  
  // === =========================================================
  // === 2D-discri
  TH2D * h2_discri[10];
  for(short a = 0 ; a < 10 ; a++){
    h2_discri[a] = new TH2D(Form("discri_A%02d",a+1),Form("discri_A%02d",a+1),3000,0,300000,300,0,3);
    h2_discri[a]->SetDirectory(0);
  } 

  // === =========================================================
  // === input data 
  
  TFile * f = new TFile(Form("../../output/conversion/raw_run%i.root",runid),"read");
  TTree * tFC = (TTree*)f->Get("EpicRawTree");
  //TChain * tFC = new TChain("EpicRawTree");
  // first part with V4B
  //tFC->Add("../../output/conversion/raw_run17.root");
  //tFC->Add("../../output/conversion/raw_run18.root");
  //tFC->Add("../../output/conversion/raw_run19.root");
  //tFC->Add("../../output/conversion/raw_run21.root");
  //tFC->Add("../../output/conversion/raw_run22.root");
  //tFC->Add("../../output/conversion/raw_run23.root");
  //tFC->Add("../../output/conversion/raw_run24.root");
  //tFC->Add("../../output/conversion/raw_run25.root");
  //tFC->Add("../../output/conversion/raw_run26.root");
  //tFC->Add("../../output/conversion/raw_run27.root");
  //tFC->Add("../../output/conversion/raw_run28.root");
  //tFC->Add("../../output/conversion/raw_run29.root");
  //tFC->Add("../../output/conversion/raw_run30.root");
  //tFC->Add("../../output/conversion/raw_run31.root");
  // second part with V4B and new parameters for HF channel
  //tFC->Add("../../output/conversion/raw_run32.root");
  //tFC->Add("../../output/conversion/raw_run33.root");
  //tFC->Add("../../output/conversion/raw_run34.root");
  //tFC->Add("../../output/conversion/raw_run35.root");
  //tFC->Add("../../output/conversion/raw_run36.root");
  //tFC->Add("../../output/conversion/raw_run37.root");
  //tFC->Add("../../output/conversion/raw_run38.root");
  //tFC->Add("../../output/conversion/raw_run39.root");
  //tFC->Add("../../output/conversion/raw_run40.root");
  
  // === =========================================================
  // === branches 
  
  epic::EpicData *epicFC = nullptr;
  int statusFC = tFC->SetBranchAddress("epic",&epicFC);


  // === =========================================================
  // === loop 
  ULong64_t nentries = tFC->GetEntries();
  cout << "nentries = " << nentries << endl;
  for(ULong64_t entry=0; entry < nentries ; entry++){
  
    if ((entry % 1000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;
    int bytes = tFC->GetEntry(entry);
    if(!epicFC) continue;
    
    int mult = epicFC->GetFCMult() ;
    
    // --- skip empty entry 
    if (mult==0) continue;
    
    // --- is there an HF in FCdata?
    short ihf  = epicFC->GetHFIndex();
    
    // --- is there an FC in FCdata?
    short ifc  = epicFC->GetQmaxIndex();
    
    // --- init data for TTree
    fFC_anode = -1;
    fFC_time = -1;
    fFC_qmax = -1;
    fFC_q1 = -1;
    fFC_q2 = -1;
    fFC_q3 = -1;
    fFC_alpha = false;
    fFC_fission = false;
    fHF_time_curr = -1000000.;
    fHF_time_next = -1000000.;
   
    // --- HFdata only
    if(ihf >= 0 && ifc == -1){

      // all HF signals
      mHF_time_curr = epicFC->GetTimeHF();
      h1_THF->Fill(mHF_time_curr*1.e-09);
      h1_DeltaTHF->Fill(epicFC->GetDeltaTHF()); 
      h1_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 

      // [noR] : remove reflexion < 1000 ns
      if(epicFC->GetDeltaTHF()>1000){
        thf_prev = thf_curr;
        thf_curr = mHF_time_curr;
        h1_noR_THF->Fill(thf_curr*1.e-09);
        h1_noR_DeltaTHF->Fill(thf_curr-thf_prev); 
        h1_noR_DeltaTHFzoom->Fill(thf_curr-thf_prev); 

        // use reflexion to recover good HF
        if(2500889 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500905) thf_curr -= 878.15;

        // [noR_Good] : remove reflexion _ keep only good DeltaT peak
        if(2500005 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500035){
          h1_noR_Good_THF->Fill(thf_curr*1.e-09); 
          thf_time_curr = thf_curr;
          // fill tree
          if(pendingFC.size()>0){
          	for(size_t i = 0 ; i<pendingFC.size() ; i++){
          		FCdata fc = pendingFC[i];
          		fFC_anode = fc.mFC_anode;
          		fFC_time = fc.mFC_time;
          		fFC_qmax = fc.mFC_qmax;
          		fFC_q1 = fc.mFC_q1;
          		fFC_q2 = fc.mFC_q2;
          		fFC_q3 = fc.mFC_q3;
          		fFC_alpha = fc.mFC_alpha;
          		fFC_fission = fc.mFC_fission;
          		fHF_time_curr = fc.mHF_curr;
          		fHF_time_next = thf_curr;
          		tout->Fill();
          	}
          }
        }// if [Good]
      }// end exclusion reflexion
      else{
    	  h1_onlyR_THF->Fill(epicFC->GetTimeHF()*1.e-09);
    	  h1_onlyR_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 
      }
      pendingFC.clear();
    } // end if HF only
   
    // --- FC data only
    else if(ifc >= 0 && ihf == -1){
      FCdata fc ;
      fc.mFC_anode = epicFC->GetAnodeNbr(ifc);
      fc.mFC_time = epicFC->GetTimeFC(ifc);
      fc.mFC_qmax = epicFC->GetQmax(ifc);
      fc.mFC_q1 = epicFC->GetQ1(ifc);
      fc.mFC_q2 = epicFC->GetQ2(ifc);
      fc.mFC_q3 = epicFC->GetQ3(ifc);
      fc.mHF_curr = thf_time_curr;
      fc.mFC_alpha = false;
      fc.mFC_fission = false;
      if(fc.mFC_q3>0 && fc.mFC_anode>0){ 
      	h2_discri[fc.mFC_anode-1]->Fill(fc.mFC_q1,fc.mFC_q2/fc.mFC_q3);
        if(cutalpha[fc.mFC_anode-1]->IsInside(fc.mFC_q1,fc.mFC_q2/fc.mFC_q3)) 
      		fc.mFC_alpha = true;
      	if(cutfission[fc.mFC_anode-1]->IsInside(fc.mFC_q1,fc.mFC_q2/fc.mFC_q3)) 
      		fc.mFC_fission = true; 
      }
      pendingFC.push_back(fc);
    }
   
    // --- FC and HF data in the same group
    else if(ifc >=0 && ihf >= 0){
      // first FC data, then HF data
      // if HF and FC: HF in the slowest neutrons, good HF
      FCdata fc ;
      fc.mFC_anode = epicFC->GetAnodeNbr(ifc);
      fc.mFC_time = epicFC->GetTimeFC(ifc);
      fc.mFC_qmax = epicFC->GetQmax(ifc);
      fc.mFC_q1 = epicFC->GetQ1(ifc);
      fc.mFC_q2 = epicFC->GetQ2(ifc);
      fc.mFC_q3 = epicFC->GetQ3(ifc);
      fc.mHF_curr = mHF_time_curr;
      fc.mFC_alpha = false;
      fc.mFC_fission = false;
      if(fc.mFC_q3>0 && fc.mFC_anode>0){ 
      	h2_discri[fc.mFC_anode-1]->Fill(fc.mFC_q1,fc.mFC_q2/fc.mFC_q3);
        if(cutalpha[fc.mFC_anode-1]->IsInside(fc.mFC_q1,fc.mFC_q2/fc.mFC_q3)) 
      		fc.mFC_alpha = true;
      	if(cutfission[fc.mFC_anode-1]->IsInside(fc.mFC_q1,fc.mFC_q2/fc.mFC_q3)) 
      		fc.mFC_fission = true; 
      }
      pendingFC.push_back(fc);

      // then, HFdata
      // all HF signals
      mHF_time_curr = epicFC->GetTimeHF();
      h1_THF->Fill(mHF_time_curr*1.e-09);
      h1_DeltaTHF->Fill(epicFC->GetDeltaTHF()); 
      h1_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 

      // [noR] : remove reflexion < 1000 ns
      if(epicFC->GetDeltaTHF()>1000){
        thf_prev = thf_curr;
        thf_curr = mHF_time_curr;
        h1_noR_THF->Fill(thf_curr*1.e-09);
        h1_noR_DeltaTHF->Fill(thf_curr-thf_prev); 
        h1_noR_DeltaTHFzoom->Fill(thf_curr-thf_prev); 

        // use reflexion to recover good HF
        if(2500889 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500905) thf_curr -= 878.15;

        // [noR_Good] : remove reflexion _ keep only good DeltaT peak
        if(2500005 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500035){
          h1_noR_Good_THF->Fill(thf_curr*1.e-09); 
          thf_time_curr = thf_curr;
          // fill tree
          if(pendingFC.size()>0){
          	for(size_t i = 0 ; i<pendingFC.size() ; i++){
          		FCdata fc = pendingFC[i];
          		fFC_anode = fc.mFC_anode;
          		fFC_time = fc.mFC_time;
          		fFC_qmax = fc.mFC_qmax;
          		fFC_q1 = fc.mFC_q1;
          		fFC_q2 = fc.mFC_q2;
          		fFC_q3 = fc.mFC_q3;
          		fFC_alpha = fc.mFC_alpha;
          		fFC_fission = fc.mFC_fission;
          		fHF_time_curr = fc.mHF_curr;
          		fHF_time_next = thf_curr;
          		tout->Fill();
          	}
          }
        }// if [Good]
      }// end exclusion reflexion
      else{
    	  h1_onlyR_THF->Fill(epicFC->GetTimeHF()*1.e-09);
    	  h1_onlyR_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 
      }
      pendingFC.clear();
    
    }// end of FC and HF in the same groupoe
    
    else{
      cout << "EVENT NOT TAKEN INTO ACCOUNT: ifc = " << ifc << ", ihf = " << ihf << endl;
    }
  }// end of loop over the entries

  cout << endl;
  
  TCanvas * can1 = new TCanvas("HF","HF",0,0,2000,3000);
  can1->Divide(1,4);
  can1->cd(1); h1_THF->Draw(); h1_noR_THF->Draw("sames"); h1_noR_Good_THF->Draw("sames"); h1_onlyR_THF->Draw("sames"); 
  can1->cd(2); h1_DeltaTHF->Draw(); h1_noR_DeltaTHF->Draw("sames"); 
  can1->cd(3); h1_DeltaTHFzoom->Draw(); h1_noR_DeltaTHFzoom->Draw("sames");
  can1->cd(4); h1_onlyR_DeltaTHFzoom->Draw(); 
  
  TCanvas * can2[10];
  for (short a = 0 ; a <10 ; a++){
  	can2[a] = new TCanvas(Form("Q2Q3vQ1_A%02d",a+1),Form("Q2Q3vQ1_A%02d",a+1),0,0,2000,1500);
  	can2[a]->cd();
  	h2_discri[a]->Draw("col"); cutalpha[a]->Draw("same"); cutfission[a]->Draw("same");
  }
  
  fout->cd();
  tout->Write();
  fout->Close();
}
