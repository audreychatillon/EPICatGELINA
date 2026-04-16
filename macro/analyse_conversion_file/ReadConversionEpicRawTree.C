#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../ClassDef/EpicRawTree.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

void run(int data_set)
{

  // ==========================================================================================
  // === VARIABLES

  char name[100];

  int nA  ;
  int Pmbar;
  int Q1max;
  int HV;
  string PA;

  // ===========================================================================================
  // === INPUT DATA
  vector<double>*  fSampler_Signal   = nullptr;

  TChain * ch = new TChain("EpicRawTree");

  switch (data_set) {
      case 1: // PA = V4b, P=1030mbar, HV=550V
          nA    = 9 ;
          HV    = 550;
          Pmbar = 1030;
          Q1max = 300000;
          PA    = "4b";
          ch->Add(Form("../../output/conversion/Test1_V4b_%imbar_%iV.root",Pmbar,HV));
          ch->Add(Form("../../output/conversion/Test2_V4b_%imbar_%iV.root",Pmbar,HV));
          ch->Add(Form("../../output/conversion/Test5_V4b_%imbar_%iV.root",Pmbar,HV));
          ch->Add(Form("../../output/conversion/Test6_V4b_%imbar_%iV.root",Pmbar,HV));
          break;
     case 2: // PA = V4b, P=1190mbar, HV = 610V
          nA    = 9 ;
          HV    = 610;
          Pmbar = 1190;
          Q1max = 300000;
          PA    = "4b";
          ch->Add(Form("../../output/conversion/Test3_V4b_%imbar_%iV.root",Pmbar,HV));
          break;
     case 3: // PA = V4b, P=1190mbar, HV = 500V
          nA    = 9 ;
          HV    = 500;
          Pmbar = 1190;
          Q1max = 300000;
          PA    = "4b";
          ch->Add(Form("../../output/conversion/Test4_V4b_%imbar_%iV.root",Pmbar,HV));
          break;
     case 4: // PA = V8.2, P=1380mbar, HV = 650V
          nA    = 4 ;
          HV    = 650;
          Pmbar = 1380;
          Q1max = 400000;
          PA    = "8-2";
          ch->Add(Form("../../output/conversion/Test9_V8-2_%imbar_%iV.root",Pmbar,HV));
          ch->Add(Form("../../output/conversion/Test13_V8-2_%imbar_%iV.root",Pmbar,HV));
          break;
     case 5: // PA = V8.2, P=1040mbar, HV = 550V
          nA    = 4 ;
          HV    = 550;
          Pmbar = 1040;
          Q1max = 400000;
          PA    = "8-2";
          ch->Add(Form("../../output/conversion/Test10_V8-2_%imbar_%iV.root",Pmbar,HV));
          ch->Add(Form("../../output/conversion/Test12_V8-2_%imbar_%iV.root",Pmbar,HV));
          break;
     case 6: // PA = V8.2, P=1210mbar, HV = 680V
          nA    = 4 ;
          HV    = 680;
          Pmbar = 1210;
          Q1max = 500000;
          PA    = "8-2";
          ch->Add(Form("../../output/conversion/Test14_V8-2_%imbar_%iV.root",Pmbar,HV));
          break;
     case 7: // PA = V8.2, P=1370mbar, HV = 720V
          nA    = 4 ;
          HV    = 720;
          Pmbar = 1370;
          Q1max = 500000;
          PA    = "8-2";
          ch->Add(Form("../../output/conversion/Test15_V8-2_%imbar_%iV.root",Pmbar,HV));
          ch->Add(Form("../../output/conversion/Test16_V8-2_%imbar_%iV.root",Pmbar,HV));
          break;
     case 8: // PA = V8.2, P=1075mbar, HV = 720V
          nA    = 4 ;
          HV    = 720;
          Pmbar = 1075;
          Q1max = 500000;
          PA    = "8-2";
          ch->Add(Form("../../output/conversion/Test17_V8-2_%imbar_%iV.root",Pmbar,HV));
          break;
     default:
          nA    = 0 ;
          HV    = 0 ;
          Pmbar = 0 ;
          Q1max = 0 ;
          PA    = "V";
          break;
  }

  unordered_map<int,int> mapindex;
  unordered_map<int,int> mapanode;
  if(nA==9){
      mapindex = { {1,0}, {2,1}, {3,2}, {4,3}, {6,4}, {8,5}, {9,6}, {10,7}, {11,8} };
  }
  else if(nA==4){
      mapindex = { {3,0}, {4,1}, {9,2}, {10,3} };
  }
  for (const auto& p : mapindex) mapanode[p.second] = p.first;
  int ncol = ceil(nA*0.5); 
  cout << "ncol: " << ncol << endl;
 
  EpicRawTree raw(ch);


  // === =========================================================
  // === histograms
  TH1F * h1_TimeHF = new TH1F("TimeHF","TimeHF",86400,0,86400);;
  TH1F * h1_DeltaTimeHF = new TH1F("DeltaTimeHF","DeltaTimeHF",10000,2,3);;
  TH1F * h1_Mult;
  TH2F * h2_MvT;
  TH1F * h1_Q1[nA];
  TH2F * h2_Q1vT[nA];
  TH1F * h1_Q2[nA];
  TH1F * h1_Q3[nA];
  TH2F * h2_Q2Q3vQ1[nA];
  //TH1F * h1_inTofRaw[nA];
  //TH1F * h1_inTofRaw_GammaPeak[nA];

  sprintf(name,"Mult");
  h1_Mult = new TH1F(name,name,13,-0.5,12.5);
  h1_Mult->SetDirectory(0);

  sprintf(name,"Mult_vs_T");
  h2_MvT = new TH2F(name,name,2400,0,72000,13,-0.5,12.5);
  h2_MvT->SetDirectory(0);

  for(unsigned short a = 0 ; a < nA; a++){
     sprintf(name,"Q1vT_A%i_%i",mapanode[a],data_set);
     h2_Q1vT[a] = new TH2F(name,name,2400,0,72000,5000,0,500000);
     h2_Q1vT[a]->GetXaxis()->SetTitle("Time [s] 30s/bin");
     h2_Q1vT[a]->SetDirectory(0);

     sprintf(name,"Q1_A%i_%i",mapanode[a],data_set);
     h1_Q1[a] = new TH1F(name,name,50000,0,500000);
     h1_Q1[a]->SetDirectory(0);

     sprintf(name,"Q2_A%i_%i",mapanode[a],data_set);
     h1_Q2[a] = new TH1F(name,name,50000,0,500000);
     h1_Q2[a]->SetLineColor(8);
     h1_Q2[a]->SetDirectory(0);

     sprintf(name,"Q3_A%i_%i",mapanode[a],data_set);
     h1_Q3[a] = new TH1F(name,name,50000,0,500000);
     h1_Q3[a]->SetLineColor(kCyan);
     h1_Q3[a]->SetDirectory(0);

     sprintf(name,"discri_A%i_%i",mapanode[a],data_set);
     h2_Q2Q3vQ1[a] = new TH2F(name,name,2500,0,500000,500,0,5);
     h2_Q2Q3vQ1[a]->SetDirectory(0);

     //sprintf(name,"inTofRaw_A%i",mapanode[a]);
     //h1_inTofRaw[a] = new TH1F(name,name,60000,-20000,40000);

     //sprintf(name,"inTofRaw_GammaPeak_A%i",mapanode[a]);
     //h1_inTofRaw_GammaPeak[a] = new TH1F(name,name,10000,500,1500);
  }

  
  // ===========================================================================================
  // === LOOP
  double Toff_runs = 0;
  ULong64_t nentries = (ULong64_t)ch->GetEntries();
  cout << "number of entries: " << nentries << endl;
  for(ULong64_t Entry=0; Entry<nentries; Entry++){
    raw.GetEntry(Entry);
    if ((Entry % 1000000)==0) cout << "\r === Entry = " << Entry << " === " << flush;
 
    double Qmax = 0.;
    short  Imax = -1;
    int    fFC_size = (int)raw.fFC_AnodeNbr.size();

    // HF dat
    if(raw.fQmax_Index == -1){
       h1_TimeHF->Fill(1.e-09*raw.fHF_Time); 
       h1_DeltaTimeHF->Fill(1.e-06*(raw.fHF_Time-raw.fHF_TimePrev)); 
    } 
    // get the channel with Qmax
    else if (fFC_size>0){
        Imax = raw.fQmax_Index;
	// read the channel with Qmax
        int anode = raw.fFC_AnodeNbr[Imax];
        double q1 = raw.fFC_Q1[Imax];
        double q2 = raw.fFC_Q2[Imax];
        double q3 = raw.fFC_Q3[Imax];
        double t_s = (raw.fFC_Time[Imax]+Toff_runs) * 1.e-09 ;
        
        int index = mapindex[anode];
        h1_Mult->Fill(fFC_size);
        h2_MvT->Fill(t_s,fFC_size);
        h1_Q1[index]->Fill(q1);
        h2_Q1vT[index]->Fill(t_s,q1);
        h1_Q2[index]->Fill(q2);
        h1_Q3[index]->Fill(q3);
        //h1_inTofRaw[index]->Fill(raw.fFC_TofRaw[i]);
        //h1_inTofRaw_GammaPeak[index]->Fill(raw.fFC_TofRaw[i]);
        if (q3>0) h2_Q2Q3vQ1[index]->Fill(q1,q2/q3);
    }// end of if mult > 0
  }//end of loop over the entries 

  cout << endl;

  TCanvas * can_T0 = new TCanvas("T0","T0",0,0,2500,1500);
  TCanvas * can_mult = new TCanvas("Mult","Mult",0,0,2500,1500);
  TCanvas * can_Q;
  TCanvas * can_Q1vT;
  TCanvas * can_discri;
  //TCanvas * can_tof_raw[nA];
  //TCanvas * can_gamma_peak[nA];
  
  can_T0->Divide(1,2);
  can_T0->cd(1); h1_TimeHF->Draw();
  can_T0->cd(2); h1_DeltaTimeHF->Draw(); 
  unsigned int effective_beam_time_s = 0;
  for (int b = 1 ; b <= h1_TimeHF->GetNbinsX(); b++){
	if(h1_TimeHF->GetBinContent(b) != 0) effective_beam_time_s++;
  }
  cout << "EFFECTIVE BEAM TIME = " << effective_beam_time_s << " s" << endl;
  cout << "NUMBER OF HF: = " << h1_DeltaTimeHF->Integral() << " " << endl;

  sprintf(name,"results/Histo_V%s_%imbar_%iV_%i.root",PA.c_str(),Pmbar,HV,data_set);
  TFile * fsave = new TFile(name,"recreate");
  fsave->cd();
  h1_TimeHF->Write();
  h1_DeltaTimeHF->Write();

  can_mult->Divide(2);
  can_mult->cd(1); h1_Mult->Draw();
  can_mult->cd(2); h2_MvT->Draw("colz");
  fsave->cd(); 
  h1_Mult->Write();
  h2_MvT->Write();

  sprintf(name,"Q_V%s_%imbar_%iV",PA.c_str(),Pmbar,HV);
  can_Q = new TCanvas(name,name,0,0,2500,1500);
  can_Q->Divide(ncol,2);   

  sprintf(name,"Q1vT_V%s_%imbar_%iV",PA.c_str(),Pmbar,HV);
  can_Q1vT = new TCanvas(name,name,0,0,2500,1500);
  can_Q1vT->Divide(ncol,2);   

  sprintf(name,"DISCRI");
  can_discri = new TCanvas(name,name,0,0,2500,1500);
  can_discri->Divide(ncol,2);   
 

  for(unsigned short a = 0 ; a < nA; a++){

     can_Q->cd(a+1);
     h1_Q1[a]->Draw();  h1_Q2[a]->Draw("same"); h1_Q3[a]->Draw("same");

     can_Q1vT->cd(a+1);
     h2_Q1vT[a]->Draw("colz");

     can_discri->cd(a+1); gPad->SetLogx();
     h2_Q2Q3vQ1[a]->Draw("colz");

     //sprintf(name,"inTofRaw_ANODE%i",a+1);
     //can_tof_raw[a] = new TCanvas(name,name,0,0,1500,1000);
     //can_tof_raw[a]->cd();
     //h1_inTofRaw[a]->Draw();

     //sprintf(name,"GammaPeak_ANODE%i",a+1);
     //can_gamma_peak[anode] = new TCanvas(name,name,0,0,1500,1000);
     //can_gamma_peak[anode]->cd();
     //h1_inTofRaw_GammaPeak[anode]->Draw();
     //
     fsave->cd(); 
     h1_Q1[a]->Write();  
     h1_Q2[a]->Write(); 
     h1_Q3[a]->Write();
     h2_Q1vT[a]->Write();
     h2_Q2Q3vQ1[a]->Write();
  }

  fsave->Close();
}
