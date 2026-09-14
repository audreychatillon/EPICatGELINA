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

#include "../ClassDef/EpicPreTreat.h"

#define PLOT_MEAN_SIGNALS   1 

#define faster_sample_size_ns 2

void run(int run_number, int P, int HV) {

    // ==========================================================================================
    // === VARIABLES

    char   name[100];
    int    Q1max = 300000;
    int    nA = 1;
    short  det;
    short  anode;
    double qmax;
    double q1;
    double q2;
    double q3;
    double q4;
    double tof_raw;
    double t_cfd;
    double t_qmax;
    int    index;
    vector<double> *fSampler_Signal = nullptr;

  // ===========================================================================================
  // === INPUT DATA
  TChain *ch = new TChain("EpicPreTreat");
  ch->Add(Form("../../output/pretreat/pretreat%i_%imbar_%iV.root",run_number,P,HV)); 
  EpicPreTreat raw(ch);


#if PLOT_MEAN_SIGNALS
    int    ncuts = 4;
    map<string,TCutG*> cuts;
    map<string,TH2F*> h2_sig;
    TFile *fcut = new TFile(Form("v4b_tcutg_discri_%imbar_%iV/fsave_tcutg_discri_A04.root",P,HV),"read");
    for(int c = 0 ; c < ncuts ; c++){
      string cutname = Form("A04_cut%02d", c+1);
      TCutG *cut = (TCutG*)fcut->Get(cutname.c_str());
      if(cut){
    	  cuts[cutname] = cut;
    	  cuts[cutname]->ls();
      }
      else     cout << "WARNING : TCutG " << cutname << "not found" << endl;
      string hisname = Form("V4b_SIGNAL_%imbar_%iV_A04_cut%02d",P,HV, c+1);
      h2_sig[hisname] = new TH2F(hisname.c_str(),hisname.c_str(),100,-20,180,2500,-300,4700);
      h2_sig[hisname]->SetDirectory(0);
    }
    fcut->Close();
#endif

  // ===========================================================================================
  // === HISTOGRAMS
  TH2F *h2_Q1_ifQmax;
  vector<TH1F *> h1_Q1_ifQmax(nA);
  vector<TH2F *> h2_DT_vs_Q1_ifQmax(nA);
  vector<TH2F *> h2_DT_vs_Qmax_ifQmax(nA);
  vector<TH2F *> h2_Discri(nA);
  vector<TH1F *> h1_TofRaw_ifQmax(nA);

  sprintf(name, "V4b_P%imbar_%iV_Q1_vs_A_ifQmax", P, HV);
  h2_Q1_ifQmax = new TH2F(name, name, 13, -0.5, 12.5, Q1max / 200, 0, Q1max);
  h2_Q1_ifQmax->SetDirectory(0);

  for (int a = 0; a < nA; a++) {

    sprintf(name, "V4b_P%imbar_%iV_Anode%i_Q1_ifQmax", P, HV, a+1);
    h1_Q1_ifQmax[a] = new TH1F(name, name, Q1max / 100, 0, Q1max);
    h1_Q1_ifQmax[a]->SetLineColor(kBlue);
    h1_Q1_ifQmax[a]->SetDirectory(0);
    h1_Q1_ifQmax[a]->ls();

    sprintf(name, "V4b_P%imbar_%iV_Anode%i_DT_vs_Q1_ifQmax", P, HV, a+1);
    h2_DT_vs_Q1_ifQmax[a] = new TH2F(name, name, Q1max / 200, 0, Q1max, 5000, -100, 400);
    h2_DT_vs_Q1_ifQmax[a]->SetDirectory(0);
    h2_DT_vs_Q1_ifQmax[a]->ls();

    sprintf(name, "V4b_P%imbar_%iV_Anode%i_DT_vs_Qmax_ifQmax", P, HV, a+1);
    h2_DT_vs_Qmax_ifQmax[a] = new TH2F(name, name, 500, 0, 5000, 5000, -100, 400);
    h2_DT_vs_Qmax_ifQmax[a]->SetDirectory(0);
    h2_DT_vs_Qmax_ifQmax[a]->ls();

    sprintf(name, "V4b_P%imbar_%iV_Anode%i_Q2Q3_vs_Q1_ifQmax", P, HV, a+1);
    h2_Discri[a] = new TH2F(name, name, Q1max / 200, 0, Q1max, 5000, 0, 5);
    h2_Discri[a]->SetDirectory(0);
    h2_Discri[a]->ls();

    sprintf(name, "V4b_P%imbar_%iV_Anode%i_TofRaw_ifQmax", P, HV, a+1);
    h1_TofRaw_ifQmax[a] = new TH1F(name, name, 30000,-500,2500);
    h1_TofRaw_ifQmax[a]->SetLineColor(kBlue);
    h1_TofRaw_ifQmax[a]->SetDirectory(0);
    h1_TofRaw_ifQmax[a]->ls();
  } // end of for(a)

  // ===========================================================================================
  // === LOOP
  Long64_t nentries = (Long64_t)ch->GetEntries();
  cout << "nentries = " << nentries << endl;
  for (Long64_t entry = 0; entry < nentries; entry++) {

    ch->GetEntry(entry);
    if ((entry % 1000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;

    det    = raw.fFC_DetNbr;
    if (det!=1)               continue; // keep only FC_1 with V4b

    anode   = raw.fFC_AnodeNbr;
    if (anode!=4) continue;
    qmax    = raw.fFC_Qmax;
    q1      = raw.fFC_Q1;
    q2      = raw.fFC_Q2;
    q3      = raw.fFC_Q3;
    q4      = raw.fFC_Q4;
    t_cfd   = raw.fFC_TimeCfd;
    t_qmax  = raw.fFC_TimeQmax;
    tof_raw = raw.fFC_TofRaw;
    fSampler_Signal = raw.fQmax_Sampler;
    index  = 0;

    if (q3 > 0) {
      h2_Q1_ifQmax->Fill(anode, q1);
      h1_Q1_ifQmax[index]->Fill(q1);
      h2_DT_vs_Q1_ifQmax[index]->Fill(q1, t_qmax - t_cfd);
      h2_DT_vs_Qmax_ifQmax[index]->Fill(qmax, t_qmax - t_cfd);
      h2_Discri[index]->Fill(q1, q2 / q3);
      h1_TofRaw_ifQmax[index]->Fill(tof_raw);
#if PLOT_MEAN_SIGNALS
      for(int c = 1 ; c <= ncuts ; c++){
        string cutname = Form("A04_cut%02d",c);
        if(cuts[cutname]->IsInside(q1,q2/q3)){ 
            string hisname = Form("V4b_SIGNAL_%imbar_%iV_A04_cut%02d",P,HV, c);
            for(int sample=0; sample<(int)fSampler_Signal->size(); sample++)
                h2_sig[hisname]->Fill(2*sample,fSampler_Signal->at(sample));
        }
      }
#endif
    } // end of if(q3)
  } // end of loop over the entries
  cout << endl;

  // ===========================================================================================
  // === DRAW
  TFile *fsave = new TFile(Form("v4b_tcutg_discri_%imbar_%iV/fsave_histosi_A04.root",P,HV),"recreate");

  sprintf(name, "V4b_P%imbar_%iV_Q2Q3vQ1_ifQmax", P, HV);
  TCanvas *can1 = new TCanvas(name, name, 0, 0, 2000, 1500);

  sprintf(name, "V4b_P%imbar_%iV_Q1_ifQmax", P, HV);
  TCanvas *can2 = new TCanvas(name, name, 0, 0, 2000, 1500);

  sprintf(name, "V4b_P%imbar_%iV_DTvQ1_ifQmax", P, HV);
  TCanvas *can3 = new TCanvas(name, name, 0, 0, 2000, 1500);

  sprintf(name, "V4b_P%imbar_%iV_DTvQmax_ifQmax", P, HV);
  TCanvas *can4 = new TCanvas(name, name, 0, 0, 2000, 1500);

  sprintf(name, "V4b_P%imbar_%iV_TofRaw_ifQmax", P, HV);
  TCanvas *can5 = new TCanvas(name, name, 0, 0, 2000, 1500);

#if PLOT_MEAN_SIGNALS
  vector<TCanvas*> can6(nA);
#endif

  for (int a = 0; a < nA; a++) {
    can1->cd();    gPad->SetLogz();  gPad->SetLogx();
    h2_Discri[a]->GetYaxis()->SetRangeUser(0.,2.1);
    h2_Discri[a]->Draw("colz");
    fsave->cd();     h2_Discri[a]->Write();

    can2->cd();    gPad->SetLogy();    gPad->SetLogx();
    h1_Q1_ifQmax[a]->Draw();
    fsave->cd();     h1_Q1_ifQmax[a]->Write();
    
    can3->cd();    gPad->SetLogz();
    h2_DT_vs_Q1_ifQmax[a]->Draw("colz");
    fsave->cd();     h2_DT_vs_Q1_ifQmax[a]->Write();
    
    can4->cd();    gPad->SetLogz();
    h2_DT_vs_Qmax_ifQmax[a]->Draw("colz");
    fsave->cd();     h2_DT_vs_Qmax_ifQmax[a]->Write();
    
    can5->cd();    gPad->SetLogz();
    h1_TofRaw_ifQmax[a]->Draw();
    fsave->cd();     h1_TofRaw_ifQmax[a]->Write();

#if PLOT_MEAN_SIGNALS
    can6[a] = new TCanvas(Form("signal_%imbar_%iV_A04",P,HV),Form("signal_%imbar_%iV_A04",P,HV),2000,1000);
    can6[a]->Divide(2,2);
    for(int c = 1 ; c <= ncuts ; c++){
        string cutname = Form("A04_cut%02d",c);
        can1->cd(); cuts[cutname]->Draw("same");
        string hisname = Form("V4b_SIGNAL_%imbar_%iV_A04_cut%02d",P,HV, c);
        can6[a]->cd(c); gPad->SetGridx(); gPad->SetGridy();
        h2_sig[hisname]->Draw("col");
    	fsave->cd();     h2_sig[hisname]->Write();
    }
    fsave->cd(); can6[a]->Write();
#endif
    fsave->cd(); 
    can1->Write();
    can3->Write();
    can4->Write();
    can5->Write();
  }
  fsave->cd();     h2_Q1_ifQmax->Write(); can2->Write();
  fsave->Close();
}
