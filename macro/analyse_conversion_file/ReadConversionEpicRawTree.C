#include <iostream>
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

void run(UShort_t run_number, int Pmbar)
{

  char name[100];

  // TODO remove hard coding
  const unsigned short m_nDets = 1;
  unsigned short m_nAnodes[m_nDets] = {11};  
  const unsigned short m_nAnodesTot = 11;


  // === =========================================================
  // === histograms
  TH1F * h1_Q1[m_nAnodesTot];
  TH2F * h2_Q1vT[m_nAnodesTot];
  TH1F * h1_Q2[m_nAnodesTot];
  TH1F * h1_Q3[m_nAnodesTot];
  TH2F * h2_Q2Q3vQ1[m_nAnodesTot];
  TH1F * h1_inTofRaw[m_nAnodesTot];
  TH1F * h1_inTofRaw_GammaPeak[m_nAnodesTot];
  int anode = 0;
  for(unsigned short d = 0 ; d < m_nDets; d++){
    for(unsigned short a = 0 ; a < m_nAnodes[d]; a++){
       sprintf(name,"Q1vT_EPIC%i_A%i",d+1,a+1);
       h2_Q1vT[anode] = new TH2F(name,name,1440,0,86400,5000,0,500000);

       sprintf(name,"Q1_EPIC%i_A%i",d+1,a+1);
       h1_Q1[anode] = new TH1F(name,name,50000,0,500000);

       sprintf(name,"Q2_EPIC%i_A%i",d+1,a+1);
       h1_Q2[anode] = new TH1F(name,name,50000,0,500000);
       h1_Q2[anode]->SetLineColor(8);

       sprintf(name,"Q3_EPIC%i_A%i",d+1,a+1);
       h1_Q3[anode] = new TH1F(name,name,50000,0,500000);
       h1_Q3[anode]->SetLineColor(kCyan);

       sprintf(name,"discri_EPIC%i_A%i",d+1,a+1);
       h2_Q2Q3vQ1[anode] = new TH2F(name,name,2500,0,500000,500,0,10);

       sprintf(name,"inTofRaw_EPIC%i_A%i",d+1,a+1);
       h1_inTofRaw[anode] = new TH1F(name,name,60000,-20000,40000);

       sprintf(name,"inTofRaw_GammaPeak_EPIC%i_A%i",d+1,a+1);
       h1_inTofRaw_GammaPeak[anode] = new TH1F(name,name,10000,500,1500);

       anode++;
    }
  }

  // === =========================================================
  // === input data 
  TChain * ch = new TChain("EpicRawTree");
  sprintf(name,"../../output/conversion/Test%i_V4b_%imbar_HV610V.root",run_number,Pmbar);
  ch->Add(name);
  ch->ls();
  EpicRawTree raw(ch);
  ULong64_t nentries = (ULong64_t)ch->GetEntries();
  cout << "number of entries: " << nentries << endl;
  for(ULong64_t Entry=0; Entry<nentries; Entry++){
    raw.GetEntry(Entry);
    if ((Entry % 5000000)==0) cout << "\r === Entry = " << Entry << " === " << flush;

    int mult = raw.fFC_DetNbr.size();
    for(int i=0; i<mult; i++){
	int det = raw.fFC_DetNbr[i];
	int anode = raw.fFC_AnodeNbr[i];
        double q1 = raw.fFC_Q1[i];
        double q2 = raw.fFC_Q2[i];
        double q3 = raw.fFC_Q3[i];
        double t_s = raw.fFC_Time[i] * 1.e-09 ;
        int index = 0;
        for(int d=0; d<det; d++){
	  index += d*m_nAnodes[d];
	}
        index += anode - 1;
        h1_Q1[index]->Fill(q1);
        h2_Q1vT[index]->Fill(t_s,q1);
        h1_Q2[index]->Fill(q2);
        h1_Q3[index]->Fill(q3);
        h1_inTofRaw[index]->Fill(raw.fFC_TofRaw[i]);
        h1_inTofRaw_GammaPeak[index]->Fill(raw.fFC_TofRaw[i]);
        if (q3>0) h2_Q2Q3vQ1[index]->Fill(q1,q2/q3);
     }


  }//end of loop over the entries 

  cout << endl;

  anode = 0 ;
  TCanvas * can_Q[m_nDets];
  TCanvas * can_Q1vT[m_nDets];
  TCanvas * can_discri[m_nDets];
  //TCanvas * can_tof_raw[m_nAnodesTot];
  //TCanvas * can_gamma_peak[m_nAnodesTot];
  for(unsigned short d = 0 ; d < m_nDets; d++){

    int ncol = ceil(0.5 * m_nAnodes[d]);
    cout << "det : " << d+1 << " ncol:" << ncol << endl;

    sprintf(name,"Q_EPIC%i",d+1);
    can_Q[d] = new TCanvas(name,name,0,0,2500,1500);
    can_Q[d]->Divide(ncol,2);   

    sprintf(name,"Q1vT_EPIC%i",d+1);
    can_Q1vT[d] = new TCanvas(name,name,0,0,2500,1500);
    can_Q1vT[d]->Divide(ncol,2);   

    sprintf(name,"DISCRI_EPIC%i",d+1);
    can_discri[d] = new TCanvas(name,name,0,0,2500,1500);
    can_discri[d]->Divide(ncol,2);   
 

    for(unsigned short a = 0 ; a < m_nAnodes[d]; a++){

       can_Q[d]->cd(a+1);
       h1_Q1[anode]->Draw();  h1_Q2[anode]->Draw("same"); h1_Q3[anode]->Draw("same");

       can_Q1vT[d]->cd(a+1);
       h2_Q1vT[anode]->Draw("colz");

       can_discri[d]->cd(a+1);
       h2_Q2Q3vQ1[anode]->Draw("colz");

       //sprintf(name,"inTofRaw_EPIC%i_ANODE%i",d+1,a+1);
       //can_tof_raw[anode] = new TCanvas(name,name,0,0,1500,1000);
       //can_tof_raw[anode]->cd();
       //h1_inTofRaw[anode]->Draw();

       //sprintf(name,"GammaPeak_EPIC%i_ANODE%i",d+1,a+1);
       //can_gamma_peak[anode] = new TCanvas(name,name,0,0,1500,1000);
       //can_gamma_peak[anode]->cd();
       //h1_inTofRaw_GammaPeak[anode]->Draw();
       anode++;
    }
  }

}
