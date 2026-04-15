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

void run(UShort_t run_number, string PA, int Pmbar, int HV)
{

  char name[100];

  // TODO remove hard coding
  const unsigned short m_nDets = 1;
  unsigned short m_nAnodes[m_nDets] = {11};  
  const unsigned short m_nAnodesTot = 11;


  // === =========================================================
  // === histograms
  TH1F * h1_TimeHF = new TH1F("TimeHF","TimeHF",86400,0,86400);;
  TH1F * h1_DeltaTimeHF = new TH1F("DeltaTimeHF","DeltaTimeHF",10000,2,3);;
  TH1F * h1_Mult[m_nDets];
  TH2F * h2_MvT[m_nDets];
  TH1F * h1_Q1[m_nAnodesTot];
  TH2F * h2_Q1vT[m_nAnodesTot];
  TH1F * h1_Q2[m_nAnodesTot];
  TH1F * h1_Q3[m_nAnodesTot];
  TH2F * h2_Q2Q3vQ1[m_nAnodesTot];
  //TH1F * h1_inTofRaw[m_nAnodesTot];
  //TH1F * h1_inTofRaw_GammaPeak[m_nAnodesTot];

  int anode = 0;
  for(unsigned short d = 0 ; d < m_nDets; d++){

    sprintf(name,"EPIC%i_Mult",d+1);
    h1_Mult[d] = new TH1F(name,name,13,-0.5,12.5);
    h1_Mult[d]->SetDirectory(0);

    sprintf(name,"EPIC%i_Mult_vs_T",d+1);
    h2_MvT[d] = new TH2F(name,name,2400,0,72000,13,-0.5,12.5);
    h2_MvT[d]->SetDirectory(0);

    for(unsigned short a = 0 ; a < m_nAnodes[d]; a++){
       sprintf(name,"Q1vT_EPIC%i_A%i",d+1,a+1);
       h2_Q1vT[anode] = new TH2F(name,name,2400,0,72000,5000,0,500000);
       h2_Q1vT[anode]->GetXaxis()->SetTitle("Time [s] 30s/bin");
       h2_Q1vT[anode]->SetDirectory(0);

       sprintf(name,"Q1_EPIC%i_A%i",d+1,a+1);
       h1_Q1[anode] = new TH1F(name,name,50000,0,500000);
       h1_Q1[anode]->SetDirectory(0);

       sprintf(name,"Q2_EPIC%i_A%i",d+1,a+1);
       h1_Q2[anode] = new TH1F(name,name,50000,0,500000);
       h1_Q2[anode]->SetLineColor(8);
       h1_Q2[anode]->SetDirectory(0);

       sprintf(name,"Q3_EPIC%i_A%i",d+1,a+1);
       h1_Q3[anode] = new TH1F(name,name,50000,0,500000);
       h1_Q3[anode]->SetLineColor(kCyan);
       h1_Q3[anode]->SetDirectory(0);

       sprintf(name,"discri_EPIC%i_A%i",d+1,a+1);
       h2_Q2Q3vQ1[anode] = new TH2F(name,name,2500,0,500000,300,0,3);
       h2_Q2Q3vQ1[anode]->SetDirectory(0);

       //sprintf(name,"inTofRaw_EPIC%i_A%i",d+1,a+1);
       //h1_inTofRaw[anode] = new TH1F(name,name,60000,-20000,40000);

       //sprintf(name,"inTofRaw_GammaPeak_EPIC%i_A%i",d+1,a+1);
       //h1_inTofRaw_GammaPeak[anode] = new TH1F(name,name,10000,500,1500);

       anode++;
    }
  }

  // === =========================================================
  // === input data 
  TChain * ch = new TChain("EpicRawTree");
  sprintf(name,"../../output/conversion/Test%i_V%s_%imbar_%iV.root",run_number,PA.c_str(),Pmbar,HV);
  ch->Add(name);
  ch->ls();
  EpicRawTree raw(ch);
  ULong64_t nentries = (ULong64_t)ch->GetEntries();
  cout << "number of entries: " << nentries << endl;
  for(ULong64_t Entry=0; Entry<nentries; Entry++){
    raw.GetEntry(Entry);
    if ((Entry % 1000000)==0) cout << "\r === Entry = " << Entry << " === " << flush;
    //if ((Entry % 2)==0) cout << "\r === Entry = " << Entry << " === " << flush;
 
    double Qmax[m_nDets];
    int    Imax[m_nDets];
    int    Mult[m_nDets];
    for(unsigned short d = 0 ; d < m_nDets; d++){
        Qmax[d] = 0.;
        Imax[d] = -1;
    	Mult[d] = 0;
    }

    // get the channel with Qmax
    int mult = raw.fFC_DetNbr.size();
    if (mult>0){
      if(raw.fFC_DetNbr[0] < 0){
         h1_TimeHF->Fill(1.e-09*raw.fHF_Time); 
         h1_DeltaTimeHF->Fill(1.e-06*(raw.fHF_Time-raw.fHF_TimePrev)); 
      } // end of if HF data
      else{
        for(int i=0; i<mult; i++){
            int det = raw.fFC_DetNbr[i];
            double qm = raw.fFC_Qmax[i];
            Mult[det-1]++;
            if(Qmax[det-1]<qm && det>=0 ){
            	Qmax[det-1] = qm;
            	Imax[det-1] = i;
            }
        }

        // read the channel with Qmax
        for(int d=0; d<m_nDets; d++){
            int det = raw.fFC_DetNbr[Imax[d]];
            int anode = raw.fFC_AnodeNbr[Imax[d]];
            double q1 = raw.fFC_Q1[Imax[d]];
            double q2 = raw.fFC_Q2[Imax[d]];
            double q3 = raw.fFC_Q3[Imax[d]];
            double t_s = raw.fFC_Time[Imax[d]] * 1.e-09 ;
            int index = 0;
            for(int i=0; i<det; i++){
              index += i*m_nAnodes[i];
            }
            index += anode - 1;
            h1_Mult[d]->Fill(Mult[d]);
            h2_MvT[d]->Fill(t_s,Mult[d]);
            h1_Q1[index]->Fill(q1);
            h2_Q1vT[index]->Fill(t_s,q1);
            h1_Q2[index]->Fill(q2);
            h1_Q3[index]->Fill(q3);
            //h1_inTofRaw[index]->Fill(raw.fFC_TofRaw[i]);
            //h1_inTofRaw_GammaPeak[index]->Fill(raw.fFC_TofRaw[i]);
            if (q3>0) h2_Q2Q3vQ1[index]->Fill(q1,q2/q3);
        }
      }// end of if FC data
    }// end of if mult > 0
  }//end of loop over the entries 

  cout << endl;

  anode = 0 ;
  TCanvas * can_T0 = new TCanvas("T0","T0",0,0,2500,1500);
  TCanvas * can_mult = new TCanvas("Mult","Mult",0,0,2500,1500);
  TCanvas * can_Q[m_nDets];
  TCanvas * can_Q1vT[m_nDets];
  TCanvas * can_discri[m_nDets];
  //TCanvas * can_tof_raw[m_nAnodesTot];
  //TCanvas * can_gamma_peak[m_nAnodesTot];
  
  can_T0->Divide(1,2);
  can_T0->cd(1); h1_TimeHF->Draw();
  can_T0->cd(2); h1_DeltaTimeHF->Draw(); 
  unsigned int effective_beam_time_s = 0;
  for (int b = 1 ; b <= h1_TimeHF->GetNbinsX(); b++){
	if(h1_TimeHF->GetBinContent(b) != 0) effective_beam_time_s++;
  }
  cout << "EFFECTIVE BEAM TIME = " << effective_beam_time_s << " s" << endl;
  cout << "NUMBER OF HF: = " << h1_DeltaTimeHF->Integral() << " " << endl;


  
  sprintf(name,"results/Histo%i_V%s_%imbar_%iV.C",run_number,PA.c_str(),Pmbar,HV);
  TFile * fsave = new TFile(name,"recreate");
  fsave->cd();
  h1_TimeHF->Write();
  h1_DeltaTimeHF->Write();

  can_mult->Divide(2,m_nDets);
  for(unsigned short d = 0 ; d < m_nDets; d++){
    can_mult->cd(2*d+1); h1_Mult[d]->Draw();
    can_mult->cd(2*d+2); h2_MvT[d]->Draw("colz");
    fsave->cd(); 
    h1_Mult[d]->Write();
    h2_MvT[d]->Write();

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

       can_discri[d]->cd(a+1); gPad->SetLogx();
       h2_Q2Q3vQ1[anode]->Draw("colz");

       //sprintf(name,"inTofRaw_EPIC%i_ANODE%i",d+1,a+1);
       //can_tof_raw[anode] = new TCanvas(name,name,0,0,1500,1000);
       //can_tof_raw[anode]->cd();
       //h1_inTofRaw[anode]->Draw();

       //sprintf(name,"GammaPeak_EPIC%i_ANODE%i",d+1,a+1);
       //can_gamma_peak[anode] = new TCanvas(name,name,0,0,1500,1000);
       //can_gamma_peak[anode]->cd();
       //h1_inTofRaw_GammaPeak[anode]->Draw();
       //
       fsave->cd(); 
       h1_Q1[anode]->Write();  
       h1_Q2[anode]->Write(); 
       h1_Q3[anode]->Write();
       h2_Q1vT[anode]->Write();
       h2_Q2Q3vQ1[anode]->Write();

       anode++;
    }
  }

  fsave->Close();


}
