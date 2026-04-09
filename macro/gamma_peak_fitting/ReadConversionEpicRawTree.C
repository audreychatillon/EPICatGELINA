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

void run(UShort_t run_number)
{

  char name[100];

  // TODO remove hard coding
  const unsigned short m_nDets = 1;
  unsigned short m_nAnodes[m_nDets] = {11};  
  const unsigned short m_nAnodesTot = 11;


  // === =========================================================
  // === histograms
  TH1F * h1_inTofRaw[m_nAnodesTot];
  TH1F * h1_inTofRaw_GammaPeak[m_nAnodesTot];
  int anode = 0;
  for(unsigned short d = 0 ; d < m_nDets; d++){
    for(unsigned short a = 0 ; a < m_nAnodes[d]; a++){
       sprintf(name,"inTofRaw_EPIC%i_A%i",d+1,a+1);
       h1_inTofRaw[anode] = new TH1F(name,name,60000,-2000,4000);
       sprintf(name,"inTofRaw_GammaPeak_EPIC%i_A%i",d+1,a+1);
       h1_inTofRaw_GammaPeak[anode] = new TH1F(name,name,10000,500,1500);
       anode++;
    }
  }

  // === =========================================================
  // === input data 
  TChain * ch = new TChain("EpicRawTree");
  sprintf(name,"../../output/conversion/raw_run%i.root",run_number);
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
        int index = 0;
        for(int d=0; d<det; d++){
	  index += d*m_nAnodes[d];
	}
        index += anode - 1;
        h1_inTofRaw[index]->Fill(raw.fFC_TofRaw[i]);
        h1_inTofRaw_GammaPeak[index]->Fill(raw.fFC_TofRaw[i]);
     }


  }//end of loop over the entries 

  cout << endl;

  anode = 0 ;
  TCanvas * can_tof_raw[m_nAnodesTot];
  TCanvas * can_gamma_peak[m_nAnodesTot];
  for(unsigned short d = 0 ; d < m_nDets; d++){
    for(unsigned short a = 0 ; a < m_nAnodes[d]; a++){
       sprintf(name,"inTofRaw_EPIC%i_ANODE%i",d+1,a+1);
       can_tof_raw[anode] = new TCanvas(name,name,0,0,1500,1000);
       can_tof_raw[anode]->cd();
       h1_inTofRaw[anode]->Draw();
       sprintf(name,"GammaPeak_EPIC%i_ANODE%i",d+1,a+1);
       can_gamma_peak[anode] = new TCanvas(name,name,0,0,1500,1000);
       can_gamma_peak[anode]->cd();
       h1_inTofRaw_GammaPeak[anode]->Draw();
       anode++;
    }
  }

}
