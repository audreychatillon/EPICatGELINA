#include <iostream>
#include <fstream>
#include <cmath>

#include "../EpicPretreatData.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include "TGraph.h"
#include <TH1.h>
#include <TH2.h>

#define FLIGHT_PATH 7320
#define MASS_U238   97.97
#define MASS_U235  622.62
#define MASS_PU242 281.88

////////////////////////////////////////////////////////////////////////////////
double TofRaw2Ene(double gammapeak, double pos, double tofraw, double &tofcal) {

  const double mn_MeV = 939.565;
  double offset = pos / 299.792458 - gammapeak;
  tofcal = tofraw + offset;
  double beta = (pos / tofcal) / 299.792458;
  double gamma = 1.0 / sqrt(1.0 - pow(beta, 2));
  return ((gamma - 1.0) * mn_MeV);
}


////////////////////////////////////////////////////////////////////////////////
void run()
{

  // ==========================================================================================
  // === VARIABLES
  
  long long int nentries = 0;
  double tofraw = -1.;

  // ==========================================================================================
  // === HISTOGRAMS
  TH1D * h1_tofraw[10];
  for(short a = 0 ; a < 10 ; a++){
	  h1_tofraw[a] = new TH1D(Form("tofraw_A%02d",a+1),Form("tofraw_A%02d",a+1),25050000,-5000,2500000);
    h1_tofraw[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
    h1_tofraw[a]->SetLineColor(kRed+1);
    h1_tofraw[a]->SetDirectory(0);
  }

  // ===========================================================================================
  // === INPUT DATA
  TChain *ch = new TChain("EpicPretreatData");
  for(short runid = 17; runid <=31; runid++){
    if(runid==20) continue;
    ch->Add(Form("../../../output/pretreat/pretreat_run%i.root",runid)); 
  }
  EpicPretreatData epic(ch);
  nentries = ch->GetEntries();
  cout << "nentries = " << nentries << endl;

  // ===========================================================================================
  // === LOOP OVER THE DATA
  for(long long int entry = 0 ; entry < nentries ; entry++){
	  epic.GetEntry(entry);
    if ((entry % 5000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;

	  if(epic.fFC_fission == true){
	    tofraw = epic.fFC_time - epic.fHF_time_next;
      if(tofraw < -1000) tofraw = epic.fFC_time-epic.fHF_time_curr;
	    h1_tofraw[epic.fFC_anode-1]->Fill(tofraw);
    }

  }// end of loop over the entries
  cout << endl;


  // ===========================================================================================
  // === DRAW
  
  TCanvas * can1 = new TCanvas("TofRaw","TofRaw",0,0,3000,2000);
  can1->Divide(5,2);
  for(short a = 1 ; a <=10 ; a++){
      can1->cd(a); gPad->SetLogy();
      h1_tofraw[a-1]->Draw();
  }
}
