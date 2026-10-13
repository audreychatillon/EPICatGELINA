#include <iostream>

#include "../EpicPretreatData.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>

void run()
{


    // ==========================================================================================
    // === VARIABLES

    long long int nentries = 0;

    // ==========================================================================================
    // === HISTOGRAMS

    TH1D * h1_tofraw_curr[10];
    TH1D * h1_tofraw_next[10];
    for(short a = 0 ; a < 10 ; a++){
	h1_tofraw_curr[a] = new TH1D(Form("tofraw_curr_A%02d",a+1),Form("tofraw_curr_A%02d",a+1),30000,2498500,2501500);
        h1_tofraw_curr[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_curr[a]->SetLineColor(kBlue+1);
        h1_tofraw_curr[a]->SetDirectory(0);

	h1_tofraw_next[a] = new TH1D(Form("tofraw_next_A%02d",a+1),Form("tofraw_next_A%02d",a+1),30000,-1500,1500);
        h1_tofraw_next[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_next[a]->SetLineColor(kBlack+1);
        h1_tofraw_next[a]->SetDirectory(0);
    }


    // ===========================================================================================
    // === INPUT DATA
    TChain *ch = new TChain("EpicPretreatData");
    ch->Add(Form("../../../output/pretreat/pretreat_run18.root")); 
    EpicPretreatData epic(ch);
    nentries = ch->GetEntries();
    cout << "nentries = " << nentries << endl;

    // ===========================================================================================
    // === LOOP OVER THE DATA
    for(long long int entry = 0 ; entry < nentries ; entry++){
	epic.GetEntry(entry);
    	if ((entry % 1000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;

    }// end of loop

    // ===========================================================================================
    // === DRAW
  
    TCanvas * can1 = new TCanvas("TofRawCurr","TofRawCurr",0,0,3000,2000);
    can1->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can1->cd(a); gPad->SetLogy();
        h1_tofraw_curr[a-1]->Draw();
    }

    TCanvas * can2 = new TCanvas("TofRawNext","TofRawNext",0,0,3000,2000);
    can2->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can2->cd(a); gPad->SetLogy();
        h1_tofraw_next[a-1]->Draw();
    }
}
