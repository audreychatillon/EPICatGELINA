#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>


void run()
{

	const int nA = 4;
	int anode[nA] = {3,4,9,10};
	char name[100];

	TFile * f82_1380_1 = new TFile("results/Histo9_V8-2_1380mbar_650V.C","read");
	TFile * f82_1040_1 = new TFile("results/Histo10_V8-2_1040mbar_550V.C","read");

	TCanvas * can_Q1 = new TCanvas("Q1","Q1",0,0,2000,1500);
	can_Q1->Divide(2,2);
	TH1F * h82_1380_Q1_1[nA];
	TH1F * h82_1040_Q1_1[nA];

	for(int a = 0 ; a < nA ; a++){

		sprintf(name,"Q1_EPIC1_A%i",anode[a]);
		h82_1380_Q1_1[a] = (TH1F*)f82_1380_1->Get(name);
		sprintf(name,"Q1_EPIC1_A%i_P1380mbar_650V_1",anode[a]);
		h82_1380_Q1_1[a]->SetTitle(name);
		h82_1380_Q1_1[a]->SetName(name);
		h82_1380_Q1_1[a]->SetLineColor(kRed);
		h82_1380_Q1_1[a]->SetDirectory(0);
		can_Q1->cd(a+1);
		h82_1380_Q1_1[a]->Draw();		

		sprintf(name,"Q1_EPIC1_A%i",anode[a]);
		h82_1040_Q1_1[a] = (TH1F*)f82_1040_1->Get(name);
		sprintf(name,"Q1_EPIC1_A%i_P1040mbar_550V_1",anode[a]);
		h82_1040_Q1_1[a]->SetTitle(name);
		h82_1040_Q1_1[a]->SetName(name);
		h82_1040_Q1_1[a]->SetLineColor(kBlue);
		h82_1040_Q1_1[a]->SetDirectory(0);
		can_Q1->cd(a+1);
		h82_1040_Q1_1[a]->Draw("same");		
		
	}


}
