#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>

void run()
{

	const int nA = 9;
	int anode[nA] = {1,2,3,4,6,8,9,10,11};
	char name[100];

	TFile * f4b_1030_550 = new TFile("results/Histo_V4b_1030mbar_550V_1.root","read");
	TFile * f4b_1190_500 = new TFile("results/Histo_V4b_1190mbar_500V_3.root","read");
	TFile * f4b_1190_610 = new TFile("results/Histo_V4b_1190mbar_610V_2.root","read");

	TCanvas * can_Q1 = new TCanvas("V4b_Q1","V4b_Q1",0,0,2000,1500);
	can_Q1->Divide(5,2);
	TH1F * h4b_1030_550_Q1[nA];
	TH1F * h4b_1190_500_Q1[nA];
	TH1F * h4b_1190_610_Q1[nA];

	TLatex t4b_1030_550;  t4b_1030_550.SetNDC();  t4b_1030_550.SetTextSize(0.04);  t4b_1030_550.SetTextColor(kBlue);
	TLatex t4b_1190_500;  t4b_1190_500.SetNDC();  t4b_1190_500.SetTextSize(0.04);  t4b_1190_500.SetTextColor(kGreen+2);
	TLatex t4b_1190_610;  t4b_1190_610.SetNDC();  t4b_1190_610.SetTextSize(0.04);  t4b_1190_610.SetTextColor(8);


	for(int a = 0 ; a < nA ; a++){

		// === ========================================= 
		// === Q1 

		sprintf(name,"Q1_A%i_1",anode[a]);
		h4b_1030_550_Q1[a] = (TH1F*)f4b_1030_550->Get(name);
		sprintf(name,"Q1_A%i_P1030mbar_550V",anode[a]);
		h4b_1030_550_Q1[a]->SetTitle(name);
		h4b_1030_550_Q1[a]->SetName(name);
		h4b_1030_550_Q1[a]->SetLineColor(kBlue);
		h4b_1030_550_Q1[a]->Rebin(20);
		h4b_1030_550_Q1[a]->SetDirectory(0);

		sprintf(name,"Q1_A%i_3",anode[a]);
		h4b_1190_500_Q1[a] = (TH1F*)f4b_1190_500->Get(name);
		sprintf(name,"Q1_A%i_P1040mbar_550V",anode[a]);
		h4b_1190_500_Q1[a]->SetTitle(name);
		h4b_1190_500_Q1[a]->SetName(name);
		h4b_1190_500_Q1[a]->SetLineColor(kGreen+2);
		h4b_1190_500_Q1[a]->Rebin(20);
		h4b_1190_500_Q1[a]->Scale(509667./118438.);
		h4b_1190_500_Q1[a]->SetDirectory(0);

		sprintf(name,"Q1_A%i_2",anode[a]);
		h4b_1190_610_Q1[a] = (TH1F*)f4b_1190_610->Get(name);
		sprintf(name,"Q1_A%i_P1380mbar_650V",anode[a]);
		h4b_1190_610_Q1[a]->SetTitle(name);
		h4b_1190_610_Q1[a]->SetName(name);
		h4b_1190_610_Q1[a]->SetLineColor(8);
		h4b_1190_610_Q1[a]->Rebin(20);
		h4b_1190_610_Q1[a]->Scale(509667./128430.);
		h4b_1190_610_Q1[a]->SetDirectory(0);

		can_Q1->cd(a+1);  gPad->SetLogy(); 
		h4b_1030_550_Q1[a]->Draw("hist");		
		h4b_1190_500_Q1[a]->Draw("hist same");		
		h4b_1190_610_Q1[a]->Draw("hist same");		

		t4b_1030_550.DrawLatex(0.5,0.85,"1040 mbar, 550 V");
		t4b_1190_500.DrawLatex(0.5,0.80,"1190 mbar, 500 V");
		t4b_1190_610.DrawLatex(0.5,0.75,"1190 mbar, 610 V");

	
	}


}
