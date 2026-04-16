#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>

void run()
{

	const int nA = 4;
	int anode[nA] = {3,4,9,10};
	char name[100];

	TFile * f82_1040 = new TFile("results/Histo_V8-2_1040mbar_550V_5.root","read");
	TFile * f82_1075 = new TFile("results/Histo_V8-2_1075mbar_720V_8.root","read");
	TFile * f82_1210 = new TFile("results/Histo_V8-2_1210mbar_680V_6.root","read");
	TFile * f82_1370 = new TFile("results/Histo_V8-2_1370mbar_720V_7.root","read");
	TFile * f82_1380 = new TFile("results/Histo_V8-2_1380mbar_650V_4.root","read");

	TCanvas * can_Q1 = new TCanvas("V82_Q1","B82_Q1",0,0,2000,1500);
	can_Q1->Divide(2,2);

	TH1F * h82_1040_Q1[nA];
	TH1F * h82_1075_Q1[nA];
	TH1F * h82_1210_Q1[nA];
	TH1F * h82_1370_Q1[nA];
	TH1F * h82_1380_Q1[nA];

	TCanvas * can_discri[nA]; 

	TH1F * h82_1040_discri[nA];
	TH1F * h82_1075_discri[nA];
	TH1F * h82_1210_discri[nA];
	TH1F * h82_1370_discri[nA];
	TH1F * h82_1380_discri[nA];

	TLatex t82_1040;  t82_1040.SetNDC();  t82_1040.SetTextSize(0.025);  t82_1040.SetTextColor(kBlue);
	TLatex t82_1075;  t82_1075.SetNDC();  t82_1075.SetTextSize(0.025);  t82_1075.SetTextColor(kCyan+2);
	TLatex t82_1210;  t82_1210.SetNDC();  t82_1210.SetTextSize(0.025);  t82_1210.SetTextColor(8);
	TLatex t82_1370;  t82_1370.SetNDC();  t82_1370.SetTextSize(0.025);  t82_1370.SetTextColor(kRed+1);
	TLatex t82_1380;  t82_1380.SetNDC();  t82_1380.SetTextSize(0.025);  t82_1380.SetTextColor(kMagenta);

	for(int a = 0 ; a < nA ; a++){

		// === ========================================= 
		// === Q1 

		sprintf(name,"Q1_A%i_6",anode[a]);
		h82_1210_Q1[a] = (TH1F*)f82_1210->Get(name);
		sprintf(name,"Q1_A%i_P1210mbar_680V",anode[a]);
		h82_1210_Q1[a]->SetTitle(name);
		h82_1210_Q1[a]->SetName(name);
		h82_1210_Q1[a]->SetLineColor(8);
		h82_1210_Q1[a]->Rebin(20);
		h82_1210_Q1[a]->SetDirectory(0);
		

		sprintf(name,"Q1_A%i_5",anode[a]);
		h82_1040_Q1[a] = (TH1F*)f82_1040->Get(name);
		sprintf(name,"Q1_A%i_P1040mbar_550V",anode[a]);
		h82_1040_Q1[a]->SetTitle(name);
		h82_1040_Q1[a]->SetName(name);
		h82_1040_Q1[a]->SetLineColor(kBlue);
		h82_1040_Q1[a]->Rebin(20);
		h82_1040_Q1[a]->Scale(42866./19333.);
		h82_1040_Q1[a]->SetDirectory(0);

		sprintf(name,"Q1_A%i_4",anode[a]);
		h82_1380_Q1[a] = (TH1F*)f82_1380->Get(name);
		sprintf(name,"Q1_A%i_P1380mbar_650V",anode[a]);
		h82_1380_Q1[a]->SetTitle(name);
		h82_1380_Q1[a]->SetName(name);
		h82_1380_Q1[a]->SetLineColor(kMagenta);
		h82_1380_Q1[a]->Rebin(20);
		h82_1380_Q1[a]->Scale(42866./27868.);
		h82_1380_Q1[a]->SetDirectory(0);

		sprintf(name,"Q1_A%i_7",anode[a]);
		h82_1370_Q1[a] = (TH1F*)f82_1370->Get(name);
		sprintf(name,"Q1_A%i_P1370mbar_720V",anode[a]);
		h82_1370_Q1[a]->SetTitle(name);
		h82_1370_Q1[a]->SetName(name);
		h82_1370_Q1[a]->SetLineColor(kRed+1);
		h82_1370_Q1[a]->Rebin(20);
		h82_1370_Q1[a]->Scale(42866./16939.);
		h82_1370_Q1[a]->SetDirectory(0);

		sprintf(name,"Q1_A%i_8",anode[a]);
		h82_1075_Q1[a] = (TH1F*)f82_1075->Get(name);
		sprintf(name,"Q1_A%i_P1075mbar_720V",anode[a]);
		h82_1075_Q1[a]->SetTitle(name);
		h82_1075_Q1[a]->SetName(name);
		h82_1075_Q1[a]->SetLineColor(kCyan+2);
		h82_1075_Q1[a]->SetLineWidth(2);
		h82_1075_Q1[a]->Rebin(40);
		h82_1075_Q1[a]->Scale(0.5*42866./(9429.));
		h82_1075_Q1[a]->SetDirectory(0);

		can_Q1->cd(a+1);  gPad->SetLogy(); 
		h82_1210_Q1[a]->Draw("hist");		
		h82_1040_Q1[a]->Draw("hist same");		
		h82_1380_Q1[a]->Draw("hist same");		
		h82_1370_Q1[a]->Draw("hist same");		
		h82_1075_Q1[a]->Draw("hist same");		
			
		t82_1040.DrawLatex(0.5,0.85,"1040 mbar, 550 V");
		t82_1210.DrawLatex(0.5,0.80,"1210 mbar, 680 V");
		t82_1370.DrawLatex(0.5,0.75,"1370 mbar, 720 V");
		t82_1380.DrawLatex(0.5,0.70,"1380 mbar, 650 V");
		t82_1075.DrawLatex(0.5,0.65,"???? mbar, ??? V");


		// === ========================================= 
		// === DISCRI 

		sprintf(name,"discri_A%i_6",anode[a]);
		h82_1210_discri[a] = (TH1F*)f82_1210->Get(name);
		sprintf(name,"discri_A%i_P1210mbar_680V",anode[a]);
		h82_1210_discri[a]->SetTitle(name);
		h82_1210_discri[a]->SetName(name);
		h82_1210_discri[a]->SetDirectory(0);
		h82_1210_discri[a]->GetXaxis()->SetRangeUser(7000,500000);
		h82_1210_discri[a]->GetYaxis()->SetRangeUser(0,3.5);

		sprintf(name,"discri_A%i_5",anode[a]);
		h82_1040_discri[a] = (TH1F*)f82_1040->Get(name);
		sprintf(name,"discri_A%i_P1040mbar_550V",anode[a]);
		h82_1040_discri[a]->SetTitle(name);
		h82_1040_discri[a]->SetName(name);
		h82_1040_discri[a]->SetLineColor(kBlue);
		h82_1040_discri[a]->SetDirectory(0);
		h82_1040_discri[a]->GetXaxis()->SetRangeUser(7000,500000);
		h82_1040_discri[a]->GetYaxis()->SetRangeUser(0,3.5);

		sprintf(name,"discri_A%i_4",anode[a]);
		h82_1380_discri[a] = (TH1F*)f82_1380->Get(name);
		sprintf(name,"discri_A%i_P1380mbar_650V",anode[a]);
		h82_1380_discri[a]->SetTitle(name);
		h82_1380_discri[a]->SetName(name);
		h82_1380_discri[a]->SetLineColor(kMagenta);
		h82_1380_discri[a]->SetDirectory(0);
		h82_1380_discri[a]->GetXaxis()->SetRangeUser(7000,500000);
		h82_1380_discri[a]->GetYaxis()->SetRangeUser(0,3.5);

		sprintf(name,"discri_A%i_7",anode[a]);
		h82_1370_discri[a] = (TH1F*)f82_1370->Get(name);
		sprintf(name,"discri_A%i_P1370mbar_720V",anode[a]);
		h82_1370_discri[a]->SetTitle(name);
		h82_1370_discri[a]->SetName(name);
		h82_1370_discri[a]->SetLineColor(kRed+1);
		h82_1370_discri[a]->SetDirectory(0);
		h82_1370_discri[a]->GetXaxis()->SetRangeUser(7000,500000);
		h82_1370_discri[a]->GetYaxis()->SetRangeUser(0,3.5);

		sprintf(name,"discri_A%i_8",anode[a]);
		h82_1075_discri[a] = (TH1F*)f82_1075->Get(name);
		sprintf(name,"discri_A%i_P1075mbar_720V",anode[a]);
		h82_1075_discri[a]->SetTitle(name);
		h82_1075_discri[a]->SetName(name);
		h82_1075_discri[a]->SetLineColor(kCyan+2);
		h82_1075_discri[a]->SetDirectory(0);
		h82_1075_discri[a]->GetXaxis()->SetRangeUser(7000,500000);
		h82_1075_discri[a]->GetYaxis()->SetRangeUser(0,3.5);
	
		sprintf(name,"discri_A%i",anode[a]);
		can_discri[a] = new TCanvas(name,name,0,0,1500,1000);
		can_discri[a]->Divide(2,3);
		can_discri[a]->cd(1); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogz(); gPad->SetLogx(); h82_1040_discri[a]->Draw("colz"); 
		can_discri[a]->cd(2); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogz(); gPad->SetLogx(); h82_1075_discri[a]->Draw("colz");
		can_discri[a]->cd(3); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogz(); gPad->SetLogx(); h82_1210_discri[a]->Draw("colz");
		can_discri[a]->cd(4); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogz(); gPad->SetLogx(); h82_1370_discri[a]->Draw("colz");
		can_discri[a]->cd(5); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogz(); gPad->SetLogx(); h82_1380_discri[a]->Draw("colz");

	
	}

	TCanvas * canSANDIA = new TCanvas("SANDIA","SANDIA",0,0,1000,1000);
	canSANDIA->cd(); gPad->SetLogy(); 
	h82_1040_Q1[3]->Draw("hist");        
	h82_1040_Q1[3]->GetXaxis()->SetRangeUser(15000,400000);        
	h82_1040_Q1[3]->GetXaxis()->SetTitle("total charge");        
	h82_1210_Q1[3]->Draw("hist same");
	h82_1370_Q1[3]->Draw("hist same");
	h82_1380_Q1[3]->Draw("hist same");
	h82_1075_Q1[3]->Draw("hist same");
	
	t82_1040.DrawLatex(0.15,0.91,"1040 mbar, 550 V");
	t82_1210.DrawLatex(0.15,0.87,"1210 mbar, 680 V");
	t82_1370.DrawLatex(0.15,0.83,"1370 mbar, 720 V");
	t82_1380.DrawLatex(0.15,0.79,"1380 mbar, 650 V");
	t82_1075.DrawLatex(0.15,0.75,"1075 mbar, 720 V");

	TPad * p1 = new TPad("discri","discri",0.5,0.6,1.0,1.0);
	canSANDIA->cd(); p1->Draw(); p1->cd(); 
	gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogz(); gPad->SetLogx(); 
	h82_1075_discri[3]->Draw("col");
	h82_1075_discri[3]->GetXaxis()->SetTitle("total charge");
	h82_1075_discri[3]->GetYaxis()->SetTitle("pulse shape discri");

	TH1F * h82clone_1040_Q1 = (TH1F*)h82_1040_Q1[3]->Clone();
	TH1F * h82clone_1075_Q1 = (TH1F*)h82_1075_Q1[3]->Clone();
	TH1F * h82clone_1210_Q1 = (TH1F*)h82_1210_Q1[3]->Clone();
	TH1F * h82clone_1370_Q1 = (TH1F*)h82_1370_Q1[3]->Clone();
	TH1F * h82clone_1380_Q1 = (TH1F*)h82_1380_Q1[3]->Clone();
	TPad * p2 = new TPad("discri","discri",0.15,0.06,0.4,0.4);
	canSANDIA->cd(); p2->Draw(); p2->cd(); 
	h82clone_1040_Q1->Draw("hist");        
	h82clone_1040_Q1->GetXaxis()->SetRangeUser(15000,40000);        
	h82clone_1040_Q1->GetXaxis()->SetTitle("total charge");        
	h82clone_1210_Q1->Draw("hist same");
	h82clone_1370_Q1->Draw("hist same");
	h82clone_1380_Q1->Draw("hist same");
	h82clone_1075_Q1->Draw("hist same");
}
