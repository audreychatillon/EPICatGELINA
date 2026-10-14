#include <iostream>
#include <fstream>
#include <map>

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


void run()
{


    // ==========================================================================================
    // === VARIABLES

    long long int nentries = 0;

    double gamma_peak_curr[10] = {
	-1746.7     , -1.74676e+03, -1.74700e+03, -1.74741e+03, -1.74305e+03, 
	-1.74109e+03, -1.73917e+03, -1.73850e+03, -1.74792e+03, -1.74648e+03 
    };

    double gamma_peak_next[10] = {
	-887.085, -887.197, -887.415, -887.778, -883.432, 
	-881.499, -879.553, -878.917, -888.23,  -886.88  
    };
    double flight_path[10] = {33.8, 27.2, 20.2, 13.6, 6.6, 0., -7.0, -13.6, -20.2, -27.2};

    double actinide_mass[10] = {97.97,  99.23,  103.8,  103.31,  92.67,  95.14,  94.07,  112.69,  107.55,  105.04 };

    double tof_raw_curr = -1.;
    double tof_raw_next = -1.;
    double tof_cal_curr = -1.;
    double tof_cal_next = -1.;
    double e_curr = -1.;
    double e_next = -1.;

    // === =========================================================
    // === TCutG for alpha and fission
    TFile * fcutg = new TFile("../2DdiscriTCutG.root","read");
    TCutG * cutalpha[10];
    TCutG * cutfission[10];
    for(short a = 0 ; a < 10 ; a++){
	cutalpha[a] = (TCutG*)fcutg->Get(Form("A%02d_alpha",a+1));
	cutfission[a] = (TCutG*)fcutg->Get(Form("A%02d_fission",a+1));
    } 
    fcutg->Close();

    // === =========================================================
    std::map<std::string,TGraph*> m_phy_gr;
    // JEFF4-0 cross sections
    vector<double> XS_JEFF_nENE;
    vector<double> XS_JEFF_ratio_235U;
    vector<double> XS_JEFF_ratio_238U;
    vector<double> XS_JEFF_ratio_242Pu;
    ifstream XSfile("../../../doc/XS_JEFF4-0_235U_238U_242Pu.dat");	
    if (!XSfile.is_open()) {
        cerr << "Erreur : impossible d'ouvrir le fichier doc/XS_JEFF4-0_235U_238U_242Pu.dat" << endl;
    }
    else{
        double e, u8, u5, pu42;  //e is in eV
        while(XSfile >> e >> u5 >> u8 >> pu42){
	        if (e<=20000000){
    	           XS_JEFF_nENE.push_back(e*1.e-06);
	           XS_JEFF_ratio_235U.push_back(u5/u5);
	           XS_JEFF_ratio_238U.push_back(u8/u5);
	           XS_JEFF_ratio_242Pu.push_back(pu42/u5);
	        }
	    }
        XSfile.close();
        m_phy_gr["XSratio_JEFF_235U_235U"]  = new TGraph(XS_JEFF_nENE.size(),XS_JEFF_nENE.data(),XS_JEFF_ratio_235U.data());   	
        m_phy_gr["XSratio_JEFF_235U_235U"]->SetName("XSratio_JEFF40_235U_235U");
        m_phy_gr["XSratio_JEFF_235U_235U"]->SetTitle("XSratio_JEFF40_235U_235U");
        m_phy_gr["XSratio_JEFF_235U_235U"]->SetLineColor(8);
        m_phy_gr["XSratio_JEFF_235U_235U"]->SetLineWidth(2);
        m_phy_gr["XSratio_JEFF_238U_235U"]  = new TGraph(XS_JEFF_nENE.size(),XS_JEFF_nENE.data(),XS_JEFF_ratio_238U.data());   	
        m_phy_gr["XSratio_JEFF_238U_235U"]->SetName("XSratio_JEFF40_238U_235U");
        m_phy_gr["XSratio_JEFF_238U_235U"]->SetTitle("XSratio_JEFF40_238U_235U");
        m_phy_gr["XSratio_JEFF_238U_235U"]->SetLineColor(8);
        m_phy_gr["XSratio_JEFF_238U_235U"]->SetLineWidth(2);
        m_phy_gr["XSratio_JEFF_242Pu_235U"] = new TGraph(XS_JEFF_nENE.size(),XS_JEFF_nENE.data(),XS_JEFF_ratio_242Pu.data());   	
        m_phy_gr["XSratio_JEFF_242Pu_235U"]->SetName("XSratio_JEFF40_242Pu_235U");
        m_phy_gr["XSratio_JEFF_242Pu_235U"]->SetTitle("XSratio_JEFF40_242Pu_235U");
        m_phy_gr["XSratio_JEFF_242Pu_235U"]->SetLineColor(8);
        m_phy_gr["XSratio_JEFF_242Pu_235U"]->SetLineWidth(2);
    }

    // BVIII cross sections
    vector<double> XS_BVIII_nENE;
    vector<double> XS_BVIII_ratio_235U;
    vector<double> XS_BVIII_ratio_238U;
    vector<double> XS_BVIII_ratio_242Pu;
    ifstream XSfile2("../../../doc/XS_BVIII_235U_238U_242Pu.dat");	
    if (!XSfile2.is_open()) {
        cerr << "Erreur : impossible d'ouvrir le fichier doc/XS_BVIII_235U_238U_242Pu.dat" << endl;
    }
    else{
        double e, u8, u5, pu42;  //e is in eV
        while(XSfile2 >> e >> u5 >> u8 >> pu42){
	        if (e<=20000000){
    	            XS_BVIII_nENE.push_back(e*1.e-06);
	            XS_BVIII_ratio_235U.push_back(u5/u5);
	            XS_BVIII_ratio_238U.push_back(u8/u5);
	            XS_BVIII_ratio_242Pu.push_back(pu42/u5);
	        }
	    }
        XSfile2.close();
        m_phy_gr["XSratio_BVIII_235U_235U"]  = new TGraph(XS_BVIII_nENE.size(),XS_BVIII_nENE.data(),XS_BVIII_ratio_235U.data());   	
        m_phy_gr["XSratio_BVIII_235U_235U"]->SetName("XSratio_BVIII40_235U_235U");
        m_phy_gr["XSratio_BVIII_235U_235U"]->SetTitle("XSratio_BVIII40_235U_235U");
        m_phy_gr["XSratio_BVIII_238U_235U"]  = new TGraph(XS_BVIII_nENE.size(),XS_BVIII_nENE.data(),XS_BVIII_ratio_238U.data());   	
        m_phy_gr["XSratio_BVIII_238U_235U"]->SetName("XSratio_BVIII40_238U_235U");
        m_phy_gr["XSratio_BVIII_238U_235U"]->SetTitle("XSratio_BVIII40_238U_235U");
        m_phy_gr["XSratio_BVIII_238U_235U"]->SetLineStyle(2);
        m_phy_gr["XSratio_BVIII_238U_235U"]->SetLineWidth(2);
        m_phy_gr["XSratio_BVIII_242Pu_235U"] = new TGraph(XS_BVIII_nENE.size(),XS_BVIII_nENE.data(),XS_BVIII_ratio_242Pu.data());   	
        m_phy_gr["XSratio_BVIII_242Pu_235U"]->SetName("XSratio_BVIII40_242Pu_235U");
        m_phy_gr["XSratio_BVIII_242Pu_235U"]->SetTitle("XSratio_BVIII40_242Pu_235U");
        m_phy_gr["XSratio_BVIII_242Pu_235U"]->SetLineStyle(2);
        m_phy_gr["XSratio_BVIII_242Pu_235U"]->SetLineWidth(2);
    }


    // ==========================================================================================
    // === HISTOGRAMS
    TH1D * h1_tofraw_curr[10];
    TH1D * h1_tofraw_next[10];
    TH1D * h1_tofraw_next_extraWindow[10];
    TH1D * h1_tofcal_curr[10];
    TH1D * h1_tofcal_next[10];
    TH1D * h1_e_curr[10];
    TH1D * h1_e_next[10];
    TH1D * h1_ew_curr[10];
    TH1D * h1_ew_next[10];
    TH2D * h2_discri[10];

    for(short a = 0 ; a < 10 ; a++){
	h1_tofraw_curr[a] = new TH1D(Form("tofraw_curr_A%02d",a+1),Form("tofraw_curr_A%02d",a+1),30000,-3000,0);
        h1_tofraw_curr[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_curr[a]->SetLineColor(kBlue+1);
        h1_tofraw_curr[a]->SetDirectory(0);

	h1_tofraw_next[a] = new TH1D(Form("tofraw_next_A%02d",a+1),Form("tofraw_next_A%02d",a+1),30000,-3000,0);
        h1_tofraw_next[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_next[a]->SetLineColor(kRed+1);
        h1_tofraw_next[a]->SetDirectory(0);

	h1_tofraw_next_extraWindow[a] = new TH1D(Form("tofraw_next_A%02d_extraWindow",a+1),Form("tofraw_next_A%02d_extraWindow",a+1),30000,-1500,1500);
        h1_tofraw_next_extraWindow[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_next_extraWindow[a]->SetLineColor(8);
        h1_tofraw_next_extraWindow[a]->SetDirectory(0);

	h1_tofcal_curr[a] = new TH1D(Form("tofcal_curr_A%02d",a+1),Form("tofcal_curr_A%02d",a+1),30000,-500,2500);
        h1_tofcal_curr[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofcal_curr[a]->SetLineColor(kBlue+1);
        h1_tofcal_curr[a]->SetDirectory(0);

	h1_tofcal_next[a] = new TH1D(Form("tofcal_next_A%02d",a+1),Form("tofcal_next_A%02d",a+1),30000,-500,2500);
        h1_tofcal_next[a]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofcal_next[a]->SetLineColor(kRed+1);
        h1_tofcal_next[a]->SetDirectory(0);

	h1_e_curr[a] = new TH1D(Form("En_curr_A%02d",a+1),Form("En_curr_A%02d",a+1),200,0,20);
        h1_e_curr[a]->GetXaxis()->SetTitle("[MeV] 100keV/bin");
        h1_e_curr[a]->SetLineColor(kBlue+1);
        h1_e_curr[a]->SetDirectory(0);

	h1_e_next[a] = new TH1D(Form("En_next_A%02d",a+1),Form("En_next_A%02d",a+1),200,0,20);
        h1_e_next[a]->GetXaxis()->SetTitle("[MeV] 100keV/bin");
        h1_e_next[a]->SetLineColor(kRed+1);
        h1_e_next[a]->SetDirectory(0);

	h1_ew_curr[a] = new TH1D(Form("En_weighted_curr_A%02d",a+1),Form("En_weighted_curr_A%02d",a+1),200,0,20);
        h1_ew_curr[a]->GetXaxis()->SetTitle("[MeV] 100keV/bin");
        h1_ew_curr[a]->SetLineColor(kBlue+1);
        h1_ew_curr[a]->SetDirectory(0);

	h1_ew_next[a] = new TH1D(Form("En_weighted_next_A%02d",a+1),Form("En_weighted_next_A%02d",a+1),200,0,20);
        h1_ew_next[a]->GetXaxis()->SetTitle("[MeV] 100keV/bin");
        h1_ew_next[a]->SetLineColor(kRed+1);
        h1_ew_next[a]->SetDirectory(0);

	h2_discri[a] = new TH2D(Form("discri_A%02d",a+1),Form("discri_A%02d",a+1),2500,0,500000,300,0,3);
	h2_discri[a]->SetDirectory(0);
    }


    // ===========================================================================================
    // === INPUT DATA
    TChain *ch = new TChain("EpicPretreatData");
    ch->Add(Form("../../../output/pretreat/pretreat_run17_run31.root")); 
    EpicPretreatData epic(ch);
    nentries = ch->GetEntries();
    cout << "nentries = " << nentries << endl;

    // ===========================================================================================
    // === LOOP OVER THE DATA
    for(long long int entry = 0 ; entry < nentries ; entry++){
	epic.GetEntry(entry);
    	if ((entry % 5000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;
	if(epic.fFC_q3>0 && epic.fFC_anode>0){ 
		h2_discri[epic.fFC_anode-1]->Fill(epic.fFC_q1,epic.fFC_q2/epic.fFC_q3);
	}

	if(epic.fFC_fission == true){
	   tof_raw_curr = epic.fFC_time - epic.fHF_time_curr - 2500000.;
	   tof_raw_next = epic.fFC_time - epic.fHF_time_next;
	   e_curr = TofRaw2Ene(gamma_peak_curr[epic.fFC_anode-1],FLIGHT_PATH + flight_path[epic.fFC_anode-1],tof_raw_curr,tof_cal_curr); 
	   e_next = TofRaw2Ene(gamma_peak_next[epic.fFC_anode-1],FLIGHT_PATH + flight_path[epic.fFC_anode-1],tof_raw_next,tof_cal_next); 


	   h1_tofraw_curr[epic.fFC_anode-1]->Fill(tof_raw_curr);
	   h1_tofraw_next[epic.fFC_anode-1]->Fill(tof_raw_next);
	   if((epic.fFC_time-epic.fHF_time_curr) < 2000.)
	   	h1_tofraw_next_extraWindow[epic.fFC_anode-1]->Fill(epic.fFC_time - epic.fHF_time_curr);
	   else
	   	h1_tofraw_next_extraWindow[epic.fFC_anode-1]->Fill(epic.fFC_time - epic.fHF_time_next);
		
	   h1_tofcal_curr[epic.fFC_anode-1]->Fill(tof_cal_curr);
	   h1_tofcal_next[epic.fFC_anode-1]->Fill(tof_cal_next);
	   h1_e_curr[epic.fFC_anode-1]->Fill(e_curr);
	   h1_e_next[epic.fFC_anode-1]->Fill(e_next);
	   h1_ew_curr[epic.fFC_anode-1]->Fill(e_curr , 1./actinide_mass[epic.fFC_anode-1]);
	   h1_ew_next[epic.fFC_anode-1]->Fill(e_next , 1./actinide_mass[epic.fFC_anode-1]);
	}


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
        h1_tofraw_next_extraWindow[a-1]->Draw();
        h1_tofraw_next[a-1]->Draw("same");
    }

    TCanvas * can3 = new TCanvas("TofCal","TofCal",0,0,3000,2000);
    can3->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can3->cd(a); gPad->SetLogy();
        h1_tofcal_next[a-1]->Draw();
        h1_tofcal_curr[a-1]->Draw("same");
    }

    TCanvas * can4 = new TCanvas("Ene","Ene",0,0,3000,2000);
    can4->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can4->cd(a); gPad->SetLogy();
        h1_e_next[a-1]->Draw();
        h1_e_curr[a-1]->Draw("same");
    }

    TCanvas * can5 = new TCanvas("Ene_weighted","Ene_weighted",0,0,3000,2000);
    can5->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can5->cd(a); gPad->SetLogy();
        h1_ew_next[a-1]->Draw();
        h1_ew_curr[a-1]->Draw("same");
    }

    // histograms per actinide

    TH1D * h1_e_U238_curr = new TH1D("En_U238_curr","En_U238_curr",200,0,20);
    TH1D * h1_e_U235_curr = new TH1D("En_U235_curr","En_U235_curr",200,0,20);
    TH1D * h1_e_Pu242_curr = new TH1D("En_Pu242_curr","En_Pu242_curr",200,0,20);
    h1_e_U238_curr->SetDirectory(0);
    h1_e_U235_curr->SetDirectory(0);
    h1_e_Pu242_curr->SetDirectory(0);
    h1_e_U238_curr->SetLineColor(kBlue+1);
    h1_e_U235_curr->SetLineColor(kBlue+1);
    h1_e_Pu242_curr->SetLineColor(kBlue+1);

    TH1D * h1_e_U238_next = new TH1D("En_U238_next","En_U238_next",200,0,20);
    TH1D * h1_e_U235_next = new TH1D("En_U235_next","En_U235_next",200,0,20);
    TH1D * h1_e_Pu242_next = new TH1D("En_Pu242_next","En_Pu242_next",200,0,20);
    h1_e_U238_next->SetDirectory(0);
    h1_e_U235_next->SetDirectory(0);
    h1_e_Pu242_next->SetDirectory(0);
    h1_e_U238_next->SetLineColor(kRed+1);
    h1_e_U235_next->SetLineColor(kRed+1);
    h1_e_Pu242_next->SetLineColor(kRed+1);
  
    h1_e_U238_curr->Add(h1_e_curr[0]);     h1_e_U238_next->Add(h1_e_next[0]); 
    h1_e_U235_curr->Add(h1_e_curr[1]);     h1_e_U235_next->Add(h1_e_next[1]); 
    h1_e_U235_curr->Add(h1_e_curr[2]);     h1_e_U235_next->Add(h1_e_next[2]); 
    h1_e_U235_curr->Add(h1_e_curr[3]);     h1_e_U235_next->Add(h1_e_next[3]); 
    h1_e_Pu242_curr->Add(h1_e_curr[4]);    h1_e_Pu242_next->Add(h1_e_next[4]); 
    h1_e_Pu242_curr->Add(h1_e_curr[5]);    h1_e_Pu242_next->Add(h1_e_next[5]); 
    h1_e_Pu242_curr->Add(h1_e_curr[6]);    h1_e_Pu242_next->Add(h1_e_next[6]); 
    h1_e_U235_curr->Add(h1_e_curr[7]);     h1_e_U235_next->Add(h1_e_next[7]); 
    h1_e_U235_curr->Add(h1_e_curr[8]);     h1_e_U235_next->Add(h1_e_next[8]); 
    h1_e_U235_curr->Add(h1_e_curr[9]);     h1_e_U235_next->Add(h1_e_next[9]); 

    TH1D * h1_ew_U238_curr = new TH1D("En_w_U238_curr","En_w_U238_curr",200,0,20);
    TH1D * h1_ew_U235_curr = new TH1D("En_w_U235_curr","En_w_U235_curr",200,0,20);
    TH1D * h1_ew_Pu242_curr = new TH1D("En_w_Pu242_curr","En_w_Pu242_curr",200,0,20);
    TH1D * h1_ew_U238_next = new TH1D("En_w_U238_next","En_w_U238_next",200,0,20);
    TH1D * h1_ew_U235_next = new TH1D("En_w_U235_next","En_w_U235_next",200,0,20);
    TH1D * h1_ew_Pu242_next = new TH1D("En_w_Pu242_next","En_w_Pu242_next",200,0,20);
    h1_ew_U238_curr->SetDirectory(0);
    h1_ew_U235_curr->SetDirectory(0);
    h1_ew_Pu242_curr->SetDirectory(0);
    h1_ew_U238_curr->SetLineColor(kBlue+1);
    h1_ew_U235_curr->SetLineColor(kBlue+1);
    h1_ew_Pu242_curr->SetLineColor(kBlue+1);
    h1_ew_U238_next->SetDirectory(0);
    h1_ew_U235_next->SetDirectory(0);
    h1_ew_Pu242_next->SetDirectory(0);
    h1_ew_U238_next->SetLineColor(kRed+1);
    h1_ew_U235_next->SetLineColor(kRed+1);
    h1_ew_Pu242_next->SetLineColor(kRed+1);
    for(short bin = 1 ; bin <= h1_ew_U238_curr->GetNbinsX(); bin++){
       h1_ew_U238_curr->SetBinContent(bin,h1_e_U238_curr->GetBinContent(bin)*1./MASS_U238);
       h1_ew_U235_curr->SetBinContent(bin,h1_e_U235_curr->GetBinContent(bin)*1./MASS_U235);
       h1_ew_Pu242_curr->SetBinContent(bin,h1_e_Pu242_curr->GetBinContent(bin)*1./MASS_PU242);
       h1_ew_U238_next->SetBinContent(bin,h1_e_U238_next->GetBinContent(bin)*1./MASS_U238);
       h1_ew_U235_next->SetBinContent(bin,h1_e_U235_next->GetBinContent(bin)*1./MASS_U235);
       h1_ew_Pu242_next->SetBinContent(bin,h1_e_Pu242_next->GetBinContent(bin)*1./MASS_PU242);
    }

    TH1D * h1_eratio_U238_curr = new TH1D("En_ratioU235_U238_curr","En_ratioU235_U238_curr",200,0,20);
    TH1D * h1_eratio_U235_curr = new TH1D("En_ratioU235_U235_curr","En_ratioU235_U235_curr",200,0,20);
    TH1D * h1_eratio_Pu242_curr = new TH1D("En_ratioU235_Pu242_curr","En_ratioU235_Pu242_curr",200,0,20);
    h1_eratio_U238_curr->SetDirectory(0);
    h1_eratio_U235_curr->SetDirectory(0);
    h1_eratio_Pu242_curr->SetDirectory(0);
    h1_eratio_U238_curr->SetLineColor(kBlue+1);
    h1_eratio_U235_curr->SetLineColor(kBlue+1);
    h1_eratio_Pu242_curr->SetLineColor(kBlue+1);
    h1_eratio_U238_curr->Add(h1_ew_U238_curr);
    h1_eratio_U235_curr->Add(h1_ew_U235_curr);
    h1_eratio_Pu242_curr->Add(h1_ew_Pu242_curr);
    h1_eratio_U238_curr->Divide(h1_ew_U235_curr);
    h1_eratio_U235_curr->Divide(h1_ew_U235_curr);
    h1_eratio_Pu242_curr->Divide(h1_ew_U235_curr);

    TH1D * h1_eratio_U238_next = new TH1D("En_ratioU235_U238_next","En_ratioU235_U238_next",200,0,20);
    TH1D * h1_eratio_U235_next = new TH1D("En_ratioU235_U235_next","En_ratioU235_U235_next",200,0,20);
    TH1D * h1_eratio_Pu242_next = new TH1D("En_ratioU235_Pu242_next","En_ratioU235_Pu242_next",200,0,20);
    h1_eratio_U238_next->SetDirectory(0);
    h1_eratio_U235_next->SetDirectory(0);
    h1_eratio_Pu242_next->SetDirectory(0);
    h1_eratio_U238_next->SetLineColor(kRed+1);
    h1_eratio_U235_next->SetLineColor(kRed+1);
    h1_eratio_Pu242_next->SetLineColor(kRed+1);
    h1_eratio_U238_next->Add(h1_ew_U238_next);
    h1_eratio_U235_next->Add(h1_ew_U235_next);
    h1_eratio_Pu242_next->Add(h1_ew_Pu242_next);
    h1_eratio_U238_next->Divide(h1_ew_U235_next);
    h1_eratio_U235_next->Divide(h1_ew_U235_next);
    h1_eratio_Pu242_next->Divide(h1_ew_U235_next);

    TCanvas * can6 = new TCanvas("Ene_per_actinide","Ene_per_actinide",0,0,3000,2000);
    can6->Divide(3,3);

    can6->cd(1); h1_e_U238_next->Draw();  h1_e_U238_curr->Draw("same");
    can6->cd(2); h1_e_U235_next->Draw();  h1_e_U235_curr->Draw("same");
    can6->cd(3); h1_e_Pu242_next->Draw(); h1_e_Pu242_curr->Draw("same");
    can6->cd(4); h1_ew_U238_next->Draw();  h1_ew_U238_curr->Draw("same");
    can6->cd(5); h1_ew_U235_next->Draw();  h1_ew_U235_curr->Draw("same");
    can6->cd(6); h1_ew_Pu242_next->Draw(); h1_ew_Pu242_curr->Draw("same");
    can6->cd(7); h1_eratio_U238_next->Draw();  h1_eratio_U238_curr->Draw("same");
    can6->cd(8); h1_eratio_U235_next->Draw();  h1_eratio_U235_curr->Draw("same");
    can6->cd(9); h1_eratio_Pu242_next->Draw(); h1_eratio_Pu242_curr->Draw("same");

    can6->cd(7); m_phy_gr["XSratio_JEFF_238U_235U"]->Draw("sameLP");
    can6->cd(9); m_phy_gr["XSratio_JEFF_242Pu_235U"]->Draw("sameLP");
    can6->cd(7); m_phy_gr["XSratio_BVIII_238U_235U"]->Draw("sameLP");
    can6->cd(9); m_phy_gr["XSratio_BVIII_242Pu_235U"]->Draw("sameLP");

    
    TCanvas * can7 = new TCanvas("Ene_Pu242","Ene_Pu242",0,0,3000,2000);
    can7->Divide(2,2);

    // histograms per anode of Pu242
    TH1D * h1_eratio_Pu242_per_anode_curr[3];
    TH1D * h1_eratio_Pu242_per_anode_next[3];
    for(short a = 1 ; a <= 3 ; a++){
    	h1_eratio_Pu242_per_anode_curr[a-1] = new TH1D(Form("En_ratio_U235_Pu242_A%02d_curr",a+4),Form("En_ratio_U235_Pu242_A%02d_curr",a+4),200,0,20);
    	h1_eratio_Pu242_per_anode_next[a-1] = new TH1D(Form("En_ratio_U235_Pu242_A%02d_next",a+4),Form("En_ratio_U235_Pu242_A%02d_next",a+4),200,0,20);
    	h1_eratio_Pu242_per_anode_curr[a-1]->SetDirectory(0);
    	h1_eratio_Pu242_per_anode_next[a-1]->SetDirectory(0);
    	h1_eratio_Pu242_per_anode_curr[a-1]->SetLineColor(kBlue+1);
    	h1_eratio_Pu242_per_anode_next[a-1]->SetLineColor(kRed+1);
    	h1_eratio_Pu242_per_anode_curr[a-1]->Add(h1_ew_curr[a+3]);
    	h1_eratio_Pu242_per_anode_next[a-1]->Add(h1_ew_next[a+3]);
    	h1_eratio_Pu242_per_anode_curr[a-1]->Divide(h1_ew_U235_curr);
    	h1_eratio_Pu242_per_anode_next[a-1]->Divide(h1_ew_U235_next);
        can7->cd(a); gPad->SetGridx(); gPad->SetGridy(); 
        h1_eratio_Pu242_per_anode_next[a-1]->Draw(); 
        h1_eratio_Pu242_per_anode_curr[a-1]->Draw("same");
        m_phy_gr["XSratio_JEFF_242Pu_235U"]->Draw("sameLP");
    }
    can7->cd(4); h1_eratio_Pu242_next->Draw();
    m_phy_gr["XSratio_JEFF_242Pu_235U"]->Draw("sameLP");



    TCanvas * can8 = new TCanvas("Discri","Discri",0,0,3000,2000);
    can8->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can8->cd(a); gPad->SetLogz();
        h2_discri[a-1]->Draw("col");
        cutfission[a-1]->Draw("same");
    }


   

}
