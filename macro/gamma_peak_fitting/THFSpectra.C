#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

//#include "/home/audrey/.local/nptool/default/include/EpicData.h"
#include "/home/chatillona/.local/nptool/default/include/EpicData.h"

using namespace std;

void run()
{

    // === =========================================================
    // === variables 

    char   name[100];

    // === for output tree

    bool   fHF_data;
    double fHF_time;
    double fHF_dthf;

    bool   fFC_data;
    short  fFC_anode;
    short  fFC_time;
    double fFC_qmax;
    double fFC_q1;
    double fFC_q2;
    double fFC_q3;
    bool   fFC_alpha;
    bool   fFC_fission;

    // === output tree
    TFile * fout = new TFile("../../output/pretreat/pretreat_run18.root","recreate");
    TTree * tout = new TTree("EpicPretreatData","HF or FC, no cross talk, no reflexion");
    tout->Branch("fHF_data",&fHF_data);
    tout->Branch("fHF_time",&fHF_time);
    tout->Branch("fHF_dthf",&fHF_dthf);
    tout->Branch("fFC_data",&fFC_data);
    tout->Branch("fFC_anode",&fFC_anode);
    tout->Branch("fFC_time",&fFC_time);
    tout->Branch("fFC_qmax",&fFC_qmax);
    tout->Branch("fFC_q1",&fFC_q1);
    tout->Branch("fFC_q2",&fFC_q2);
    tout->Branch("fFC_q3",&fFC_q3);
    tout->Branch("fFC_alpha",&fFC_alpha);
    tout->Branch("fFC_fission",&fFC_fission);

    // === =========================================================
    // === TCutG for alpha and fission
    TFile * fcutg = new TFile("2DdiscriTCutG.root","read");
    TCutG * cutalpha[10];
    TCutG * cutfission[10];
    for(short a = 0 ; a < 10 ; a++){
	cutalpha[a] = (TCutG*)fcutg->Get(Form("A%02d_alpha",a+1));
	cutfission[a] = (TCutG*)fcutg->Get(Form("A%02d_fission",a+1));
    } 
    fcutg->Close();
 
    // === =========================================================
    // === histograms to understand the HF 

    TH1D * h1_THF = new TH1D("THF","THF",86400,0,86400);
    TH1D * h1_DeltaTHF = new TH1D("DeltaTHF","DeltaTHF",5100,0,5100000);
    TH1D * h1_DeltaTHFzoom = new TH1D("DeltaTHFzoom","DeltaTHFzoom",20000,2499000,2501000);
    TH1D * h1_noR_THF = new TH1D("THF_noReflex","THF_noReflex",86400,0,86400);
    TH1D * h1_noR_DeltaTHF = new TH1D("DeltaTHF_noReflex","DeltaTHF_noReflex",5100,0,5100000);
    TH1D * h1_noR_HFonly_DeltaTHF = new TH1D("DeltaTHFzoom_noReflex_HFonly","DeltaTHFzoom_noReflex_HFonly",5100,0,5100000);
    TH1D * h1_noR_DeltaTHFzoom = new TH1D("DeltaTHFzoom_noReflex","DeltaTHFzoom_noReflex",20000,2499000,2501000);
    TH1D * h1_noR_Good_THF = new TH1D("THF_noReflex_GoodDeltaT","THF_noReflex_GoodDeltaT",86400,0,86400);

    h1_THF->GetXaxis()->SetTitle("Time [s]");
    h1_DeltaTHF->GetXaxis()->SetTitle("Time [ns]");
    h1_DeltaTHFzoom->GetXaxis()->SetTitle("Time [ns], 100 ps/bin");

    h1_noR_THF->SetLineColor(kRed);
    h1_noR_DeltaTHF->SetLineColor(kRed);
    h1_noR_DeltaTHFzoom->SetLineColor(kRed);

    h1_noR_Good_THF->SetLineColor(8);
    h1_noR_Good_THF->SetLineStyle(2);

    h1_noR_HFonly_DeltaTHF->SetLineColor(kOrange+7);
    h1_noR_HFonly_DeltaTHF->SetLineStyle(3);
    h1_noR_HFonly_DeltaTHF->SetLineWidth(3);

    h1_THF->SetDirectory(0);
    h1_DeltaTHF->SetDirectory(0);
    h1_DeltaTHFzoom->SetDirectory(0);
    h1_noR_THF->SetDirectory(0);
    h1_noR_DeltaTHF->SetDirectory(0);
    h1_noR_HFonly_DeltaTHF->SetDirectory(0);
    h1_noR_DeltaTHFzoom->SetDirectory(0);
    h1_noR_Good_THF->SetDirectory(0);

    // === =========================================================
    // === 2D-discri
    TH2D * h2_discri[10];
    for(short a = 0 ; a < 10 ; a++){
	h2_discri[a] = new TH2D(Form("discri_A%02d",a+1),Form("discri_A%02d",a+1),3000,0,300000,300,0,3);
	h2_discri[a]->SetDirectory(0);
    } 
    

    // === =========================================================
    // === input data 
  
    //TFile * f = new TFile(Form("../../output/conversion/raw%i.root",run_number),"read");
    //TTree * tFC = (TTree*)f->Get("EpicRawTree");
    TChain * tFC = new TChain("EpicRawTree");
    // first part with V4B
    tFC->Add("../../output/conversion/test18.root"); // raw data with reflexion
    //tFC->Add("../../output/conversion/raw_run18.root"); // raw data without reflexion
    //tFC->Add("../../output/conversion/raw_run17.root");
    //tFC->Add("../../output/conversion/raw_run18.root");
    //tFC->Add("../../output/conversion/raw_run19.root");
    //tFC->Add("../../output/conversion/raw_run21.root");
    //tFC->Add("../../output/conversion/raw_run22.root");
    //tFC->Add("../../output/conversion/raw_run23.root");
    //tFC->Add("../../output/conversion/raw_run24.root");
    //tFC->Add("../../output/conversion/raw_run25.root");
    //tFC->Add("../../output/conversion/raw_run26.root");
    //tFC->Add("../../output/conversion/raw_run27.root");
    //tFC->Add("../../output/conversion/raw_run28.root");
    //tFC->Add("../../output/conversion/raw_run29.root");
    //tFC->Add("../../output/conversion/raw_run30.root");
    //tFC->Add("../../output/conversion/raw_run31.root");
    // second part with V4B and new parameters for HF channel
    //tFC->Add("../../output/conversion/raw_run32.root");
    //tFC->Add("../../output/conversion/raw_run33.root");
    //tFC->Add("../../output/conversion/raw_run34.root");
    //tFC->Add("../../output/conversion/raw_run35.root");
    //tFC->Add("../../output/conversion/raw_run36.root");
    //tFC->Add("../../output/conversion/raw_run37.root");
    //tFC->Add("../../output/conversion/raw_run38.root");
    //tFC->Add("../../output/conversion/raw_run39.root");
    //tFC->Add("../../output/conversion/raw_run40.root");

    // === =========================================================
    // === branches 

    epic::EpicData *epicFC = nullptr;
    int statusFC = tFC->SetBranchAddress("epic",&epicFC);


    // === =========================================================
    // === loop 
    ULong64_t nentries = tFC->GetEntries();
    cout << "nentries = " << nentries << endl;
    double thf_curr = -1;
    double thf_prev = -1;
    for(ULong64_t entry=0; entry < nentries ; entry++){

        if ((entry % 1000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;

        int bytes = tFC->GetEntry(entry);
        if(!epicFC) continue;

        int mult = epicFC->GetFCMult() ;

        // --- skip empty entry 
        if (mult==0) continue;

	// --- is there an HF in FCdata?
        short ihf  = epicFC->GetHFIndex();

	// --- is there an FC in FCdata?
        short ifc  = epicFC->GetQmaxIndex();

	// --- init data for TTree
	fHF_data = false;
        fHF_time = -1;
        fHF_dthf = -1;
	fFC_data = false;
	fFC_anode = -1;
	fFC_time = -1;
	fFC_qmax = -1;
	fFC_q1 = -1;
	fFC_q2 = -1;
	fFC_q3 = -1;
	fFC_alpha = false;
	fFC_fission = false;

        if(ihf >= 0 && ifc == -1){
		h1_THF->Fill(epicFC->GetTimeHF()*1.e-09);
		h1_DeltaTHF->Fill(epicFC->GetDeltaTHF()); 
		h1_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 
        	if(epicFC->GetDeltaTHF()>1000){
			thf_prev = thf_curr;
			thf_curr = epicFC->GetTimeHF();
			h1_noR_THF->Fill(thf_curr*1.e-09);
			h1_noR_DeltaTHF->Fill(thf_curr-thf_prev); 
			if(ifc==-1) h1_noR_HFonly_DeltaTHF->Fill(thf_curr-thf_prev); 
			h1_noR_DeltaTHFzoom->Fill(thf_curr-thf_prev); 
			if(2500005 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500035){
				h1_noR_Good_THF->Fill(thf_curr*1.e-09); 
				fHF_data = true;
				fHF_time = thf_curr;
				fHF_dthf = thf_curr - thf_prev;
				tout->Fill();
				// reinit
				fHF_data = false;
				fHF_time = -1;
				fHF_dthf = -1;
			}
		}
 	}

	if(ifc >= 0 && ihf == -1){
		fFC_data = true;
		fFC_anode = epicFC->GetAnodeNbr(ifc);
		fFC_time = epicFC->GetTimeFC(ifc);
		fFC_qmax = epicFC->GetQmax(ifc);
		fFC_q1 = epicFC->GetQ1(ifc);
		fFC_q2 = epicFC->GetQ2(ifc);
		fFC_q3 = epicFC->GetQ3(ifc);
		if(fFC_q3>0 && fFC_anode>0){ 
			h2_discri[fFC_anode-1]->Fill(fFC_q1,fFC_q2/fFC_q3);
                	if(cutalpha[fFC_anode-1]->IsInside(fFC_q1,fFC_q2/fFC_q3)) fFC_alpha = true;
			if(cutfission[fFC_anode-1]->IsInside(fFC_q1,fFC_q2/fFC_q3)) fFC_fission = true; 
		}
		tout->Fill();
		// reinit
		fFC_data = false;
		fFC_anode = -1;
		fFC_time = -1;
		fFC_qmax = -1;
		fFC_q1 = -1;
		fFC_q2 = -1;
		fFC_q3 = -1;
		fFC_alpha = false;
		fFC_fission = false;
	}

	//TODO if(ifc >=0 && ihf >= 0)

    }// end of loop over the entries
    
    TCanvas * can1 = new TCanvas("HF","HF",0,0,2000,3000);
    can1->Divide(1,3);
    can1->cd(1); h1_THF->Draw(); h1_noR_THF->Draw("sames"); h1_noR_Good_THF->Draw("sames"); 
    can1->cd(2); h1_DeltaTHF->Draw(); h1_noR_DeltaTHF->Draw("sames"); h1_noR_HFonly_DeltaTHF->Draw("sames");
    can1->cd(3); h1_DeltaTHFzoom->Draw(); h1_noR_DeltaTHFzoom->Draw("sames");

    TCanvas * can2[10];
    for (short a = 0 ; a <10 ; a++){
	can2[a] = new TCanvas(Form("Q2Q3vQ1_A%02d",a+1),Form("Q2Q3vQ1_A%02d",a+1),0,0,2000,1500);
	can2[a]->cd();
	h2_discri[a]->Draw("col"); cutalpha[a]->Draw("same"); cutfission[a]->Draw("same");
    }

    fout->cd();
    tout->Write();
    fout->Close();

}
