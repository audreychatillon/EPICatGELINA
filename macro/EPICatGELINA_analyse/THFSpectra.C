#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "/home/audrey/.local/nptool/default/include/EpicData.h"
//#include "/home/chatillona/.local/nptool/default/include/EpicData.h"

using namespace std;

void run()
{
    // === =========================================================
    // === variables 

    char   name[100];

    // to calculate DeltaT removing reflexion
    double thf_curr = -1.;
    double thf_prev = -1.;

    // to keep track of last good HF data
    double thf_curr_good = -1.;

    // to check TofRaw
    struct FCdata{
        short  FC_anode;
    	double FC_qmax;
    	double FC_q1;
    	double FC_q2;
    	double FC_q3;
    	bool   FC_alpha;
    	bool   FC_fission;
        double FC_time;
        double HF_curr;
    };
    vector<FCdata> pendingFC;
    double tofraw_curr;
    double tofraw_next;

    // === =========================================================
    // === histograms to understand the HF 
    
    TH1D * h1_THF = new TH1D("THF","THF",86400,0,86400);
    TH1D * h1_DeltaTHF = new TH1D("DeltaTHF","DeltaTHF",5100,0,5100000);
    TH1D * h1_DeltaTHFzoom = new TH1D("DeltaTHFzoom","DeltaTHFzoom",20000,2499000,2501000);
    TH1D * h1_onlyR_THF = new TH1D("THF_onlyReflex","THF_onlyReflex",86400,0,86400);
    TH1D * h1_onlyR_DeltaTHFzoom = new TH1D("DeltaTHFzoom_onlyReflex","DeltaTHFzoom_onlyReflex",10000,850,950);
    TH1D * h1_noR_THF = new TH1D("THF_noReflex","THF_noReflex",86400,0,86400);
    TH1D * h1_noR_DeltaTHF = new TH1D("DeltaTHF_noReflex","DeltaTHF_noReflex",5100,0,5100000);
    TH1D * h1_noR_HFonly_DeltaTHF = new TH1D("DeltaTHFzoom_noReflex_HFonly","DeltaTHFzoom_noReflex_HFonly",5100,0,5100000);
    TH1D * h1_noR_DeltaTHFzoom = new TH1D("DeltaTHFzoom_noReflex","DeltaTHFzoom_noReflex",20000,2499000,2501000);
    TH1D * h1_noR_Good_THF = new TH1D("THF_noReflex_GoodDeltaT","THF_noReflex_GoodDeltaT",86400,0,86400);
    
    h1_THF->GetXaxis()->SetTitle("Time [s]");
    h1_DeltaTHF->GetXaxis()->SetTitle("Time [ns]");
    h1_DeltaTHFzoom->GetXaxis()->SetTitle("Time [ns], 100 ps/bin");
   
    h1_onlyR_THF->SetLineColor(kBlack);
    h1_onlyR_DeltaTHFzoom->SetLineColor(kBlack);

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
    // === TofRaw histograms
    TH1D * h1_TofRaw_curr[10];
    TH1D * h1_TofRaw_next[10];
    for(short a = 0 ; a < 10 ; a++){
        sprintf(name,"TofRaw_curr_A%02d",a+1);
        h1_TofRaw_curr[a] = new TH1D(name,name,30000,-2500,500);
        h1_TofRaw_curr[a]->SetDirectory(0);
        h1_TofRaw_curr[a]->SetLineColor(kBlue);
        sprintf(name,"TofRaw_next_A%02d",a+1);
        h1_TofRaw_next[a] = new TH1D(name,name,30000,-2500,500);
        h1_TofRaw_next[a]->SetDirectory(0);
        h1_TofRaw_next[a]->SetLineColor(kRed);

    }



    // === =========================================================
    // === input data 
    
    TChain * tFC = new TChain("EpicRawTree");
    // first part with V4B
    tFC->Add("../../output/conversion/raw_run17.root");
    tFC->Add("../../output/conversion/raw_run18.root");
    tFC->Add("../../output/conversion/raw_run19.root");
    tFC->Add("../../output/conversion/raw_run21.root");
    tFC->Add("../../output/conversion/raw_run22.root");
    tFC->Add("../../output/conversion/raw_run23.root");
    tFC->Add("../../output/conversion/raw_run24.root");
    tFC->Add("../../output/conversion/raw_run25.root");
    tFC->Add("../../output/conversion/raw_run26.root");
    tFC->Add("../../output/conversion/raw_run27.root");
    tFC->Add("../../output/conversion/raw_run28.root");
    tFC->Add("../../output/conversion/raw_run29.root");
    tFC->Add("../../output/conversion/raw_run30.root");
    tFC->Add("../../output/conversion/raw_run31.root");
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
    
        // --- HFdata only
        if(ihf >= 0 && ifc == -1){

            // all HF signals
    	    h1_THF->Fill(epicFC->GetTimeHF()*1.e-09);
    	    h1_DeltaTHF->Fill(epicFC->GetDeltaTHF()); 
    	    h1_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 
        
            // [noR] : remove relexion < 1000 ns
        	if(epicFC->GetDeltaTHF()>1000){
    		    thf_prev = thf_curr;
    		    thf_curr = epicFC->GetTimeHF();
    		    h1_noR_THF->Fill(thf_curr*1.e-09);
    		    h1_noR_DeltaTHF->Fill(thf_curr-thf_prev); 
    		    h1_noR_DeltaTHFzoom->Fill(thf_curr-thf_prev); 
    		    if(ifc==-1) h1_noR_HFonly_DeltaTHF->Fill(thf_curr-thf_prev); 

                // [noR_Good] : remove reflexion + keep only good DeltaT peak
    		    if(2500005 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500035){
    		    	h1_noR_Good_THF->Fill(thf_curr*1.e-09);

                    // calculate TofRaw if pendingFC and HFonly data
                    if(pendingFC.size() > 0){
                        for(size_t i = 0 ; i < pendingFC.size(); i++){
                            FCdata fc = pendingFC[i];
                            tofraw_curr = fc.FC_time - fc.HF_curr - 2500000.;
                            tofraw_next = fc.FC_time - thf_curr;
                            h1_TofRaw_curr[fc.FC_anode-1]->Fill(tofraw_curr);
                            h1_TofRaw_next[fc.FC_anode-1]->Fill(tofraw_next);
                        }
                        pendingFC.clear();
                    }
                    
                    // keep the current HF
                    thf_curr_good = thf_curr;
    		    }
    	    }// end of if reflexion
            else{
    		    h1_onlyR_THF->Fill(epicFC->GetTimeHF()*1.e-09);
    		    h1_onlyR_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 
            }
    	}

        // --- FCdata only
        else if(ifc >= 0 && ihf == -1){
            FCdata fc;
            fc.FC_anode = epicFC->GetAnodeNbr(ifc);
            fc.FC_time  = epicFC->GetTimeFC(ifc);
			fc.FC_qmax = epicFC->GetQmax(ifc);
			fc.FC_q1 = epicFC->GetQ1(ifc);
			fc.FC_q2 = epicFC->GetQ2(ifc);
			fc.FC_q3 = epicFC->GetQ3(ifc);
            fc.HF_curr  = thf_curr_good;
            pendingFC.push_back(fc);
        }

        // --- HFdata and FCdata
        else if(ifc >= 0 && ihf >= 0){
            // first FC data, then HF data and fill TofRaw
            // if HF and FC: HF in the slowest neutrons, good HF
            FCdata fc;
            fc.FC_anode = epicFC->GetAnodeNbr(ifc);
            fc.FC_time  = epicFC->GetTimeFC(ifc);
            fc.HF_curr  = thf_curr_good;
            pendingFC.push_back(fc);

            // all HF signals
    	    h1_THF->Fill(epicFC->GetTimeHF()*1.e-09);
    	    h1_DeltaTHF->Fill(epicFC->GetDeltaTHF()); 
    	    h1_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 

            // [noR] : remove relexion < 1000 ns
        	if(epicFC->GetDeltaTHF()>1000){
    		    thf_prev = thf_curr;
    		    thf_curr = epicFC->GetTimeHF();
    		    h1_noR_THF->Fill(thf_curr*1.e-09);
    		    h1_noR_DeltaTHF->Fill(thf_curr-thf_prev); 
    		    h1_noR_DeltaTHFzoom->Fill(thf_curr-thf_prev); 
    		    if(ifc==-1) h1_noR_HFonly_DeltaTHF->Fill(thf_curr-thf_prev); 

                // [noR_Good] : remove reflexion + keep only good DeltaT peak
    		    if(2500005 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500035){
    		    	h1_noR_Good_THF->Fill(thf_curr*1.e-09);
                    // calculate TofRaw if pendingFC and HFonly data
                    if(pendingFC.size() > 0){
                        for(size_t i = 0 ; i < pendingFC.size(); i++){
                            FCdata fc = pendingFC[i];
                            tofraw_curr = fc.FC_time - fc.HF_curr - 2500000.;
                            tofraw_next = fc.FC_time - thf_curr;
                            h1_TofRaw_curr[fc.FC_anode-1]->Fill(tofraw_curr);
                            h1_TofRaw_next[fc.FC_anode-1]->Fill(tofraw_next);
                        }
                        pendingFC.clear();
                    }
                    // keep the current HF
                    thf_curr_good = thf_curr;
    		    } // end of [noR_Good]

                // [curr = reflexion / next = Good]
                if(24199134 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2499148 ){

                }
        
                // [curr = Good / next = relexion]
                if(2500889 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500905){

                }


    	    }// end of if reflexion
            else{
    		    h1_onlyR_THF->Fill(epicFC->GetTimeHF()*1.e-09);
    		    h1_onlyR_DeltaTHFzoom->Fill(epicFC->GetDeltaTHF()); 
            }
        
        }

        // --- else ifc == -1 && ihf == -1 shoud never happen !!!
        else{
            cout << "ifc = " << ifc << " and ihf = " << ihf << " : event not processed" << endl; 
        }
    
    }// end of loop over the entries
    
    cout << endl;
    
    TCanvas * can1 = new TCanvas("HF","HF",0,0,2000,3000);
    can1->Divide(1,4);
    can1->cd(1); h1_THF->Draw(); h1_noR_THF->Draw("sames"); h1_noR_Good_THF->Draw("sames"); h1_onlyR_THF->Draw("sames"); 
    can1->cd(2); h1_DeltaTHF->Draw(); h1_noR_DeltaTHF->Draw("sames"); h1_noR_HFonly_DeltaTHF->Draw("sames");
    can1->cd(3); h1_DeltaTHFzoom->Draw(); h1_noR_DeltaTHFzoom->Draw("sames");
    can1->cd(4); h1_onlyR_DeltaTHFzoom->Draw(); 


    TCanvas * can2 = new TCanvas("TofRaw","TofRaw",2000,3000);
    can2->Divide(5,2);
    for(short a = 0 ; a < 10 ; a++){
        can2->cd(a+1); h1_TofRaw_next[a]->Draw(); h1_TofRaw_curr[a]->Draw("same");
    }

}
