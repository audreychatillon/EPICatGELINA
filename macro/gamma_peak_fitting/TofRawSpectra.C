#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "/home/audrey/.local/nptool/default/include/EpicData.h"
//#include "/home/chatillona/.local/nptool/default/include/EpicData.h"

using namespace std;

void run()
{

    // === =========================================================
    // === variables 

    char   name[100];
    double thf_prev = -1;
    double thf_curr = -1;

    vector<short>   *fFC_DetNbr = nullptr;

    struct TofRawInfo{
        short  anode;
        double tFC ;
        double tHF_prev;
        double tHF_curr;
    };
    vector<TofRawInfo> pendingFC;
    vector<TofRawInfo> pendingFCx;
   

    // === =========================================================
    // === histograms 

    TH1D * h1_tofraw_curr_zoom[11];
    TH1D * h1_tofraw_next_zoom[11];
    TH1D * h1_tofraw_curr[11];
    TH1D * h1_tofraw_next[11];
    TH1D * h1_tofraw_curr_offset[11];
    TH1D * h1_tofraw_next_offset[11];
    TH1D * h1_tofraw_wExtraWindow_curr_zoom[11];
    TH1D * h1_tofraw_wExtraWindow_next_zoom[11];
    for(short a = 1 ; a <= 11 ; a++){
        sprintf(name,"gnf_tofraw_curr_zoom_A%02d",a);
        h1_tofraw_curr_zoom[a-1] = new TH1D(name,name,30000,2498500,2501500);
        h1_tofraw_curr_zoom[a-1]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_curr_zoom[a-1]->SetLineColor(kBlue);
        h1_tofraw_curr_zoom[a-1]->SetDirectory(0);

        sprintf(name,"gnf_tofraw_next_zoom_A%02d",a);
        h1_tofraw_next_zoom[a-1] = new TH1D(name,name,30000,-1500,1500);
        h1_tofraw_next_zoom[a-1]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_next_zoom[a-1]->SetLineColor(kRed);
        h1_tofraw_next_zoom[a-1]->SetDirectory(0);

        sprintf(name,"gnf_tofraw_curr_A%02d",a);
        h1_tofraw_curr[a-1] = new TH1D(name,name,50400,-10000,2510000);
        h1_tofraw_curr[a-1]->GetXaxis()->SetTitle("[ns] 50ns/bin");
        h1_tofraw_curr[a-1]->SetLineColor(kBlue);
        h1_tofraw_curr[a-1]->SetDirectory(0);

        sprintf(name,"gnf_tofraw_next_zoom_A%02d",a);
        h1_tofraw_next[a-1] = new TH1D(name,name,50400,-2510000,10000);
        h1_tofraw_next[a-1]->GetXaxis()->SetTitle("[ns] 50ns/bin");
        h1_tofraw_next[a-1]->SetLineColor(kRed);
        h1_tofraw_next[a-1]->SetDirectory(0);

        sprintf(name,"gnf_tofraw_curr_offset_A%02d",a);
        h1_tofraw_curr_offset[a-1] = new TH1D(name,name,30000,0,3000);
        h1_tofraw_curr_offset[a-1]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_curr_offset[a-1]->SetLineColor(kBlue);
        h1_tofraw_curr_offset[a-1]->SetDirectory(0);

        sprintf(name,"gnf_tofraw_next_offset_zoom_A%02d",a);
        h1_tofraw_next_offset[a-1] = new TH1D(name,name,30000,-2501000,-2498000);
        h1_tofraw_next_offset[a-1]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_next_offset[a-1]->SetLineColor(kRed);
        h1_tofraw_next_offset[a-1]->SetDirectory(0);

        sprintf(name,"gnf_tofraw_curr_zoom_A%02d",a);
        h1_tofraw_wExtraWindow_curr_zoom[a-1] = new TH1D(name,name,30000,2498500,2501500);
        h1_tofraw_wExtraWindow_curr_zoom[a-1]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_wExtraWindow_curr_zoom[a-1]->SetLineColor(kBlue);
        h1_tofraw_wExtraWindow_curr_zoom[a-1]->SetDirectory(0);

        sprintf(name,"gnf_tofraw_next_zoom_A%02d",a);
        h1_tofraw_wExtraWindow_next_zoom[a-1] = new TH1D(name,name,30000,-1500,1500);
        h1_tofraw_wExtraWindow_next_zoom[a-1]->GetXaxis()->SetTitle("[ns] 100ps/bin");
        h1_tofraw_wExtraWindow_next_zoom[a-1]->SetLineColor(kRed);
        h1_tofraw_wExtraWindow_next_zoom[a-1]->SetDirectory(0);
    }



    // === =========================================================
    // === input data 
  
    //TFile * f = new TFile(Form("../../output/conversion/raw%i.root",run_number),"read");
    //TTree * tFC = (TTree*)f->Get("EpicRawTree");
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
    //for(ULong64_t entry=0; entry < nentries ; entry++){
    for(ULong64_t entry=0; entry < 50000000 ; entry++){

        if ((entry % 1000000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;
        

        int bytes = tFC->GetEntry(entry);
        if(!epicFC) continue;

        int mult = epicFC->GetFCMult() ;

        //--- skip empty entry 
        if (mult==0) continue;

        //--- get t_hf infos
        if(epicFC->GetDetNbr(0) == -1){
          
            double thf_next = epicFC->GetTimeHF();

            // Fill histograms if pendingFC has data
            for(const TofRawInfo &tof : pendingFC){
                h1_tofraw_curr_zoom[tof.anode-1]->Fill(tof.tFC - tof.tHF_curr);
                h1_tofraw_next_zoom[tof.anode-1]->Fill(tof.tFC - thf_next);
                h1_tofraw_curr[tof.anode-1]->Fill(tof.tFC - tof.tHF_curr);
                h1_tofraw_next[tof.anode-1]->Fill(tof.tFC - thf_next);
                h1_tofraw_curr_offset[tof.anode-1]->Fill(tof.tFC - tof.tHF_curr);
                h1_tofraw_next_offset[tof.anode-1]->Fill(tof.tFC - thf_next);
            }
            for(const TofRawInfo &tof : pendingFCx){
                h1_tofraw_wExtraWindow_curr_zoom[tof.anode-1]->Fill(tof.tFC - tof.tHF_curr);
                h1_tofraw_wExtraWindow_next_zoom[tof.anode-1]->Fill(tof.tFC - thf_next);
            }

            // Clear
            pendingFC.clear();
            pendingFCx.clear();

            // newvalue
            thf_curr = thf_next;
            thf_prev = thf_curr - epicFC->GetDeltaTHF();
            continue;
        }
        else{

            //--- process FC entry
            short index_qmax = epicFC->GetQmaxIndex();

            //    skip alpha decay
            if(!epicFC->GetIsFission(index_qmax))  continue;

            //    get FC data
            TofRawInfo fc;
            fc.anode    = epicFC->GetAnodeNbr(index_qmax);
            fc.tFC      = epicFC->GetTimeFC(index_qmax);
            fc.tHF_prev = thf_prev;
            fc.tHF_curr = thf_curr;
            pendingFC.push_back(fc);
            if((fc.tFC-fc.tHF_curr) < 1000){
                h1_tofraw_wExtraWindow_curr_zoom[fc.anode-1]->Fill(fc.tFC - fc.tHF_prev);
                h1_tofraw_wExtraWindow_next_zoom[fc.anode-1]->Fill(fc.tFC - fc.tHF_curr);
            }
            else
                pendingFCx.push_back(fc);
        }
    }// end of loop over the entries
    

    TCanvas * can1 = new TCanvas("TofRawCurrZoom","TofRawCurrZoom",0,0,3000,2000);
    can1->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can1->cd(a); gPad->SetLogy();
        h1_tofraw_curr_zoom[a-1]->Draw();
    }

    TCanvas * can2 = new TCanvas("TofRawNextZoom","TofRawNextZoom",0,0,3000,2000);
    can2->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can2->cd(a); gPad->SetLogy();
        h1_tofraw_next_zoom[a-1]->Draw();
    }

    TCanvas * can3 = new TCanvas("TofRawCurr","TofRawCurr",0,0,3000,2000);
    can3->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can3->cd(a); gPad->SetLogy();
        h1_tofraw_curr[a-1]->Draw();
    }

    TCanvas * can4 = new TCanvas("TofRawNext","TofRawNext",0,0,3000,2000);
    can4->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can4->cd(a); gPad->SetLogy();
        h1_tofraw_next[a-1]->Draw();
    }

    TCanvas * can5 = new TCanvas("TofRawCurrOffset","TofRawCurrOffset",0,0,3000,2000);
    can5->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can5->cd(a); gPad->SetLogy();
        h1_tofraw_curr_offset[a-1]->Draw();
    }

    TCanvas * can6 = new TCanvas("TofRawNextOffset","TofRawNextOffset",0,0,3000,2000);
    can6->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can6->cd(a); gPad->SetLogy();
        h1_tofraw_next_offset[a-1]->Draw();
    }


    TCanvas * can7 = new TCanvas("TofRawCurrZoom_extraW","TofRawCurrZoom_extraW",0,0,3000,2000);
    can7->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can7->cd(a); gPad->SetLogy();
        h1_tofraw_wExtraWindow_curr_zoom[a-1]->Draw();
    }

    TCanvas * can8 = new TCanvas("TofRawNextZoom_extraW","TofRawNextZoom_extraW",0,0,3000,2000);
    can8->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can8->cd(a); gPad->SetLogy();
        h1_tofraw_wExtraWindow_next_zoom[a-1]->Draw();
    }


}
