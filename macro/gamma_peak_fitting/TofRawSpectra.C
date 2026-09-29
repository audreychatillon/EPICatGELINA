#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

//#include "/home/audrey/.local/nptool/default/include/EpicData.h"
#include "/home/chatillona/.local/nptool/default/include/EpicData.h"

using namespace std;

void run()
{

    // === =========================================================
    // === variables 

    char   name[100];
    double thf_curr = -1;

    vector<short>   *fFC_DetNbr = nullptr;

    struct TofRawInfo{
        short  anode;
        double tFC ;
        double tHF_curr;
    };
    vector<TofRawInfo> pendingFC;
   
    double tofraw_curr_offset = -2500000.;
    double tofraw_next_offset = 0.;




    // === =========================================================
    // === histograms 

    TH1D * h1_tofraw_curr[11];
    TH1D * h1_tofraw_next[11];
    for(short a = 1 ; a <= 11 ; a++){
        sprintf(name,"tofraw_curr_A%02d",a);
        h1_tofraw_curr[a-1] = new TH1D(name,name,20000,-1500,500);
        h1_tofraw_curr[a-1]->SetLineColor(kBlue);
        h1_tofraw_curr[a-1]->SetDirectory(0);

        sprintf(name,"tofraw_next_A%02d",a);
        h1_tofraw_next[a-1] = new TH1D(name,name,20000,-1500,500);
        h1_tofraw_next[a-1]->SetLineColor(kRed);
        h1_tofraw_next[a-1]->SetDirectory(0);
    }



    // === =========================================================
    // === input data 
  
    //TFile * f = new TFile(Form("../../output/conversion/raw%i.root",run_number),"read");
    //TTree * tFC = (TTree*)f->Get("EpicRawTree");
    TChain * tFC = new TChain("EpicRawTree");
    // first part with V4B
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
    tFC->Add("../../output/conversion/raw_run32.root");
    tFC->Add("../../output/conversion/raw_run33.root");
    tFC->Add("../../output/conversion/raw_run34.root");
    tFC->Add("../../output/conversion/raw_run35.root");
    tFC->Add("../../output/conversion/raw_run36.root");
    tFC->Add("../../output/conversion/raw_run37.root");
    tFC->Add("../../output/conversion/raw_run38.root");
    tFC->Add("../../output/conversion/raw_run39.root");
    tFC->Add("../../output/conversion/raw_run40.root");

    // === =========================================================
    // === branches 

    epic::EpicData *epicFC = nullptr;
    int statusFC = tFC->SetBranchAddress("epic",&epicFC);


    // === =========================================================
    // === loop 
    ULong64_t nentries = tFC->GetEntries();
    cout << "nentries = " << nentries << endl;
    for(ULong64_t entry=0; entry < nentries ; entry++){
    //for(ULong64_t entry=0; entry < 15000000 ; entry++){

        if ((entry % 500000) == 0)   cout << "\r === Entry = " << entry << " / " << nentries << " === " << flush;
        

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
                h1_tofraw_curr[tof.anode-1]->Fill(tof.tFC - tof.tHF_curr + tofraw_curr_offset);
                h1_tofraw_next[tof.anode-1]->Fill(tof.tFC - thf_next + tofraw_next_offset);
            }

            // Clear
            pendingFC.clear();

            // newvalue
            thf_curr = epicFC->GetTimeHF();
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
            fc.tHF_curr = thf_curr;
            pendingFC.push_back(fc);
        }
    }// end of loop over the entries
    

    TCanvas * can = new TCanvas("TofRawComparison","TofRawComparison",0,0,3000,2000);
    can->Divide(4,3);
    for(short a = 1 ; a <=11 ; a++){
        can->cd(a);
        h1_tofraw_next[a-1]->Draw();
        h1_tofraw_curr[a-1]->Draw("same");

    }


}
