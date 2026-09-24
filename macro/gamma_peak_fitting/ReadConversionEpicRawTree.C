#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include "/home/audrey/.local/nptool/default/include/EpicData.h"

using namespace std;

void run(UShort_t run_number)
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
   
    double tofraw_prev_offset = -5000000.+1000.;
    double tofraw_curr_offset = -2500000.+1000.;
    double tofraw_next_offset = 1000.;




    // === =========================================================
    // === histograms 

    TH1D * h1_tofraw_prev[11];
    TH1D * h1_tofraw_curr[11];
    TH1D * h1_tofraw_next[11];
    for(short a = 1 ; a <= 11 ; a++){
        sprintf(name,"tofraw_prev_A%02d",a);
        h1_tofraw_prev[a-1] = new TH1D(name,name,20000,-500,1500);
        h1_tofraw_prev[a-1]->SetLineColor(kBlack);
        h1_tofraw_prev[a-1]->SetDirectory(0);

        sprintf(name,"tofraw_curr_A%02d",a);
        h1_tofraw_curr[a-1] = new TH1D(name,name,20000,-500,1500);
        h1_tofraw_curr[a-1]->SetLineColor(kBlue);
        h1_tofraw_curr[a-1]->SetDirectory(0);

        sprintf(name,"tofraw_next_A%02d",a);
        h1_tofraw_next[a-1] = new TH1D(name,name,20000,-500,1500);
        h1_tofraw_next[a-1]->SetLineColor(kRed);
        h1_tofraw_next[a-1]->SetDirectory(0);
    }



    // === =========================================================
    // === input data 
  
    TFile * f = new TFile(Form("../../output/conversion/raw%i.root",run_number),"read");
    TTree * tFC = (TTree*)f->Get("EpicRawTree");



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
                h1_tofraw_prev[tof.anode-1]->Fill(tof.tFC - tof.tHF_prev + tofraw_prev_offset);
                h1_tofraw_curr[tof.anode-1]->Fill(tof.tFC - tof.tHF_curr + tofraw_curr_offset);
                h1_tofraw_next[tof.anode-1]->Fill(tof.tFC - thf_next + tofraw_next_offset);
            }

            // Clear
            pendingFC.clear();

            // newvalue
            thf_prev = epicFC->GetTimePrevHF();
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
            fc.tHF_prev = thf_prev;
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
        h1_tofraw_prev[a-1]->Draw("same");

    }


}
