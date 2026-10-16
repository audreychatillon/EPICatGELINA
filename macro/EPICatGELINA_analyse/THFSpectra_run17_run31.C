#include <iostream>
#include <fstream>
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



////////////////////////////////////////////////////////////////////////////////
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

    double gammapeak[10] = { 
        -887.04      , -8.87141e+02 , -8.87399e+02 , -8.87769e+02 , -8.83413e+02 ,
        -8.81467e+02 , -8.79544e+02 , -8.78904e+02 , -8.88189e+02 , -8.86883e+02 
    };
    double flight_path[10] = {33.8, 27.2, 20.2, 13.6, 6.6, 0., -7.0, -13.6, -20.2, -27.2};

    double actinide_mass[10] = {97.97,  99.23,  103.8,  103.31,  92.67,  95.14,  94.07,  112.69,  107.55,  105.04 };


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
    TH1D * h1_DeltaTHFzoom = new TH1D("DeltaTHFzoom","DeltaTHFzoom",30000,2498000,2501000);
    TH1D * h1_onlyR_THF = new TH1D("THF_onlyReflex","THF_onlyReflex",86400,0,86400);
    TH1D * h1_onlyR_DeltaTHFzoom = new TH1D("DeltaTHFzoom_onlyReflex","DeltaTHFzoom_onlyReflex",20000,870,890);
    TH1D * h1_noR_THF = new TH1D("THF_noReflex","THF_noReflex",86400,0,86400);
    TH1D * h1_noR_DeltaTHF = new TH1D("DeltaTHF_noReflex","DeltaTHF_noReflex",5100,0,5100000);
    TH1D * h1_noR_HFonly_DeltaTHF = new TH1D("DeltaTHFzoom_noReflex_HFonly","DeltaTHFzoom_noReflex_HFonly",5100,0,5100000);
    TH1D * h1_noR_DeltaTHFzoom = new TH1D("DeltaTHFzoom_noReflex","DeltaTHFzoom_noReflex",30000,2498000,2501000);
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
        h1_TofRaw_curr[a] = new TH1D(name,name,20000,-1500,500);
        h1_TofRaw_curr[a]->SetDirectory(0);
        h1_TofRaw_curr[a]->SetLineColor(kBlue);
        sprintf(name,"TofRaw_next_A%02d",a+1);
        h1_TofRaw_next[a] = new TH1D(name,name,20000,-1500,500);
        h1_TofRaw_next[a]->SetDirectory(0);
        h1_TofRaw_next[a]->SetLineColor(kRed);
    }

    // === =========================================================
    // === TofCal histograms
    TH1D * h1_TofCal_next[10];
    for(short a = 0 ; a < 10 ; a++){
        sprintf(name,"TofCal_next_A%02d",a+1);
        h1_TofCal_next[a] = new TH1D(name,name,20000,-100,1900);
        h1_TofCal_next[a]->SetDirectory(0);
        h1_TofCal_next[a]->SetLineColor(kRed);

    }

    // === =========================================================
    // === Ene histograms
    TH1D * h1_e_next[10];
    TH1D * h1_ew_next[10];
    for(short a = 0 ; a < 10 ; a++){
	    h1_e_next[a] = new TH1D(Form("En_next_A%02d",a+1),Form("En_next_A%02d",a+1),200,0,20);
        h1_e_next[a]->GetXaxis()->SetTitle("[MeV] 100keV/bin");
        h1_e_next[a]->SetLineColor(kRed+1);
        h1_e_next[a]->SetDirectory(0);

	    h1_ew_next[a] = new TH1D(Form("En_weighted_next_A%02d",a+1),Form("En_weighted_next_A%02d",a+1),200,0,20);
        h1_ew_next[a]->GetXaxis()->SetTitle("[MeV] 100keV/bin");
        h1_ew_next[a]->SetLineColor(kRed+1);
        h1_ew_next[a]->SetDirectory(0);
    }

    // === =========================================================
    // === TGraph ratio XS
    std::map<std::string,TGraph*> m_phy_gr;
    // JEFF4-0 cross sections
    vector<double> XS_JEFF_nENE;
    vector<double> XS_JEFF_ratio_235U;
    vector<double> XS_JEFF_ratio_238U;
    vector<double> XS_JEFF_ratio_242Pu;
    ifstream XSfile("../../doc/XS_JEFF4-0_235U_238U_242Pu.dat");	
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
    ifstream XSfile2("../../doc/XS_BVIII_235U_238U_242Pu.dat");	
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
        
                // thf_prev = Good / thf_curr = reflexion] : this one should be enough
                if(2500889 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500905){
                    //cout << endl;
                    //cout << "case :  [thf_prev = Good / thf_curr = reflexion] " << endl;
                    //cout << setprecision(25) << " thf_prev = " << thf_prev << ", th_curr = " << thf_curr << endl;
                    thf_curr -= 878.1; 
                    //cout << setprecision(25) << " thf_prev = " << thf_prev << ", th_curr = " << thf_curr << endl;
                    //cout << "new DeltaT = " << thf_curr - thf_prev << endl;
                }

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
                            if(fc.FC_fission){
                                tofraw_curr = fc.FC_time - fc.HF_curr - 2500000.;
                                tofraw_next = fc.FC_time - thf_curr;
                                h1_TofRaw_curr[fc.FC_anode-1]->Fill(tofraw_curr);
                                h1_TofRaw_next[fc.FC_anode-1]->Fill(tofraw_next);
                                double tofcal;
                                double e = TofRaw2Ene(gammapeak[fc.FC_anode-1],FLIGHT_PATH + flight_path[fc.FC_anode-1],tofraw_curr,tofcal);
                                h1_TofCal_next[fc.FC_anode-1]->Fill(tofcal);
                                h1_e_next[fc.FC_anode-1]->Fill(e);
                                h1_ew_next[fc.FC_anode-1]->Fill(e, 1./actinide_mass[fc.FC_anode-1]);
                            }
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
		    fc.FC_alpha = false;
		    fc.FC_fission = false;
		    if(fc.FC_q3>0 && fc.FC_anode>0){ 
                if(cutalpha[fc.FC_anode-1]->IsInside(fc.FC_q1,fc.FC_q2/fc.FC_q3))		fc.FC_alpha = true;
		    	if(cutfission[fc.FC_anode-1]->IsInside(fc.FC_q1,fc.FC_q2/fc.FC_q3))		fc.FC_fission = true; 
		    }
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
			fc.FC_qmax = epicFC->GetQmax(ifc);
			fc.FC_q1 = epicFC->GetQ1(ifc);
			fc.FC_q2 = epicFC->GetQ2(ifc);
			fc.FC_q3 = epicFC->GetQ3(ifc);
		    fc.FC_alpha = false;
		    fc.FC_fission = false;
		    if(fc.FC_q3>0 && fc.FC_anode>0){ 
                if(cutalpha[fc.FC_anode-1]->IsInside(fc.FC_q1,fc.FC_q2/fc.FC_q3))		fc.FC_alpha = true;
		    	if(cutfission[fc.FC_anode-1]->IsInside(fc.FC_q1,fc.FC_q2/fc.FC_q3))		fc.FC_fission = true; 
		    }
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
        
                // thf_prev = Good / thf_curr = reflexion]
                if(2500889 < (thf_curr-thf_prev) && (thf_curr-thf_prev) < 2500905)  thf_curr -= 878.1; 

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
                            if(fc.FC_fission){
                                tofraw_curr = fc.FC_time - fc.HF_curr - 2500000.;
                                tofraw_next = fc.FC_time - thf_curr;
                                h1_TofRaw_curr[fc.FC_anode-1]->Fill(tofraw_curr);
                                h1_TofRaw_next[fc.FC_anode-1]->Fill(tofraw_next);
                                double tofcal;
                                double e = TofRaw2Ene(gammapeak[fc.FC_anode-1],FLIGHT_PATH + flight_path[fc.FC_anode-1],tofraw_curr,tofcal);
                                h1_TofCal_next[fc.FC_anode-1]->Fill(tofcal);
                                h1_e_next[fc.FC_anode-1]->Fill(e);
                                h1_ew_next[fc.FC_anode-1]->Fill(e, 1./actinide_mass[fc.FC_anode-1]);
                            }
                        }
                        pendingFC.clear();
                    }
                    // keep the current HF
                    thf_curr_good = thf_curr;
    		    } // end of [noR_Good]



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

    TCanvas * can3 = new TCanvas("TofCal","TofCal",2000,3000);
    can3->Divide(5,2);
    for(short a = 0 ; a < 10 ; a++){
        can3->cd(a+1); h1_TofCal_next[a]->Draw();
    }

    TCanvas * can4 = new TCanvas("Ene","Ene",0,0,3000,2000);
    can4->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can4->cd(a); gPad->SetLogy();
        h1_e_next[a-1]->Draw();
    }

    TCanvas * can5 = new TCanvas("Ene_weighted","Ene_weighted",0,0,3000,2000);
    can5->Divide(5,2);
    for(short a = 1 ; a <=10 ; a++){
        can5->cd(a); gPad->SetLogy();
        h1_ew_next[a-1]->Draw();
    }

    // histograms per actinide

    TH1D * h1_e_U238_next = new TH1D("En_U238_next","En_U238_next",200,0,20);
    TH1D * h1_e_U235_next = new TH1D("En_U235_next","En_U235_next",200,0,20);
    TH1D * h1_e_Pu242_next = new TH1D("En_Pu242_next","En_Pu242_next",200,0,20);
    h1_e_U238_next->SetDirectory(0);
    h1_e_U235_next->SetDirectory(0);
    h1_e_Pu242_next->SetDirectory(0);
    h1_e_U238_next->SetLineColor(kRed+1);
    h1_e_U235_next->SetLineColor(kRed+1);
    h1_e_Pu242_next->SetLineColor(kRed+1);
  
    h1_e_U238_next->Add(h1_e_next[0]); 
    h1_e_U235_next->Add(h1_e_next[1]); 
    h1_e_U235_next->Add(h1_e_next[2]); 
    h1_e_U235_next->Add(h1_e_next[3]); 
    h1_e_Pu242_next->Add(h1_e_next[4]); 
    h1_e_Pu242_next->Add(h1_e_next[5]); 
    h1_e_Pu242_next->Add(h1_e_next[6]); 
    h1_e_U235_next->Add(h1_e_next[7]); 
    h1_e_U235_next->Add(h1_e_next[8]); 
    h1_e_U235_next->Add(h1_e_next[9]); 

    TH1D * h1_ew_U238_next = new TH1D("En_w_U238_next","En_w_U238_next",200,0,20);
    TH1D * h1_ew_U235_next = new TH1D("En_w_U235_next","En_w_U235_next",200,0,20);
    TH1D * h1_ew_Pu242_next = new TH1D("En_w_Pu242_next","En_w_Pu242_next",200,0,20);
    h1_ew_U238_next->SetDirectory(0);
    h1_ew_U235_next->SetDirectory(0);
    h1_ew_Pu242_next->SetDirectory(0);
    h1_ew_U238_next->SetLineColor(kRed+1);
    h1_ew_U235_next->SetLineColor(kRed+1);
    h1_ew_Pu242_next->SetLineColor(kRed+1);
    for(short bin = 1 ; bin <= h1_ew_U238_next->GetNbinsX(); bin++){
       h1_ew_U238_next->SetBinContent(bin,h1_e_U238_next->GetBinContent(bin)*1./MASS_U238);
       h1_ew_U235_next->SetBinContent(bin,h1_e_U235_next->GetBinContent(bin)*1./MASS_U235);
       h1_ew_Pu242_next->SetBinContent(bin,h1_e_Pu242_next->GetBinContent(bin)*1./MASS_PU242);
    }

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

    can6->cd(1); h1_e_U238_next->Draw();  
    can6->cd(2); h1_e_U235_next->Draw();  
    can6->cd(3); h1_e_Pu242_next->Draw(); 
    can6->cd(4); h1_ew_U238_next->Draw(); 
    can6->cd(5); h1_ew_U235_next->Draw(); 
    can6->cd(6); h1_ew_Pu242_next->Draw();
    can6->cd(7); h1_eratio_U238_next->Draw(); 
    can6->cd(8); h1_eratio_U235_next->Draw(); 
    can6->cd(9); h1_eratio_Pu242_next->Draw();

    can6->cd(7); m_phy_gr["XSratio_JEFF_238U_235U"]->Draw("sameLP");
    can6->cd(9); m_phy_gr["XSratio_JEFF_242Pu_235U"]->Draw("sameLP");
    can6->cd(7); m_phy_gr["XSratio_BVIII_238U_235U"]->Draw("sameLP");
    can6->cd(9); m_phy_gr["XSratio_BVIII_242Pu_235U"]->Draw("sameLP");

    
    TCanvas * can7 = new TCanvas("Ene_Pu242","Ene_Pu242",0,0,3000,2000);
    can7->Divide(2,2);

    // histograms per anode of Pu242
    TH1D * h1_eratio_Pu242_per_anode_next[3];
    for(short a = 1 ; a <= 3 ; a++){
    	h1_eratio_Pu242_per_anode_next[a-1] = new TH1D(Form("En_ratio_U235_Pu242_A%02d_next",a+4),Form("En_ratio_U235_Pu242_A%02d_next",a+4),200,0,20);
    	h1_eratio_Pu242_per_anode_next[a-1]->SetDirectory(0);
    	h1_eratio_Pu242_per_anode_next[a-1]->SetLineColor(kRed+1);
    	h1_eratio_Pu242_per_anode_next[a-1]->Add(h1_ew_next[a+3]);
    	h1_eratio_Pu242_per_anode_next[a-1]->Divide(h1_ew_U235_next);
        can7->cd(a); gPad->SetGridx(); gPad->SetGridy(); 
        h1_eratio_Pu242_per_anode_next[a-1]->Draw(); 
        m_phy_gr["XSratio_JEFF_242Pu_235U"]->Draw("sameLP");
    }
    can7->cd(4); h1_eratio_Pu242_next->Draw();
    m_phy_gr["XSratio_JEFF_242Pu_235U"]->Draw("sameLP");


}
