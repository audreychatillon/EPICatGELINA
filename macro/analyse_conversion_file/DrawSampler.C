#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <unordered_map>

#include <Riostream.h>
#include <TROOT.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TTree.h>

#include "../ClassDef/EpicRawTree.h"

// legend: T_CFD: blue thick line from plugin
//         T_max: blue thin line  from plugin
//         newT_CFD: red thick line
//         new_Tmax: red thin line

#define PLOT_MEAN_SIGNALS    1
#define PLOT_REJECTED_SIGNAL 0
#define DISPLAY_STYLE        0 // 0: no online-display, 1: press enter, 2: latence display, 3: double-click

#if DISPLAY_STYLE == 2
#define LATENCE_DISPLAY_MS 1500
#endif

#if PLOT_REJECTED_SIGNAL
#define can_rej_Xstart 0
#define can_rej_Ystart 0
#define can_rej_Xwidth 500
#define can_rej_Ywidth 500
#endif

#define faster_sample_size_ns 2

//================================================================================================
//================================================================================================
//================================================================================================
double calculateCFD(vector<double> Signal, vector<double> &samplesX,
                    vector<double> &signalCFD, const double fract,
                    const int delay, const double thres, // CFD constants
                    double &Qmax,
                    double &tQmax, // Qmax of the Signal vector and its time
                    bool &bTriggered, bool &bThreshold) {

  const int nsamples_delay = delay / faster_sample_size_ns;
  double CFDmax = -1000.;
  int index_Qmax = -1;
  int index_CFDmax = -1;
  int index_thres = -1;
  int index_first = -1;
  double Tcfd = -10000.;
  if (Signal.size() == 0)
    return Tcfd;

  for (int i = 0; i < (int)Signal.size(); ++i) {

    // === fill time scale
    samplesX[i] = i * (double)faster_sample_size_ns + 1.;

    // === get Qmax
    double Qi = Signal[i];
    if (Qi > Qmax) {
      Qmax = Qi;
      index_Qmax = i;
      tQmax = samplesX[i];
    }

    // === Convert to CFD signal
    if (i >= nsamples_delay) {
      //    "The input signal is attenuated and inversely added to the delayed
      //    input signal" therefore: CFDsignal = delayed_signal -
      //    fraction*signal for these positive PA signals: signals should first
      //    be <0 before to cross 0
      signalCFD[i] = Signal[i - nsamples_delay] - fract * Qi;
      if (signalCFD[i] > CFDmax) {
        CFDmax = signalCFD[i];
        index_CFDmax = i;
      }
    } // end of CFD conversion for bin i
    // else{
    //   signalCFD[i] = 0 ;
    // }
  } // end of loop over Signal

  // === Starting from the maximum value of the CFD
  //     find the first negative value
  for (int i = index_CFDmax; i > 0; i--) {
    if (signalCFD[i] < 0) {
      index_first = i;
      bTriggered = true;
      break;
    }
  }
  // === Starting from the first negative value before max
  //     find the first bin above threshold
  for (int i = index_first; i <= index_CFDmax; i++)
    if (signalCFD[i] > thres) {
      index_thres = i;
      bThreshold = true;
      break;
    }

  // === Check that there is at least 3 consecutive points above threshold
  if (bThreshold) {
    int npts_above_thres = 0;
    for (int i = index_thres; i <= index_CFDmax + 2; i++) {
      if (signalCFD[i] >= thres)
        npts_above_thres++;
      else
        npts_above_thres = 0;
      if (npts_above_thres >= 3)
        break;
    }
    if (npts_above_thres < 3)
      bThreshold = false;
  }

  /*
    // === Starting from the first negative value before max
    //     find the first bin above threshold
    //     check that there are at least 3 points above threshold
    // === check that there is at least 3 consecutive points above threshold
    for (int i = index_first ; i <= index_CFDmax; i++){
      if (signalCFD[i] > thres){
          index_thres = i;
          if (index_thres < (int)(signalCFD.size()-3) && index_thres <=
    (index_CFDmax-2)) bThreshold = true; break;
      }
    }
    if(bThreshold){
      for(int i = index_CFDmax; i >= index_CFDmax-2 ; i--)
          if(signalCFD[i]<thres) bThreshold = false;
    }
  */

#if PLOT_REJECTED_SIGNAL
  if (!bThreshold) {
    cout << endl;
    cout << "--- bThreshold = false " << thres << endl;
    cout << "    maxCFD = " << signalCFD[index_CFDmax]
         << " at bin = " << index_CFDmax
         << " t_cfd_max = " << samplesX[index_CFDmax] << endl;
    cout << "    first negative value  at bin =  " << index_first
         << " t_first = " << samplesX[index_first] << endl;
    for (int i = index_first; i <= index_CFDmax + 2; i++)
      cout << "          bin i = " << i << ", signalCFD[" << i
           << "] = " << signalCFD[i] << " , signal[" << i << "] " << Signal[i]
           << endl;
  }
#endif

  // === Calculate Tcfd = signalCFD cross 0
  //     quadratic function cfd = a.t^2 + b.t + c
  if (bThreshold && bTriggered) {
    double x1 = (double)samplesX[index_first];
    double x2 = (double)samplesX[index_first + 1];
    double x3 = (double)samplesX[index_first + 2];
    double y1 = (double)signalCFD[index_first];
    double y2 = (double)signalCFD[index_first + 1];
    double y3 = (double)signalCFD[index_first + 2];
    // Tcfd = x1 - (double)faster_sample_size_ns * y1 / (y2 - y1) ;

    // cfd => y and t => x
    // we have 3 equations with 3 unknowns :  yi = a * xi*xi + b * xi + c
    //    from (3):       c = y3 - a * x3*x3 - b * x3
    // => c in (1): y1 - y3 = a * (x1*x1 - x3*x3) + b * (x1 - x3)
    //    c in (2): y2 - y3 = a * (x2*x2 - x3*x3) + b * (x2 - x3)
    // => from (1):       a = (y1 - y3) / (x1*x1 - x3*x3) - b / (x1 + x3)
    // => a in (2): y2 - y3 = (y1 - y3)*(x2*x2 - x3*x3) / (x1*x1 - x3*x3) - b *
    // (x2*x2 - x3*x3) / (x1 + x3) + b * (x2 - x3)
    //              y2 - y3 = (y1 - y3)*(x2*x2 - x3*x3) / (x1*x1 - x3*x3) + b *
    //              (x2 - x3 - (x2*x2 - x3*x3) / (x1 + x3))
    //                    b = ( y2-y3 - (x2*x2-x3*x3)*(y1-y3)/(x1*x1-x3*x3)) ) /
    //                    ( x2-x3 - (x2*x2-x3*x3)/(x1+x3) )
    double b = ((y2 - y3) -
                (((x2 * x2 - x3 * x3) * (y1 - y3)) / (x1 * x1 - x3 * x3))) /
               ((x2 - x3) - ((x2 * x2 - x3 * x3) / (x1 + x3)));
    double a = ((y1 - y2) / (x1 * x1 - x2 * x2)) - (b / (x1 + x2));
    double c = y1 - a * x1 * x1 - b * x1;
    double d = b * b - 4 * a * c;

    if (d < 0) {
      std::cout << " no CFD time found !!!" << std::endl;
      std::cout << " quadratic fit gives : " << std::endl;
      std::cout << " a = " << a << ", b = " << b << ", c = " << c
                << " => DELTA = " << d << std::endl;
    } else if (d == 0) {
      Tcfd = -b / (2 * a);
    } else {
      if (std::abs(a) > 1.e-12) {
        double Tcfd_1 = (-b - sqrt(d)) / (2. * a);
        double Tcfd_2 = (-b + sqrt(d)) / (2. * a);
        if (x1 < Tcfd_1 && Tcfd_1 < x2) {
          Tcfd = Tcfd_1;
          if (x1 < Tcfd_2 && Tcfd_2 < x3)
            std::cout << "2 solutions for Tcfd : " << Tcfd_1 << ", and "
                      << Tcfd_2 << std::endl;
        } else if (x1 < Tcfd_2 && Tcfd_2 < x3) {
          Tcfd = Tcfd_2;
        } else {
          // this was always happening when a~0
          std::cout << "----no Tcfd found : min = " << x1 << ", max = " << x3
                    << " ns" << std::endl;
          std::cout << "    time_x1 = " << x1 << " , cfd_y1 = " << y1
                    << std::endl;
          std::cout << "    time_x2 = " << x2 << " , cfd_y2 = " << y2
                    << std::endl;
          std::cout << "    time_x3 = " << x3 << " , cfd_y3 = " << y3
                    << std::endl;
          std::cout << "       quadratic fit gives : a = " << a << ", b = " << b
                    << ", c = " << c << " => DELTA = " << d << std::endl;
          std::cout << "       Tcfd_1 = " << Tcfd_1 << std::endl;
          std::cout << "       Tcfd_2 = " << Tcfd_2 << std::endl;
        }
      } // end of else abs(a) > 1.e-12
      else { // linear extrapolation
        Tcfd = x1 - (double)faster_sample_size_ns * y1 / (y2 - y1);
      }
    }
  }
  return Tcfd;
}

//================================================================================================
//================================================================================================
//================================================================================================
double integrateSignal(double t_start, double t_stop, vector<double> signal) {

  // int    i_start = TMath::Max(0,(int)floor(t_start / faster_sample_size_ns))
  // ; int    i_stop  = TMath::Min((int)signal.size()-1,(int)floor(t_stop  /
  // faster_sample_size_ns)) ;
  int i_start =
      TMath::Max(0, static_cast<int>(t_start / faster_sample_size_ns));
  int i_stop = TMath::Min((int)signal.size() - 1,
                          static_cast<int>(t_stop / faster_sample_size_ns));

  if (t_start >= t_stop || i_start >= i_stop - 1)
    return -1;

  double frac_start = (((i_start + 1) * faster_sample_size_ns) - t_start) /
                      faster_sample_size_ns;
  double frac_stop =
      (t_stop - (i_stop * faster_sample_size_ns)) / faster_sample_size_ns;
  double q = 0;

  q = frac_start * signal[i_start];
  q += frac_stop * signal[i_stop];
  for (int i = i_start + 1; i <= i_stop - 1; i++)
    q += signal[i];

  q = q * faster_sample_size_ns;

  return q;
}

//================================================================================================
//================================================================================================
//================================================================================================
void run(int data_set) {

  // ==========================================================================================
  // === VARIABLES

  char name[100];

  int nA;
  int Pmbar;
  int Q1max;
  int HV;
#if PLOT_MEAN_SIGNALS
  TFile *fcut;
  TCutG *cutg[8];
  TH2F *h2_sig[8];
#endif

  short det;
  short anode;
  double qmax;
  double q1;
  double q2;
  double q3;
  double tof_raw;
  double t_cfd;
  double t_qmax;

  int index;

  int i_qmax;
  double q_qmax;

  // ===========================================================================================
  // === INPUT DATA
  vector<double> *fSampler_Signal = nullptr;

  TChain *ch = new TChain("EpicRawTree");

  switch (data_set) {
  case 1: // PA = V4b, P=1030mbar, HV=550V
    nA = 9;
    HV = 550;
    Pmbar = 1030;
    Q1max = 300000;
    ch->Add(Form("../../output/conversion/Test1_V4b_1030mbar_550V.root"));
    ch->Add(Form("../../output/conversion/Test2_V4b_1030mbar_550V.root"));
    ch->Add(Form("../../output/conversion/Test5_V4b_1030mbar_550V.root"));
    ch->Add(Form("../../output/conversion/Test6_V4b_1030mbar_550V.root"));
    break;
  case 2: // PA = V4b, P=1190mbar, HV = 610V
    nA = 9;
    HV = 610;
    Pmbar = 1190;
    Q1max = 300000;
    ch->Add(Form("../../output/conversion/Test3_V4b_1030mbar_610V.root"));
    break;
  case 3: // PA = V4b, P=1190mbar, HV = 500V
    nA = 9;
    HV = 500;
    Pmbar = 1190;
    Q1max = 300000;
    ch->Add(Form("../../output/conversion/Test4_V4b_1030mbar_500V.root"));
    break;
  case 4: // PA = V8.2, P=1380mbar, HV = 650V
    nA = 4;
    HV = 650;
    Pmbar = 1380;
    Q1max = 300000;
#if PLOT_MEAN_SIGNALS
    fcut = new TFile(Form("TCutG_DISCRI_DataSet%i.root", data_set), "read");
    fcut->ls();
    for (int i = 1; i <= 3; i++) {
      cutg[i - 1] = (TCutG *)fcut->Get(Form("A10_DISCRI_ALPHA%i", i));
      h2_sig[i - 1] =
          new TH2F(Form("A10_cutALPHA%i", i), Form("A10_cutALPHA%i", i), 70, 0,
                   140, 1600, -400, 1200);
      h2_sig[i - 1]->SetDirectory(0);
    }
    for (int i = 1; i <= 5; i++) {
      cutg[i + 2] = (TCutG *)fcut->Get(Form("A10_DISCRI_FF%i", i));
      h2_sig[i + 2] = new TH2F(Form("A10_cutFF%i", i), Form("A10_cutFF%i", i),
                               70, 0, 140, 1300, -1000, 12000);
      h2_sig[i + 2]->SetDirectory(0);
    }
    fcut->Close();
    for (int i = 0; i < 8; i++) {
      cutg[i]->ls();
      h2_sig[i]->ls();
    }
#endif
    ch->Add(Form("../../output/conversion/Test9_V8-2_1380mbar_650V.root"));
    ch->Add(Form("../../output/conversion/Test13_V8-2_1380mbar_650V.root"));
    break;
  case 5: // PA = V8.2, P=1040mbar, HV = 550V
    nA = 4;
    HV = 550;
    Pmbar = 1040;
    Q1max = 300000;
    ch->Add(Form("../../output/conversion/Test10_V8-2_1040mbar_550V.root"));
    ch->Add(Form("../../output/conversion/Test12_V8-2_1040mbar_550V.root"));
    break;
  case 6: // PA = V8.2, P=1210mbar, HV = 680V
    nA = 4;
    HV = 680;
    Pmbar = 1210;
    Q1max = 500000;
    ch->Add(Form("../../output/conversion/Test14_V8-2_1210mbar_680V.root"));
    break;
  case 7: // PA = V8.2, P=1370mbar, HV = 720V
    nA = 4;
    HV = 720;
    Pmbar = 1370;
    Q1max = 500000;
    ch->Add(Form("../../output/conversion/Test15_V8-2_1370mbar_720V.root"));
    ch->Add(Form("../../output/conversion/Test16_V8-2_1370mbar_720V.root"));
    break;
  case 8: // PA = V8.2, P=1075mbar, HV = 720V
    nA = 4;
    HV = 720;
    Pmbar = 1075;
    Q1max = 500000;
    ch->Add(
        Form("../../output/conversion/Test17_V8-2_%imbar_%iV.root", Pmbar, HV));
    ch->Add(
        Form("../../output/conversion/Test18_V8-2_%imbar_%iV.root", Pmbar, HV));
    break;
  case 9: // PA = V8.2, P=1065mbar, HV = 750V
    nA = 4;
    HV = 750;
    Pmbar = 1065;
    Q1max = 500000;
    ch->Add(
        Form("../../output/conversion/Test19_V8-2_%imbar_%iV.root", Pmbar, HV));
    break;
  default:
    nA = 0;
    HV = 0;
    Pmbar = 0;
    Q1max = 0;
    break;
  }

  unordered_map<int, int> mapindex;
  unordered_map<int, int> mapanode;
  if (nA == 9) {
    mapindex = {{1, 0}, {2, 1}, {3, 2},  {4, 3}, {6, 4},
                {8, 5}, {9, 6}, {10, 7}, {11, 8}};
  } else if (nA == 4) {
    mapindex = {{3, 0}, {4, 1}, {9, 2}, {10, 3}};
  }
  for (const auto &p : mapindex)
    mapanode[p.second] = p.first;
  int ncol = ceil(nA * 0.5);
  int ncolnew = ceil((nA + 1) * 0.5);

  EpicRawTree raw(ch);

  // CFD parameters
  double CFDthres = 10.;
  double CFDfract = 1. / 3.;
  double CFDdelay = 10.;

  // discri parameters
  vector<double> Q1gate_start(nA, -16);
  vector<double> Q1gate_stop(nA, 30);
  vector<double> Q2gate_start(nA, -16);
  vector<double> Q2gate_stop(nA, 2);
  vector<double> Q3gate_start(nA, 8);
  vector<double> Q3gate_stop(nA, 30);
  vector<double> Q4gate_start(nA, 40);
  vector<double> Q4gate_stop(nA, 80);

  // ===========================================================================================
  // === HISTOGRAMS

  TH1I *h1_multFC = new TH1I("MultFC", "MultFC", 10, -0.5, 9.5);

  TH2F *h2_Q1_ifQmax;
  vector<TH1F *> h1_Q1_ifQmax(nA);
  vector<TH2F *> h2_DT_vs_Q1_ifQmax(nA);
  vector<TH2F *> h2_DT_vs_Qmax_ifQmax(nA);
  vector<TH2F *> h2_newDT_vs_Qmax_ifQmax(nA);
  vector<TH2F *> h2_newDT_vs_Q1_ifQmax(nA);
  vector<TH2F *> h2_Discri(nA);
  vector<TH2F *> h2_newDiscri(nA);
  vector<TH2F *> h2_newQ4vQ1(nA);
  vector<TH2F *> h2_newQ4Q2vQ1(nA);
  vector<TH2F *> h2_newQ4Q3vQ1(nA);
  vector<TH2F *> h2_MvQ1_ifQmax(nA);

#if DISPLAY_STYLE
  TCanvas *can_samples;
  vector<TGraph *> gr_Signal(nA);
  vector<TGraph *> gr_CFD(nA);
  vector<TLine *> l_Tcfd(nA);
  vector<TLine *> l_Tmax(nA);
  vector<TLine *> l_newTcfd(nA);
  vector<TLine *> l_newTmax(nA);
  vector<TLine *> l_Q1start(nA);
  vector<TLine *> l_Q1stop(nA);
  vector<TLine *> l_Q2stop(nA);
  vector<TLine *> l_Q3start(nA);
#endif

#if PLOT_REJECTED_SIGNAL
  TCanvas *can_rejected_CFD = new TCanvas(
      "rejected_signal_from_CFD", "rejected_signal_from_CFD", can_rej_Xstart,
      can_rej_Ystart, can_rej_Xwidth, can_rej_Ywidth);
  can_rejected_CFD->cd();
  gPad->SetGridy();
  TGraph *gr_rejected_Signal = nullptr;
  TGraph *gr_rejected_CFD = nullptr;
#endif

  sprintf(name, "P%imbar_Q1_vs_A_ifQmax_%i", Pmbar, data_set);
  h2_Q1_ifQmax = new TH2F(name, name, 13, -0.5, 12.5, Q1max / 200, 0, Q1max);
  h2_Q1_ifQmax->SetDirectory(0);

#if DISPLAY_STYLE
  sprintf(name, "P%imbar_SAMPLES_ifQmax_%i", Pmbar, data_set);
  can_samples = new TCanvas(name, name, 0, 0, 2000, 1000);
  can_samples->Divide(ncolnew, 2);
  can_samples->cd(nA + 1);
  gPad->SetLogz();
#endif

  for (int a = 0; a < nA; a++) {

    sprintf(name, "P%imbar_Anode%i_Q1_ifQmax_%i", Pmbar, mapanode[a], data_set);
    h1_Q1_ifQmax[a] = new TH1F(name, name, Q1max / 100, 0, Q1max);
    h1_Q1_ifQmax[a]->SetLineColor(kBlue);
    h1_Q1_ifQmax[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_DT_vs_Q1_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_DT_vs_Q1_ifQmax[a] =
        new TH2F(name, name, Q1max / 200, 0, Q1max, 5000, -100, 400);
    h2_DT_vs_Q1_ifQmax[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_DT_vs_Qmax_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_DT_vs_Qmax_ifQmax[a] =
        new TH2F(name, name, 500, 0, 5000, 5000, -100, 400);
    h2_DT_vs_Qmax_ifQmax[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_newDT_vs_Qmax_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_newDT_vs_Qmax_ifQmax[a] =
        new TH2F(name, name, 500, 0, 5000, 5000, -100, 400);
    h2_newDT_vs_Qmax_ifQmax[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_newDT_vs_Q1_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_newDT_vs_Q1_ifQmax[a] =
        new TH2F(name, name, Q1max / 200, 0, Q1max, 5000, -100, 400);
    h2_newDT_vs_Q1_ifQmax[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_Q2Q3_vs_Q1_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_Discri[a] = new TH2F(name, name, Q1max / 200, 0, Q1max, 5000, 0, 5);
    h2_Discri[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_new_Q2Q3_vs_Q1_ifQmax_%i", Pmbar,
            mapanode[a], data_set);
    h2_newDiscri[a] = new TH2F(name, name, Q1max / 200, 0, Q1max, 5000, 0, 5);
    h2_newDiscri[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_newQ4_vs_newQ1_ifQmax_%i", Pmbar,
            mapanode[a], data_set);
    h2_newQ4vQ1[a] =
        new TH2F(name, name, Q1max / 200, 0, Q1max, 1200, -50000, 10000);
    h2_newQ4vQ1[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_Q4Q2_vs_Q1_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_newQ4Q2vQ1[a] = new TH2F(name, name, Q1max / 200, 0, Q1max, 1200, -5, 1);
    h2_newQ4Q2vQ1[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_Q4Q3_vs_Q1_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_newQ4Q3vQ1[a] = new TH2F(name, name, Q1max / 200, 0, Q1max, 1200, -5, 1);
    h2_newQ4Q3vQ1[a]->SetDirectory(0);

    sprintf(name, "P%imbar_Anode%i_Mult_vs_Q1_ifQmax_%i", Pmbar, mapanode[a],
            data_set);
    h2_MvQ1_ifQmax[a] =
        new TH2F(name, name, Q1max / 100, 0, Q1max, 15, -0.4, 14.5);
    h2_MvQ1_ifQmax[a]->SetDirectory(0);

#if DISPLAY_STYLE
    can_samples->cd(a + 1);
    gPad->SetGridy();
    gr_Signal[a] = nullptr;
    gr_CFD[a] = nullptr;
    l_Tcfd[a] = nullptr;
    l_Q1start[a] = nullptr;
    l_Q1stop[a] = nullptr;
    l_Tmax[a] = nullptr;
    l_newTcfd[a] = nullptr;
    l_newTmax[a] = nullptr;
#endif
  } // end of for(a)

  // ===========================================================================================
  // === LOOP
  double TimeOffset = 0;
  Long64_t nentries = (Long64_t)ch->GetEntries();
  cout << "nentries = " << nentries << endl;
  for (Long64_t entry = 0; entry < nentries; entry++) {

    ch->GetEntry(entry);
#if DISPLAY_STYLE == 2
    if ((entry % 50) == 0)
      cout << "\r === Entry = " << entry << " / " << nentries
           << " === " << flush;
#elif DISPLAY_STYLE != 2
    if ((entry % 1000000) == 0)
      cout << "\r === Entry = " << entry << " / " << nentries
           << " === " << flush;
#endif

    int fFC_size = (int)raw.fFC_AnodeNbr.size();
    if (fFC_size <= 0)
      continue;
    if (raw.fQmax_Index < 0)
      continue; // fHF data only
    h1_multFC->Fill(raw.fFC_AnodeNbr.size());

    i_qmax = raw.fQmax_Index;
    anode = raw.fFC_AnodeNbr[i_qmax];
    qmax = raw.fFC_Qmax[i_qmax];
    q1 = raw.fFC_Q1[i_qmax];
    q2 = raw.fFC_Q2[i_qmax];
    q3 = raw.fFC_Q3[i_qmax];
    t_cfd = raw.fFC_TimeCfd[i_qmax];
    t_qmax = raw.fFC_TimeQmax[i_qmax];
    fSampler_Signal = &raw.fQmax_Sampler;
    index = mapindex[anode];

    if (q3 > 0) {
      h2_Q1_ifQmax->Fill(anode, q1);
      h1_Q1_ifQmax[index]->Fill(q1);
      h2_DT_vs_Q1_ifQmax[index]->Fill(q1, t_qmax - t_cfd);
      h2_DT_vs_Qmax_ifQmax[index]->Fill(qmax, t_qmax - t_cfd);
      h2_Discri[index]->Fill(q1, q2 / q3);
      h2_MvQ1_ifQmax[index]->Fill(q1, raw.fFC_AnodeNbr.size());
#if PLOT_MEAN_SIGNALS
      if (anode == 10) {
        for (int c = 0; c < 8; c++) {
          if (cutg[c]->IsInside(q1, q2 / q3)) {
            for (int s = 0; s < (int)fSampler_Signal->size(); s++) {
              h2_sig[c]->Fill(s * 2 + 1, fSampler_Signal->at(s));
            }
            break;
          }
        }
      }
#endif
    }

    // calculate CFD
    int faster_signal_nsamples = fSampler_Signal->size();
    if (faster_signal_nsamples == 0) {
      cout << endl << "faster_signal_nsamples==0 but q1_max = " << q1 << endl;
    } else {
      vector<double> time(faster_signal_nsamples, 0);
      vector<double> vCFD(faster_signal_nsamples, 0);
      double newQmax = -10000.;
      double newTmax = -20000.;
      double newQ1 = -1;
      double newQ2 = -1;
      double newQ3 = -1;
      double newQ4 = -10000;
      bool isTrig = false;
      bool isThrs = false;
      double newTcfd =
          calculateCFD(*fSampler_Signal, time, vCFD, CFDfract, CFDdelay,
                       CFDthres, newQmax, newTmax, isTrig, isThrs);
      if ((!isTrig || !isThrs) && (faster_signal_nsamples > 0)) {
#if PLOT_REJECTED_SIGNAL
        if (gr_rejected_Signal)
          delete gr_rejected_Signal;
        if (gr_rejected_CFD)
          delete gr_rejected_CFD;
        gr_rejected_Signal = new TGraph(faster_signal_nsamples, time.data(),
                                        fSampler_Signal->data());
        gr_rejected_Signal->SetName(Form("REJECTED_FASTER_SIGNAL_A%i", anode));
        gr_rejected_Signal->SetTitle(Form("REJECTED_FASTER_SIGNAL_A%i", anode));
        gr_rejected_Signal->SetMarkerStyle(20);
        gr_rejected_Signal->SetMarkerSize(0.5);
        gr_rejected_Signal->SetMarkerColor(kBlue);
        gr_rejected_Signal->SetLineColor(kBlue);
        gr_rejected_Signal->SetLineWidth(2);
        gr_rejected_CFD =
            new TGraph(faster_signal_nsamples, time.data(), vCFD.data());
        gr_rejected_CFD->SetName(
            Form("REJECTED_CFD_SIGNAL_A%i_P%i", anode, Pmbar));
        gr_rejected_CFD->SetTitle(
            Form("REJECTED_CFD_SIGNAL_A%i_P%i", anode, Pmbar));
        gr_rejected_CFD->SetMarkerStyle(20);
        gr_rejected_CFD->SetMarkerSize(0.15);
        gr_rejected_CFD->SetMarkerColor(kBlack);
        gr_rejected_CFD->SetLineColor(kBlack);
        can_rejected_CFD->cd();
        if (newQmax > 200)
          gr_rejected_Signal->GetHistogram()->GetYaxis()->SetRangeUser(
              -500, 1.1 * newQmax);
        else
          gr_rejected_Signal->GetHistogram()->GetYaxis()->SetRangeUser(-100,
                                                                       200);
        gr_rejected_Signal->Draw("ALP");
        gr_rejected_CFD->Draw("sameLP");
        gPad->Update();
        if (newQmax > 200) {
          cout << "\r Press enter for the next event..." << endl;
          cin.get(); // wait up to enter
        } else
          this_thread::sleep_for(chrono::milliseconds(1));
#endif
      } else {
        if (newTcfd < 0)
          cout << "valid event but new Tcfd = " << newTcfd << endl;
        h2_newDT_vs_Qmax_ifQmax[index]->Fill(newQmax, newTmax - newTcfd);
        newQ1 = integrateSignal(newTcfd + Q1gate_start[index],
                                newTcfd + Q1gate_stop[index], *fSampler_Signal);
        newQ2 = integrateSignal(newTcfd + Q2gate_start[index],
                                newTcfd + Q2gate_stop[index], *fSampler_Signal);
        newQ3 = integrateSignal(newTcfd + Q3gate_start[index],
                                newTcfd + Q3gate_stop[index], *fSampler_Signal);
        newQ4 = integrateSignal(newTcfd + Q4gate_start[index],
                                newTcfd + Q4gate_stop[index], *fSampler_Signal);
        h2_newDT_vs_Q1_ifQmax[index]->Fill(newQ1, newTmax - newTcfd);
        if (newQ2 > 0) {
          h2_newQ4Q2vQ1[index]->Fill(newQ1, newQ4 / newQ2);
        }
        if (newQ3 > 0) {
          h2_newDiscri[index]->Fill(newQ1, newQ2 / newQ3);
          h2_newQ4Q3vQ1[index]->Fill(newQ1, newQ4 / newQ3);
        }
        h2_newQ4vQ1[index]->Fill(newQ1, newQ4);
      }
#if DISPLAY_STYLE
      // TGraph : CFD

      if (gr_CFD[index])
        delete gr_CFD[index];
      gr_CFD[index] =
          new TGraph(faster_signal_nsamples, time.data(), vCFD.data());
      gr_CFD[index]->SetName(Form("CFD_SIGNAL_A%i_P%i", anode, Pmbar));
      gr_CFD[index]->SetTitle(
          Form("CFD: A%i P%imbar Q1=%.0f", anode, Pmbar, raw.fFC_Q1[i_qmax]));
      gr_CFD[index]->SetMarkerStyle(20);
      gr_CFD[index]->SetMarkerSize(0.15);
      gr_CFD[index]->SetMarkerColor(kBlack);
      gr_CFD[index]->SetLineColor(kBlack);

      // TGraph : FASTER SIGNAL

      if (gr_Signal[index])
        delete gr_Signal[index];
      gr_Signal[index] = new TGraph(faster_signal_nsamples, time.data(),
                                    fSampler_Signal->data());
      gr_Signal[index]->SetName(Form("FASTER_SIGNAL_A%i_P%i", anode, Pmbar));
      gr_Signal[index]->SetTitle(
          Form("A%i P%imbar Q1=%.0f", anode, Pmbar, raw.fFC_Q1[i_qmax]));
      gr_Signal[index]->SetMarkerStyle(20);
      gr_Signal[index]->SetMarkerSize(0.5);
      gr_Signal[index]->SetMarkerColor(kBlue);
      gr_Signal[index]->SetLineColor(kBlue);
      gr_Signal[index]->SetLineWidth(2);

      // DRAW

      auto ymax = gr_Signal[index]->GetY()[TMath::LocMax(
          gr_Signal[index]->GetN(), gr_Signal[index]->GetY())];
      auto ymin = gr_CFD[index]->GetY()[TMath::LocMin(gr_CFD[index]->GetN(),
                                                      gr_CFD[index]->GetY())];

      t_cfd = raw.fFC_TimeCfd[i_qmax];
      if (l_Tcfd[index])
        delete l_Tcfd[index];
      l_Tcfd[index] = new TLine(t_cfd, 1.3 * ymin, t_cfd, 1.1 * ymax);
      l_Tcfd[index]->SetLineWidth(4);
      l_Tcfd[index]->SetLineStyle(2);
      l_Tcfd[index]->SetLineColor(kBlue);

      t_qmax = raw.fFC_TimeQmax[i_qmax];
      if (l_Tmax[index])
        delete l_Tmax[index];
      l_Tmax[index] = new TLine(t_qmax, ymin, t_qmax, ymax);
      l_Tmax[index]->SetLineWidth(2);
      l_Tmax[index]->SetLineStyle(2);
      l_Tmax[index]->SetLineColor(kBlue);

      if (l_newTcfd[index])
        delete l_newTcfd[index];
      l_newTcfd[index] = new TLine(newTcfd, 1.3 * ymin, newTcfd, 1.1 * ymax);
      l_newTcfd[index]->SetLineWidth(4);
      l_newTcfd[index]->SetLineColor(kRed);

      if (l_newTmax[index])
        delete l_newTmax[index];
      l_newTmax[index] = new TLine(newTmax, ymin, newTmax, ymax);
      l_newTmax[index]->SetLineWidth(2);
      l_newTmax[index]->SetLineColor(kRed);

      can_samples->cd(index + 1);
      gr_Signal[index]->GetHistogram()->GetYaxis()->SetRangeUser(1.3 * ymin,
                                                                 1.1 * ymax);
      gr_Signal[index]->Draw("ALP");
      gr_CFD[index]->Draw("sameLP");
      if (isTrig && isThrs) {
        if (l_Q1start[index])
          delete l_Q1start[index];
        l_Q1start[index] = new TLine(newTcfd + Q1gate_start[index], ymin,
                                     newTcfd + Q1gate_start[index], ymax);
        l_Q1start[index]->SetLineWidth(2);
        l_Q1start[index]->SetLineStyle(1);
        l_Q1start[index]->SetLineColor(8);
        if (l_Q1stop[index])
          delete l_Q1stop[index];
        l_Q1stop[index] = new TLine(newTcfd + Q1gate_stop[index], ymin,
                                    newTcfd + Q1gate_stop[index], ymax);
        l_Q1stop[index]->SetLineWidth(2);
        l_Q1stop[index]->SetLineStyle(2);
        l_Q1stop[index]->SetLineColor(8);
        if (l_Q2stop[index])
          delete l_Q2stop[index];
        l_Q2stop[index] = new TLine(newTcfd + Q2gate_stop[index], ymin,
                                    newTcfd + Q2gate_stop[index], ymax);
        l_Q2stop[index]->SetLineWidth(2);
        l_Q2stop[index]->SetLineStyle(3);
        l_Q2stop[index]->SetLineColor(8);
        if (l_Q3start[index])
          delete l_Q3start[index];
        l_Q3start[index] = new TLine(newTcfd + Q2gate_stop[index], ymin,
                                     newTcfd + Q2gate_stop[index], ymax);
        l_Q3start[index]->SetLineWidth(2);
        l_Q3start[index]->SetLineStyle(4);
        l_Q3start[index]->SetLineColor(8);
        l_Q1start[index]->Draw("same");
        l_Q1stop[index]->Draw("same");
        l_Q2stop[index]->Draw("same");
        l_Q3start[index]->Draw("same");
        l_newTcfd[index]->Draw("same");
        l_newTmax[index]->Draw("same");
      } else {
        cout << "newCFD OFF:  A " << anode << " : bTriggered " << isTrig
             << ", bThreshold " << isThrs << endl;
      }
      l_Tcfd[index]->Draw("same");
      l_Tmax[index]->Draw("same");
      gPad->Update();

      can_samples->cd(nA + 1);
      h2_Q1_ifQmax->Draw("colz");
      gPad->Update();
#if DISPLAY_STYLE == 1
      // === enter for the next event
      cout << "\r Press enter for the next event..." << fflush;
      cin.get(); // wait up to enter
#elif DISPLAY_STYLE == 2
      // === wait the defined time
      this_thread::sleep_for(chrono::milliseconds(LATENCE_DISPLAY_MS));
#else
      // === double-click for next event
      gPad->GetCanvas()->WaitPrimitive();
#endif
#endif
    } // end of else (faster_signal_nsamples > 0 )
  } // end of loop over the entries
  cout << endl;

  // ===========================================================================================
  // === DRAW

  TCanvas *can1 =
      new TCanvas("multFC_per_evt", "multFC_per_evt", 0, 0, 1000, 1000);
  can1->cd();
  gPad->SetLogy();
  h1_multFC->Draw();

  sprintf(name, "P%imbar_Q1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can2 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can2->Divide(ncolnew, 2);

  sprintf(name, "P%imbar_DTvQ1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can3 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can3->Divide(ncol, 2);

  sprintf(name, "P%imbar_DTvQmax_ifQmax_%i", Pmbar, data_set);
  TCanvas *can4 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can4->Divide(ncol, 2);

  sprintf(name, "P%imbar_newDTvQmax_ifQmax_%i", Pmbar, data_set);
  TCanvas *can5 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can5->Divide(ncol, 2);

  sprintf(name, "P%imbar_newDTvQ1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can6 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can6->Divide(ncol, 2);

  sprintf(name, "P%imbar_Q2Q3vQ1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can7 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can7->Divide(ncol, 2);

  sprintf(name, "P%imbar_newQ2Q3vQ1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can8 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can8->Divide(ncol, 2);

#if PLOT_MEAN_SIGNALS
  sprintf(name, "P%imbar_fSampler_Signal_A10_%i", Pmbar, data_set);
  TCanvas *can9 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can9->Divide(4, 2);
#endif

  sprintf(name, "P%imbar_newQ4vQ1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can10 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can10->Divide(ncol, 2);

  sprintf(name, "P%imbar_newQ4Q2vQ1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can11 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can11->Divide(ncol, 2);

  sprintf(name, "P%imbar_newQ4Q3vQ1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can12 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can12->Divide(ncol, 2);

  sprintf(name, "P%imbar_mult_vs_Q1_ifQmax_%i", Pmbar, data_set);
  TCanvas *can13 = new TCanvas(name, name, 0, 0, 2000, 1500);
  can13->Divide(ncol, 2);

  for (int a = 0; a < nA; a++) {
    can2->cd(a + 1);
    gPad->SetLogy();
    gPad->SetLogx();
    h1_Q1_ifQmax[a]->Draw();
    can3->cd(a + 1);
    gPad->SetLogz();
    h2_DT_vs_Q1_ifQmax[a]->Draw("colz");
    can4->cd(a + 1);
    gPad->SetLogz();
    h2_DT_vs_Qmax_ifQmax[a]->Draw("colz");
    can5->cd(a + 1);
    gPad->SetLogz();
    h2_newDT_vs_Qmax_ifQmax[a]->Draw("colz");
    can6->cd(a + 1);
    gPad->SetLogz();
    h2_newDT_vs_Q1_ifQmax[a]->Draw("colz");
    can7->cd(a + 1);
    gPad->SetLogz();
    h2_Discri[a]->Draw("colz");
    can8->cd(a + 1);
    gPad->SetLogz();
    h2_newDiscri[a]->Draw("colz");
    can10->cd(a + 1);
    gPad->SetLogz();
    h2_newQ4vQ1[a]->Draw("colz");
    can11->cd(a + 1);
    gPad->SetLogz();
    h2_newQ4Q2vQ1[a]->Draw("colz");
    can12->cd(a + 1);
    gPad->SetLogz();
    h2_newQ4Q3vQ1[a]->Draw("colz");
    can13->cd(a + 1);
    gPad->SetLogz();
    h2_MvQ1_ifQmax[a]->Draw("colz");

#if PLOT_MEAN_SIGNALS
    if (mapanode[a] == 10) {
      can7->cd(a + 1);
      gPad->SetLogx();
      for (int c = 0; c < 8; c++) {
        if (cutg[c]) {
          can7->cd(a + 1);
          cutg[c]->SetLineWidth(2);
          cutg[c]->SetLineColor(c + 1);
          cutg[c]->Draw("same");
        }
        can9->cd(c + 1);
        gPad->SetLogz();
        gPad->SetGridx();
        gPad->SetGridy();
        h2_sig[c]->Draw("colz");
      }
    }
#endif
  }
  can2->cd(nA + 1);
  gPad->SetLogy();
  gPad->SetLogz();
  h2_Q1_ifQmax->Draw("colz");
}
