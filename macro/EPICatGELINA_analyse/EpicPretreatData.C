#include "EpicPretreatData.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

EpicPretreatData::EpicPretreatData(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../output/pretreat/pretreat_run17.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../output/pretreat/pretreat_run17.root");
      }
      f->GetObject("EpicPretreatData",tree);

   }
   Init(tree);
}

EpicPretreatData::~EpicPretreatData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EpicPretreatData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EpicPretreatData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EpicPretreatData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fFC_anode", &fFC_anode, &b_fFC_anode);
   fChain->SetBranchAddress("fFC_time", &fFC_time, &b_fFC_time);
   fChain->SetBranchAddress("fFC_qmax", &fFC_qmax, &b_fFC_qmax);
   fChain->SetBranchAddress("fFC_q1", &fFC_q1, &b_fFC_q1);
   fChain->SetBranchAddress("fFC_q2", &fFC_q2, &b_fFC_q2);
   fChain->SetBranchAddress("fFC_q3", &fFC_q3, &b_fFC_q3);
   fChain->SetBranchAddress("fFC_alpha", &fFC_alpha, &b_fFC_alpha);
   fChain->SetBranchAddress("fFC_fission", &fFC_fission, &b_fFC_fission);
   fChain->SetBranchAddress("fHF_time_curr", &fHF_time_curr, &b_fHF_time_curr);
   fChain->SetBranchAddress("fHF_time_next", &fHF_time_next, &b_fHF_time_next);
   Notify();
}

bool EpicPretreatData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be for a new TTree in a TChain. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void EpicPretreatData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EpicPretreatData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void EpicPretreatData::Loop()
{
//   In a ROOT session, you can do:
//      root> .L EpicPretreatData.C
//      root> EpicPretreatData t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
