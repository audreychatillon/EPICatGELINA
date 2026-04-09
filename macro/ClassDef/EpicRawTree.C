#include "EpicRawTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

EpicRawTree::EpicRawTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../output/conversion/raw_run19.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../output/conversion/raw_run19.root");
      }
      f->GetObject("EpicRawTree",tree);

   }
   Init(tree);
}

EpicRawTree::~EpicRawTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EpicRawTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EpicRawTree::LoadTree(Long64_t entry)
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

void EpicRawTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fFC_DetNbr", &fFC_DetNbr, &b_epic_fFC_DetNbr);
   fChain->SetBranchAddress("fFC_AnodeNbr", &fFC_AnodeNbr, &b_epic_fFC_AnodeNbr);
   fChain->SetBranchAddress("fFC_PulserTrig", &fFC_PulserTrig, &b_epic_fFC_PulserTrig);
   fChain->SetBranchAddress("fFC_Time", &fFC_Time, &b_epic_fFC_Time);
   fChain->SetBranchAddress("fFC_TofRaw", &fFC_TofRaw, &b_epic_fFC_TofRaw);
   fChain->SetBranchAddress("fFC_TimeCfd", &fFC_TimeCfd, &b_epic_fFC_TimeCfd);
   fChain->SetBranchAddress("fFC_TimeQmax", &fFC_TimeQmax, &b_epic_fFC_TimeQmax);
   fChain->SetBranchAddress("fFC_Qmax", &fFC_Qmax, &b_epic_fFC_Qmax);
   fChain->SetBranchAddress("fFC_Q1", &fFC_Q1, &b_epic_fFC_Q1);
   fChain->SetBranchAddress("fFC_Q2", &fFC_Q2, &b_epic_fFC_Q2);
   fChain->SetBranchAddress("fFC_Q3", &fFC_Q3, &b_epic_fFC_Q3);
   fChain->SetBranchAddress("fQmax_Index", &fQmax_Index, &b_epic_fQmax_Index);
   fChain->SetBranchAddress("fQmax_Sampler", &fQmax_Sampler, &b_epic_fQmax_Sampler);
   fChain->SetBranchAddress("fHF_Time", &fHF_Time, &b_epic_fHF_Time);
   fChain->SetBranchAddress("fHF_TimePrev", &fHF_TimePrev, &b_epic_fHF_TimePrev);
   Notify();
}

bool EpicRawTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void EpicRawTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EpicRawTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void EpicRawTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L EpicRawTree.C
//      root> EpicRawTree t
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
