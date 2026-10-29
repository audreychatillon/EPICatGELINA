//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 29 18:18:52 2026 by ROOT version 6.40.04
// from TTree EpicPretreatData/no cross talk, no reflexion, HF_curr, HF_next
// found on file: ../../output/pretreat/pretreat_run17.root
//////////////////////////////////////////////////////////

#ifndef EpicPretreatData_h
#define EpicPretreatData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class EpicPretreatData {
public :
   TTree          *fChain;   ///<!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; ///<!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Short_t         fFC_anode;
   Double_t        fFC_time;
   Double_t        fFC_qmax;
   Double_t        fFC_q1;
   Double_t        fFC_q2;
   Double_t        fFC_q3;
   Bool_t          fFC_alpha;
   Bool_t          fFC_fission;
   Double_t        fHF_time_curr;
   Double_t        fHF_time_next;

   // List of branches
   TBranch        *b_fFC_anode;   ///<!
   TBranch        *b_fFC_time;   ///<!
   TBranch        *b_fFC_qmax;   ///<!
   TBranch        *b_fFC_q1;   ///<!
   TBranch        *b_fFC_q2;   ///<!
   TBranch        *b_fFC_q3;   ///<!
   TBranch        *b_fFC_alpha;   ///<!
   TBranch        *b_fFC_fission;   ///<!
   TBranch        *b_fHF_time_curr;   ///<!
   TBranch        *b_fHF_time_next;   ///<!

   EpicPretreatData(TTree *tree=0);
   virtual ~EpicPretreatData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
