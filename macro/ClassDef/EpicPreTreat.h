//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 10 15:44:46 2026 by ROOT version 6.32.06
// from TTree EpicPreTreat/pretreated raw data, cross talk removed run3 P=1190mbar, HV=610V
// found on file: ../../output/pretreat/pretreat3_1190mbar_610V.root
//////////////////////////////////////////////////////////

#ifndef EpicPreTreat_h
#define EpicPreTreat_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class EpicPreTreat {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Short_t         fFC_DetNbr;
   Short_t         fFC_AnodeNbr;
   Double_t        fFC_TofRaw;
   Double_t        fFC_TimeCfd;
   Double_t        fFC_TimeQmax;
   Double_t        fFC_Qmax;
   Double_t        fFC_Q1;
   Double_t        fFC_Q2;
   Double_t        fFC_Q3;
   Double_t        fFC_Q4;
   vector<double>  *fQmax_Sampler;

   // List of branches
   TBranch        *b_fFC_DetNbr;   //!
   TBranch        *b_fFC_AnodeNbr;   //!
   TBranch        *b_fFC_TofRaw;   //!
   TBranch        *b_fFC_TimeCfd;   //!
   TBranch        *b_fFC_TimeQmax;   //!
   TBranch        *b_fFC_Qmax;   //!
   TBranch        *b_fFC_Q1;   //!
   TBranch        *b_fFC_Q2;   //!
   TBranch        *b_fFC_Q3;   //!
   TBranch        *b_fFC_Q4;   //!
   TBranch        *b_fQmax_Sampler;   //!

   EpicPreTreat(TTree *tree=0);
   virtual ~EpicPreTreat();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
