//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr  5 16:37:53 2026 by ROOT version 6.32.06
// from TTree EpicRawTree/nptool tree
// found on file: ../../output/conversion/raw_run19.root
//////////////////////////////////////////////////////////

#ifndef EpicRawTree_h
#define EpicRawTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class EpicRawTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<short>   fFC_DetNbr;
   vector<short>   fFC_AnodeNbr;
   vector<bool>    fFC_PulserTrig;
   vector<double>  fFC_Time;
   vector<double>  fFC_TofRaw;
   vector<double>  fFC_TimeCfd;
   vector<double>  fFC_TimeQmax;
   vector<double>  fFC_Qmax;
   vector<double>  fFC_Q1;
   vector<double>  fFC_Q2;
   vector<double>  fFC_Q3;
   Short_t         fQmax_Index;
   vector<double>  fQmax_Sampler;
   Double_t        fHF_Time;
   Double_t        fHF_TimePrev;

   // List of branches
   TBranch        *b_epic_fFC_DetNbr;   //!
   TBranch        *b_epic_fFC_AnodeNbr;   //!
   TBranch        *b_epic_fFC_PulserTrig;   //!
   TBranch        *b_epic_fFC_Time;   //!
   TBranch        *b_epic_fFC_TofRaw;   //!
   TBranch        *b_epic_fFC_TimeCfd;   //!
   TBranch        *b_epic_fFC_TimeQmax;   //!
   TBranch        *b_epic_fFC_Qmax;   //!
   TBranch        *b_epic_fFC_Q1;   //!
   TBranch        *b_epic_fFC_Q2;   //!
   TBranch        *b_epic_fFC_Q3;   //!
   TBranch        *b_epic_fQmax_Index;   //!
   TBranch        *b_epic_fQmax_Sampler;   //!
   TBranch        *b_epic_fHF_Time;   //!
   TBranch        *b_epic_fHF_TimePrev;   //!

   EpicRawTree(TTree *tree=0);
   virtual ~EpicRawTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
