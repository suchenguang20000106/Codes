//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Nov 17 13:55:44 2019 by ROOT version 5.34/38
// from TTree ppi0/Lc+ Study
// found on file: Lcppi0_loose2.root
//////////////////////////////////////////////////////////

#ifndef loose_h
#define loose_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class loose {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           n_Lc;
   Int_t           n_P;
   Int_t           n_Pi0;
   Int_t           n_Gam;
   Double_t        p4_Lc[1000][4];   //[n_Lc]
   Double_t        p4_P[1000][4];   //[n_P]
   Double_t        p4_Pi0[1000][4];   //[n_Pi0]
   Double_t        p4_Gam[1000][4];   //[n_Gam]
   Double_t        Mbc[1000];   //[n_Lc]
   Double_t        deltaE[1000];   //[n_Lc]
   Double_t        chisquare[1000];   //[n_Pi0]
   Double_t        Mpi0[1000];   //[n_Pi0]
   Double_t        p_Vr[1000];   //[n_P]
   Double_t        mdang[1000];   //[n_Gam]
   Double_t        mpdang[1000];   //[n_Gam]

   // List of branches
   TBranch        *b_n_Lc;   //!
   TBranch        *b_n_P;   //!
   TBranch        *b_n_Pi0;   //!
   TBranch        *b_n_Gam;   //!
   TBranch        *b_p4_Lc;   //!
   TBranch        *b_p4_P;   //!
   TBranch        *b_p4_Pi0;   //!
   TBranch        *b_p4_Gam;   //!
   TBranch        *b_Mbc;   //!
   TBranch        *b_deltaE;   //!
   TBranch        *b_chisquare;   //!
   TBranch        *b_Mpi0;   //!
   TBranch        *b_p_Vr;   //!
   TBranch        *b_mdang;   //!
   TBranch        *b_mpdang;   //!

   loose(TTree *tree=0);
   virtual ~loose();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef loose_cxx
loose::loose(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Lcppi0_loose3.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Lcppi0_loose3.root");
      }
      f->GetObject("ppi0",tree);

   }
   Init(tree);
}

loose::~loose()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t loose::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t loose::LoadTree(Long64_t entry)
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

void loose::Init(TTree *tree)
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

   fChain->SetBranchAddress("n_Lc", &n_Lc, &b_n_Lc);
   fChain->SetBranchAddress("n_P", &n_P, &b_n_P);
   fChain->SetBranchAddress("n_Pi0", &n_Pi0, &b_n_Pi0);
   fChain->SetBranchAddress("n_Gam", &n_Gam, &b_n_Gam);
   fChain->SetBranchAddress("p4_Lc", p4_Lc, &b_p4_Lc);
   fChain->SetBranchAddress("p4_P", p4_P, &b_p4_P);
   fChain->SetBranchAddress("p4_Pi0", p4_Pi0, &b_p4_Pi0);
   fChain->SetBranchAddress("p4_Gam", p4_Gam, &b_p4_Gam);
   fChain->SetBranchAddress("Mbc", Mbc, &b_Mbc);
   fChain->SetBranchAddress("deltaE", deltaE, &b_deltaE);
   fChain->SetBranchAddress("chisquare", chisquare, &b_chisquare);
   fChain->SetBranchAddress("Mpi0", Mpi0, &b_Mpi0);
   fChain->SetBranchAddress("p_Vr", p_Vr, &b_p_Vr);
   fChain->SetBranchAddress("mdang", mdang, &b_mdang);
   fChain->SetBranchAddress("mpdang", mpdang, &b_mpdang);
   Notify();
}

Bool_t loose::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void loose::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t loose::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef loose_cxx
