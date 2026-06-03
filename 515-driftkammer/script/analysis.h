//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  1 12:17:14 2015 by ROOT version 5.34/03
// from TTree t/t
// found on file: run_150401_120548.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include "RtypesCore.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <vector>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

struct hit_chunk {
   unsigned char n_hits;
   unsigned char wire[8];
   unsigned char tot[8];
   unsigned char time_le[8];
};


class analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          event;
   Double_t        eventTime;
   Double_t        deltaT;
   UInt_t          hwEvent;
   UInt_t          nBytes;
   UInt_t          c3trigTime;
   UInt_t          c3trigNr;
   UInt_t          nhits_le;
   UInt_t          wire_le[1000];   //[nhits_le]
   UInt_t          time_le[1000];   //[nhits_le]
   UInt_t          nhits_te;
   UInt_t          wire_te[1000];   //[nhits_te]
   UInt_t          time_te[1000];   //[nhits_te]
   UInt_t          tot[1000];   //[nhits_le]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_eventTime;   //!
   TBranch        *b_deltaT;   //!
   TBranch        *b_hwEvent;   //!
   TBranch        *b_nBytes;   //!
   TBranch        *b_c3trigTime;   //!
   TBranch        *b_c3trigNr;   //!
   TBranch        *b_nhits_le;   //!
   TBranch        *b_wire_le;   //!
   TBranch        *b_time_le;   //!
   TBranch        *b_nhits_te;   //!
   TBranch        *b_wire_te;   //!
   TBranch        *b_time_te;   //!
   TBranch        *b_tot;   //!


   unsigned char sorted[100];
   unsigned char valid[100];
   unsigned char n_valid;
   bool filter_enabled = true;
   Double_t max_le_time = 300;
   Double_t min_le_time = 0;
   Double_t min_tot = 100;
   std::vector<UInt_t> wire_lut = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48};

   analysis(TTree *tree=0);
   virtual ~analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop() {};
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool get_next_entry();
   virtual bool filter_exclude(unsigned int);
   virtual void reset_entry_count();
   virtual TH1D dt_hist(TString);
   virtual TH2D wire_correlation();
   virtual TH2D dt_wire_hist();
   virtual TH2D dt_tot_relation();
   virtual TH1D basic_angle_distrib();
   virtual TH2D dist_plot(TH1D&, unsigned int, unsigned int);
   virtual void argsort();
   Double_t get_runtime();

   private:
   Long64_t current_entry = 0;
   Long64_t n_entries = -1;
};

#endif

#ifdef analysis_cxx
analysis::analysis(TTree *tree) : fChain(0) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run_150401_120548.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("run_150401_120548.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}


analysis::~analysis() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


Int_t analysis::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}


Long64_t analysis::LoadTree(Long64_t entry) {
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


void analysis::Init(TTree *tree) {
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("eventTime", &eventTime, &b_eventTime);
   fChain->SetBranchAddress("deltaT", &deltaT, &b_deltaT);
   fChain->SetBranchAddress("hwEvent", &hwEvent, &b_hwEvent);
   fChain->SetBranchAddress("nBytes", &nBytes, &b_nBytes);
   fChain->SetBranchAddress("c3trigTime", &c3trigTime, &b_c3trigTime);
   fChain->SetBranchAddress("c3trigNr", &c3trigNr, &b_c3trigNr);
   fChain->SetBranchAddress("nhits_le", &nhits_le, &b_nhits_le);
   fChain->SetBranchAddress("wire_le", &wire_le, &b_wire_le);
   fChain->SetBranchAddress("time_le", &time_le, &b_time_le);
   fChain->SetBranchAddress("nhits_te", &nhits_te, &b_nhits_te);
   fChain->SetBranchAddress("wire_te", &wire_te, &b_wire_te);
   fChain->SetBranchAddress("time_te", &time_te, &b_time_te);
   fChain->SetBranchAddress("tot", &tot, &b_tot);
   Notify();
}

Bool_t analysis::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}


Int_t analysis::Cut(Long64_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_cxx
