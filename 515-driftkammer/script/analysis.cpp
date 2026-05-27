#include "RtypesCore.h"
#include "TH1.h"
#include <vector>
#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TRint.h>


const double TIME_LB  = -1.25;
const double TIME_UB  = 250 * 2.5 + 2.5 / 2.;
const unsigned int TIME_N = 251;
const double WIRE_LB = 0.5;
const double WIRE_UB = 48.5;
const unsigned int WIRE_N = 48;

// i feel guilty for this but i am also too lazy to repeat myself
#define TIME_BINS TIME_N, TIME_LB, TIME_UB
#define WIRE_BINS WIRE_N, WIRE_LB, WIRE_UB


void analysis::reset_entry_count() {
   this->current_entry = 0;
}


bool analysis::filter_exclude(unsigned int hit) {
   bool hit_too_late = time_le[hit] * 2.5 > 300;
   bool tot_too_short = tot[hit] * 2.5 < 100;
   return filter_enabled && (hit_too_late || tot_too_short);
}


bool analysis::get_next_entry() {
   // gets next entry and returns whether or not there are more
   if (this->n_entries == -1) 
      this->n_entries = fChain->GetEntriesFast();

   if (this->current_entry < this->n_entries) {
      this->GetEntry(this->current_entry);
      return this->current_entry++ < this->n_entries;
   }
   return false;
}


TH1D analysis::dt_relation() {
   TH1D drift_time_hist = TH1D("Driftzeiten", "Driftzeiten", TIME_BINS);

   for (reset_entry_count(); get_next_entry();) {
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         Double_t time = time_le[hit] * 2.5;
	       drift_time_hist.Fill(time);
	    }
   }

   return drift_time_hist;
}


TH2D analysis::wire_correlation() {
   TH2D wire_correlation = TH2D("wireCorrelation", "wire correlations", WIRE_BINS, WIRE_BINS);
   for (reset_entry_count(); get_next_entry();) {
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         if (filter_exclude(hit)) continue;
         for (UInt_t j = 0; j<nhits_le; j++) {
            if (hit == j) continue;
            if (wire_le[hit] == wire_le[j]) continue;
            wire_correlation.Fill(wire_le[hit], wire_le[j]);
         	}
	    }
   }
   return wire_correlation;
}


TH2D analysis::tot_wire_hist() {
   TH2D tot_hist = TH2D("tot_wire_hist", "Time over Treshhold per wire", WIRE_BINS, TIME_BINS);

   for (reset_entry_count(); get_next_entry();) {
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         if (filter_exclude(hit)) continue;
         Double_t time = tot[hit] * 2.5;
         int wire = wire_le[hit];
         if (time < 5) continue;
	       tot_hist.Fill(wire, time);
	    }
   }
   return tot_hist;
}


TH2D analysis::dt_tot_relation() {
   TH2D hist = TH2D("dt_tot_relation", "Driftzeit / TOT Relation", TIME_BINS, TIME_BINS);

   for (reset_entry_count(); get_next_entry();) {
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         if (filter_exclude(hit)) continue;
         Double_t time = this->tot[hit] * 2.5;
         Double_t dt = this->time_le[hit] * 2.5;
	       hist.Fill(dt, time);
	    }
   }
   return hist;
}


void analysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L analysis.C
//      Root > analysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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
   TH1D* driftTimesHisto = new TH1D("Driftzeiten", "Driftzeiten", 251, -2.5 / 2., 250 * 2.5 + 2.5 / 2.);
   TH2 *wireCorrHisto = new TH2D("wireCorrelation", "wire correlations", 48, 0.5, 48.5, 48, 0.5, 48.5);

   if (fChain == 0) return;

   for (;this->get_next_entry();) {
      
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         Double_t time = time_le[hit] * 2.5;
	       driftTimesHisto->Fill(time);

         for (UInt_t j = 0; j<nhits_le; j++) {
            if (hit == j) {
               continue;
            }
            wireCorrHisto->Fill(wire_le[hit],wire_le[j]);
         	}
      }

      // if (Cut(ientry) < 0) continue;
   }
   driftTimesHisto->GetXaxis()->SetTitle("Zeit / ns");
   driftTimesHisto->GetYaxis()->SetTitle("Trefferanzahl");
   //gStyle->SetOptStat(0);
   driftTimesHisto->Draw();
}


std::vector<int> make_bin_lut(TH2D& wire_correlation) {
   auto res = std::vector<int>(48);
   for (int i = 0; i < 48; ++i) {
      Long64_t max_before_i, max_after_i = 0;
      for (int j = 0; j < i; ++j) {
         auto elem = wire_correlation.GetBinContent(i, j);
         auto current_max = wire_correlation.GetBinContent(i, max_before_i);
         max_before_i = elem > current_max ? j : max_before_i;
      }

      for (int j = i; j < 48; ++j) {
         auto elem = wire_correlation.GetBinContent(i, j);
         auto current_max = wire_correlation.GetBinContent(i, max_after_i);
         max_after_i = elem > current_max ? j : max_after_i;
      }

      printf("for row %d: %lld and %lld\n", i, max_before_i, max_after_i);
   }
   return res;
}


int main(int argc, char** argv) {
   TROOT root("app","app");
   Int_t dargc=1;
   char** dargv = &argv[0];
   TRint app = TRint("app", &dargc, dargv);
   TCanvas c1 = TCanvas("c", "c", 800, 600);
   TFile f = TFile(argv[1]);
   TTree* tree = (TTree*) f.FindObjectAny("t");
   analysis* ana = new analysis(tree);
   // ana->Loop();
   auto dt_rel = ana->dt_relation();
   // dt_rel.Draw();
   auto tot_plot = ana->tot_wire_hist();
   // tot_plot.Draw();
   auto wire_correlation = ana->wire_correlation();
   wire_correlation.Draw();
   // auto _ = make_bin_lut(wire_correlation);

   auto dt_tot = ana->dt_tot_relation();
   // dt_tot.Draw();

   app.Run(kTRUE);
}
