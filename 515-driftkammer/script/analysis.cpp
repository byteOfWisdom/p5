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

#define for_range(I, A, B) for (auto I = A; I < B; ++I)
#define forr for_range

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

      for_range(j, 0, nhits_le) this->wire_le[j] = this->wire_lut[wire_le[j]];
      for_range(j, 0, nhits_te) this->wire_te[j] = this->wire_lut[wire_te[j]];
      return this->current_entry++ < this->n_entries;
   }
   return false;
}


TH1D analysis::dt_relation() {
   TH1D drift_time_hist = TH1D("Driftzeiten", "Driftzeiten", TIME_BINS);

   for (reset_entry_count(); get_next_entry();) {
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         if (filter_exclude(hit)) continue;
         Double_t time = time_le[hit] * 2.5;
	       drift_time_hist.Fill(time);
	    }
   }

   return drift_time_hist;
}


TH2D analysis::wire_correlation() {
   TH2D wire_correlation = TH2D("wireCorrelation", "wire correlations", WIRE_BINS, WIRE_BINS);
   for (reset_entry_count(); get_next_entry();) {
      for_range(hit, 0, nhits_le) {
         if (filter_exclude(hit)) continue;
         for_range(j, 0, nhits_le) {
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


std::vector<UInt_t> make_wire_lut() {
   auto res = std::vector<UInt_t>(48);
   forr(i, 0, 48) res[i] = (i % 2? i + 48 + 1: i + 48 - 1) % 48;
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
   ana->wire_lut = make_wire_lut();
   auto dt_rel = ana->dt_relation();
   auto tot_plot = ana->tot_wire_hist();
   auto wire_correlation = ana->wire_correlation();
   auto dt_tot = ana->dt_tot_relation();


   // dt_rel.Draw();
   // tot_plot.Draw();
   wire_correlation.Draw();
   // dt_tot.Draw();

   app.Run(kTRUE);
}
