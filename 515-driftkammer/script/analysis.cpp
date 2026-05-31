#include "RtypesCore.h"
#include "TH1.h"
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <string>
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
const unsigned int SPACE_N = TIME_N;

// i feel guilty for this but i am also too lazy to repeat myself
#define TIME_BINS TIME_N, TIME_LB, TIME_UB
#define WIRE_BINS WIRE_N, WIRE_LB, WIRE_UB
#define SPACE_BINS SPACE_N, -8.5, 8.5


struct checklist_64{
   uint64_t content = 0;

   bool check(unsigned char i) {
      bool was_checked = (1 << i) & content;
      content |= (1 << i);
      return !was_checked;
   }

   bool check_static(unsigned char i) {
      return (1 << i) & content;
   }
};


void analysis::reset_entry_count() {
   this->current_entry = 0;
}


bool analysis::filter_exclude(unsigned int hit) {
   bool hit_too_late = time_le[hit] * 2.5 > max_le_time;
   bool hit_too_early = time_le[hit] * 2.5 < min_le_time;
   bool tot_too_short = tot[hit] * 2.5 < min_tot;
   bool not_first_hit = false;
   forr (i, 0, hit) not_first_hit |= wire_le[hit] == wire_le[i];
   return filter_enabled && (hit_too_late || tot_too_short || not_first_hit || hit_too_early);
   // bool is_longest_tot_for_wire = true;
   // forr (i, 0, nhits_le) {
   //    if ((tot[hit] < tot[i]) && (i != hit) && (wire_le[hit] == wire_le[i]))
   //       is_longest_tot_for_wire = false;
   // }
   // bool not_longest = !is_longest_tot_for_wire;
   // return filter_enabled && (hit_too_late || tot_too_short || not_longest || hit_too_early);
}


bool analysis::get_next_entry() {
   // gets next entry and returns whether or not there are more
   if (this->n_entries == -1) 
      this->n_entries = fChain->GetEntriesFast();

   if (this->current_entry < this->n_entries) {
      this->GetEntry(this->current_entry);

      for_range(j, 0, nhits_le) this->wire_le[j] = this->wire_lut[wire_le[j]];
      for_range(j, 0, nhits_te) this->wire_te[j] = this->wire_lut[wire_te[j]];
      argsort();
      return this->current_entry++ < this->n_entries;
   }
   return false;
}


template<typename T>
void c_argsort(const T* array, unsigned char* indices, size_t len) {
    std::iota(indices, indices + len, 0);
    std::sort(indices, indices + len, [&array](int left, int right) -> bool { return array[left] < array[right]; });
}


void analysis::argsort() {
   c_argsort(wire_le, sorted, nhits_le);
}


TH1D analysis::dt_hist() {
   TH1D drift_time_hist = TH1D("Driftzeiten", "Driftzeiten", TIME_BINS);

   reset_entry_count();
   while (get_next_entry()) {
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
   reset_entry_count();
   while (get_next_entry()) {
      for_range(hit, 0, nhits_le) {
         if (filter_exclude(hit)) continue;
         for_range(j, 0, nhits_le) {
            if (hit == j) continue;
            // if (wire_le[hit] == wire_le[j]) continue;
            wire_correlation.Fill(wire_le[hit], wire_le[j]);
         	}
	    }
   }
   return wire_correlation;
}


TH2D analysis::dt_wire_hist() {
   TH2D tot_hist = TH2D("dt_wire_hist", "Driftzeit pro Draht", WIRE_BINS, TIME_BINS);

   reset_entry_count();
   while (get_next_entry()) {
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         if (filter_exclude(hit)) continue;
         Double_t time = time_le[hit] * 2.5;
         int wire = wire_le[hit];
         if (time < 5) continue;
	       tot_hist.Fill(wire, time);
	    }
   }
   return tot_hist;
}


TH2D analysis::dt_tot_relation() {
   TH2D hist = TH2D("dt_tot_relation", "Driftzeit / TOT Relation", TIME_BINS, TIME_BINS);

   reset_entry_count();
   while (get_next_entry()) {
      for(UInt_t hit = 0; hit < nhits_le; hit++) {
         if (filter_exclude(hit)) continue;
         Double_t time = this->tot[hit] * 2.5;
         Double_t dt = this->time_le[hit] * 2.5;
	       hist.Fill(dt, time);
	    }
   }
   return hist;
}

TH1D analysis::basic_angle_distrib() {
   auto cell_edges_to_angles = [](Double_t x) -> Double_t {
      const auto scint_pos = 0;
      x -= scint_pos;
      x *= 17.; // to position
      return 17. * x;
   };

   Double_t bin_edges[25];
   std::iota(bin_edges, bin_edges + 25, 0);
   std::transform(bin_edges, bin_edges + 25, bin_edges, cell_edges_to_angles);

   TH1D angle_distribution = TH1D("basic_angle_distribution", "Winkelverteilung der Kosmischen Strahlung", 24, bin_edges);

   reset_entry_count();
   while(get_next_entry()) {
      bool in_seq = false;
      unsigned char block_start;
      forr (i, 0, nhits_le - 1) {
         auto hit = sorted[i];
         auto hit_next = sorted[i + 1];
         if (filter_exclude(hit)) continue;
         auto j = 2;
         while (filter_exclude(hit_next)) {
            hit_next = sorted[i + j];
            if (i + j < nhits_le) j ++;
            else break;
         }

         bool sequential = wire_le[hit] + 1 == wire_le[hit_next];
         bool sequential_same_layer = wire_le[hit] + 2 == wire_le[hit_next];
         bool seq = sequential || sequential_same_layer;

         if (in_seq && seq) continue;
         if (!in_seq && seq) {
            in_seq = true;
            block_start = wire_le[hit];
         }
         if (in_seq && !seq) {
            in_seq = false;
            printf("block from %u to %u\n", block_start, wire_le[hit]);
         }
      }
   }

   return angle_distribution;
}


std::vector<UInt_t> make_wire_lut() {
   auto res = std::vector<UInt_t>(49);
   forr(i, 0, 49) res[i] = (i % 2? i + 49 + 1: i + 49 - 1) % 49;
   // forr(i, 0, 49) printf("%u\n", res[i]);
   return res;
}


TH1D make_odb(TH1D& drift_time_spectrum) {
   TH1D odb = TH1D("odb", "Orts- Driftzeitbeziehung", TIME_BINS);

   Double_t sum = 0;
   forr(i, 1, SPACE_N + 1) {
      sum += drift_time_spectrum.GetBinContent(i);
      odb.SetBinContent(i, sum);
   }

   odb.Scale(8.5 / sum);
   return odb;
}


TH2D analysis::dist_plot(TH1D& odb, unsigned int wire_lb, unsigned int wire_ub) {
   TH2D dists = TH2D ("dists_plot", "TODO", 40, -4.5, 4.5, 40, -0.2, 17.2);
   reset_entry_count();
   while(get_next_entry()) {
      checklist_64 hit_wires;
      forr(hit, 0, nhits_le) hit_wires.check(wire_le[hit]);
      forr (hit, 0, nhits_le) {
         if (filter_exclude(hit)) continue;
         if ((wire_le[hit] < wire_lb) || (wire_le[hit] > wire_ub)) continue;
         forr (i, 0, nhits_le) {
            if (filter_exclude(i)) continue;
            if (wire_le[hit] + 1 == wire_le[i]) {
               unsigned int time_a, time_b;
               time_a = time_le[hit];
               time_b = time_le[i];
               Double_t dist_a = odb.At(time_a);
               Double_t dist_b = odb.At(time_b);


               dists.Fill(0.5 * (dist_a - dist_b), (dist_a + dist_b));

               // if (0.5 * (dist_a - dist_b) < 0. && (dist_a + dist_b) < 0.8)
               //    printf("%u %u %u %u %u %u %d %d\n", wire_le[hit], wire_le[i], time_le[hit], time_le[i], tot[hit], tot[i], hit_wires.check_static(wire_le[hit] - 2), hit_wires.check_static(wire_le[hit] + 2));
               break;
            }
         }
      }

   }
   return dists;
}


template<typename Plotable>
void plot(Plotable& p, TCanvas* C) {
   C->cd();
   p.Draw();
   C->Update();
}


int main(int argc, char** argv) {
   // weniger statistik mit guten parametern ..161632.root
   // viel statistik mit guten parametern ...161826.root
   TROOT root("app","app");
   Int_t dargc=1;
   char** dargv = &argv[0];
   TRint app = TRint("app", &dargc, dargv);
   TFile f = TFile(argv[1]);
   TTree* tree = (TTree*) f.FindObjectAny("t");
   analysis* ana = new analysis(tree);

   TFile calib_file = TFile("../data/B103/run_260513_161826.root"); 
   TTree* calib_tree = (TTree*) calib_file.FindObjectAny("t");
   analysis* calib_data = new analysis(calib_tree);
   calib_data->max_le_time = 500;


   ana->wire_lut = make_wire_lut();
   // auto dt_rel = ana->dt_hist();
   auto dt_rel = calib_data->dt_hist();
   auto dt_wire_plot = ana->dt_wire_hist();
   auto wire_correlation = ana->wire_correlation();
   auto dt_tot = ana->dt_tot_relation();

   std::vector<TCanvas*> canvas_vec = std::vector<TCanvas*>();
   forr(i, 0, 7) canvas_vec.push_back(new TCanvas(("c" + std::to_string(i)).c_str(), "c", 800, 600));


   plot(dt_rel, canvas_vec[0]);
   plot(dt_wire_plot, canvas_vec[1]);
   plot(wire_correlation, canvas_vec[2]);
   plot(dt_tot, canvas_vec[3]);

   // TODO: maybe this needs to be done before filtering??
   auto odb = make_odb(dt_rel);
   plot(odb, canvas_vec[4]);

   auto basic_angles = ana->basic_angle_distrib();
   plot(basic_angles, canvas_vec[5]);

   auto sum_vs_diff = ana->dist_plot(odb, 12, 30);
   plot(sum_vs_diff, canvas_vec[6]);

   app.Run(kTRUE);
}
