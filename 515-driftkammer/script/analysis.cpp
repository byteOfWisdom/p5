#include "RtypesCore.h"
#include "TColor.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TString.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
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
#include <TMath.h>
#include <TF1.h>
#include <TFitResult.h>

#define for_range(I, A, B) for (auto I = A; I < B; ++I)
#define iterate(I, A, B) for (auto I = A; I != A + B; ++I)
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


class dataset{
   private:
   TTree* tree;
   TFile* f;

   public:
   analysis* ana;
   dataset(TString fname) {
      this->f = new TFile(fname);
      this->tree = (TTree*) f->FindObjectAny("t");
      this->ana = new analysis(this->tree);
   };

  ~dataset() {
      delete this->ana;
      delete this->tree;     
      delete this->f;
  }; 
};


class global_object_store{
   public:
      global_object_store() {};
      TH1D* hold(TH1D obj) {
         this->hists.push_back(obj);
         return &this->hists.back();
      };

      TH2D* hold(TH2D obj){
         this->hists2.push_back(obj);
         return &this->hists2.back();
      };

      dataset* hold(dataset obj){
         this->datasets.push_back(obj);
         return &this->datasets.back();
      };

   private:
      std::vector<TH1D> hists;
      std::vector<TH2D> hists2;
      std::vector<dataset> datasets;
};



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
      forr (i, 0, nhits_le) {
         if (filter_exclude(i)) continue;
         valid[n_valid] = i;
         n_valid ++;
      }
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


TH1D analysis::dt_hist(TString name = "Driftzeiten") {
   TH1D drift_time_hist = TH1D(name, "Driftzeiten", TIME_BINS);

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

inline unsigned char abs_diff(unsigned char a, unsigned char b) {
   return abs((int) a - (int) b); }


TH1D analysis::basic_angle_distrib() {
   const auto middle_bin = 20;
   auto cell_edges_to_angles = [](Double_t x) -> Double_t {
      const auto scint_pos = 0.5 * middle_bin - 0.5;
      const auto height_diff = 12.5e-2;
      auto d = (x - scint_pos) * 17.e-3;
      auto theta = atan(d / height_diff);
      const auto rad_to_deg = 180 / TMath::Pi();
      return theta * rad_to_deg;
   };

   Double_t bin_edges[25];
   std::iota(bin_edges, bin_edges + 25, 0);
   std::transform(bin_edges, bin_edges + 25, bin_edges, cell_edges_to_angles);

   // TH1D angle_distribution = TH1D("basic_angle_distribution", "Winkelverteilung der Kosmischen Strahlung", 24, bin_edges);
   TH1D block_starts = TH1D("block_starts", "Winkelverteilung der Kosmischen Strahlung", WIRE_BINS);

   reset_entry_count();
   while(get_next_entry()) {
      bool in_seq = false;
      unsigned char block_start;
      iterate(hit, valid, n_valid - 1){
         bool sequential = wire_le[*hit] + 1 == wire_le[*(hit + 1)];
         bool sequential_same_layer = wire_le[*hit] + 2 == wire_le[*(hit + 1)];
         bool seq = sequential || sequential_same_layer;

         if (in_seq && seq) continue;
         if (!in_seq && seq) {
            in_seq = true;
            block_start = wire_le[*hit];
         }
         if (in_seq && !seq) {
            in_seq = false;
            auto first_hit = abs_diff(block_start, middle_bin) <= abs_diff(wire_le[*hit], middle_bin)? block_start: wire_le[*hit];
            block_starts.Fill(first_hit);
            // printf("block from %u to %u, first hit at %u\n", block_start, wire_le[*hit], first_hit);
         }
      }
   }


   TH1D angle_distribution = TH1D("basic_angle_distribution", "Winkelverteilung der Kosmischen Strahlung", 24, bin_edges);
   forr (wire, 1, 49) {
      auto count = block_starts.GetBinContent(wire);
      int direction = wire < middle_bin? 1: -1;
      int target_bin = wire % 2? wire / 2 + direction: wire / 2;
      if (wire % 2) continue;
      angle_distribution.AddBinContent(target_bin, count);
   }

   return angle_distribution;
}


std::vector<UInt_t> make_wire_lut() {
   auto res = std::vector<UInt_t>(49);
   forr(i, 0, 49) res[i] = (i % 2? i + 49 + 1: i + 49 - 1) % 49;
   // forr(i, 0, 49) printf("%u %u\n", i, res[i]);
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


void plot_set(std::vector<dataset*>& files, TCanvas* canv, global_object_store* gob) {
   canv->cd();
   EColor colors[] = {kBlue, kPink, kRed, kOrange};
   int i = 0;
   printf("entering func\n");
   for (dataset* ds : files) {
      ds->ana->filter_enabled = false;
      auto hist = gob->hold(ds->ana->dt_hist("hist" + std::to_string(i)));

      hist->SetLineColor(colors[i]);
      if (i == 0) hist->DrawCopy("HIST");
      else hist->DrawCopy("HIST same");

      i++;
   }

   canv->Update();
}


int main(int argc, char** argv) {
   // weniger statistik mit guten parametern ..161632.root
   // viel statistik mit guten parametern ...161826.root
   TROOT root("app","app");
   Int_t dargc=1;
   char** dargv = &argv[0];
   TRint app = TRint("app", &dargc, dargv);
   // TFile f = TFile(argv[1]);
   // TTree* tree = (TTree*) f.FindObjectAny("t");
   // analysis* ana = new analysis(tree);
   dataset* data = new dataset(argv[1]);
   analysis* ana = data->ana;
   // analysis* ana = load_file(argv[1]);

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

   TF1 angle_dist_func = TF1("angle_dist_func", "[0] * cos([1] * x * pi / 180)^2", -90, 90);

   auto basic_angles = ana->basic_angle_distrib();
   TFitResultPtr fit = basic_angles.Fit(&angle_dist_func, "S");
   plot(basic_angles, canvas_vec[5]);
   plot(*fit, canvas_vec[5]);

   auto sum_vs_diff = ana->dist_plot(odb, 12, 30);
   plot(sum_vs_diff, canvas_vec[6]);


   TCanvas* voltage_canvas = new TCanvas("voltage_sweep", "Verschiedene Beschleunigungsspannungen");
   std::vector<dataset*> voltage_data = std::vector<dataset*> ({
      new dataset("../data/B103/run_260513_145435.root"),
      new dataset("../data/B103/run_260513_145148.root"),
      new dataset("../data/B103/run_260513_144851.root"),
   });

   global_object_store* gob = new global_object_store();
   plot_set(voltage_data, voltage_canvas, gob);


   app.Run(kTRUE);
}
