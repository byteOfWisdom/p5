#include "RtypesCore.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TString.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <numeric>
#include <string>
#include <tuple>
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
#include <TLegend.h>

#define for_range(I, A, B) for (auto I = A; I < B; ++I)
#define iterate(I, A, B) for (auto I = A; I != A + B; ++I)
#define forr for_range
#define IS_EVEN % 2 == 0
#define IS_ODD % 2 == 1

const double TIME_LB  = -1.25;
const double TIME_UB  = 250 * 2.5 + 2.5 / 2.;
const unsigned int TIME_N = 251;
const double WIRE_LB = 0.5;
const double WIRE_UB = 48.5;
const unsigned int WIRE_N = 48;
const unsigned int SPACE_N = TIME_N;
const auto SCINT_CENTERED_OVER = 22;

// i feel guilty for this but i am also too lazy to repeat myself
#define TIME_BINS TIME_N, TIME_LB, TIME_UB
#define WIRE_BINS WIRE_N, WIRE_LB, WIRE_UB
#define SPACE_BINS SPACE_N, -8.5, 8.5

typedef std::vector<unsigned int> wire_chunk_t;

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
   bool content[64];

   checklist_64() {
      forr(i, 0, 64) content[i] = false;
   };

   bool check(unsigned char i) {
      bool was_checked = content[i];
      content[i] = true;
      return !was_checked;
   }

   bool check_static(unsigned char i) {
      return content[i];
   }
};


void analysis::reset_entry_count() {
   this->current_entry = 0;
}


typedef enum {
   left,
   center,
   right
} position_in_cell_t;

constexpr Double_t cell_edges_to_angles(unsigned int n, position_in_cell_t edge) {
   // const auto middle_bin = 20;
   const auto middle_bin = SCINT_CENTERED_OVER;
   const auto height_bottom_row = 12.5e-2;
   const auto height_top_row = height_bottom_row - 1.7e-2; // TODO: check!!!
   const auto scint_edge_l = middle_bin * 8.5e-3;
   const auto scint_edge_r = (middle_bin + 1) * 8.5e-3;
   const auto scint_edge_c = (middle_bin + 0.5) * 8.5e-3;
   const auto rad_to_deg = 180 / M_PI;
   const auto tan_30_deg = 0.4492185466271425; // sufficient and then stuff can stay constexpr
   const auto delta_h = 8.5e-3 * tan_30_deg;

   auto a = n IS_ODD? height_bottom_row: height_top_row; 
   auto d = 0.;
   if (edge == left) {
      d = (n * 8.5e-3) - scint_edge_l;
      a = a - delta_h;
   }
   if (edge == right) {
      d = (n * 8.5e-3) - scint_edge_r;
      a = a - delta_h;
   }
   if (edge == center) {
      d = (n * 8.5e-3) - scint_edge_c;
   }

   const auto theta = atan(d / a);
   return theta * rad_to_deg;
}

Double_t delta_theta(unsigned int cell, Double_t r) {
   const auto middle_bin = SCINT_CENTERED_OVER;
   const auto height_bottom_row = 12.5e-2;
   const auto height_top_row = height_bottom_row - 1.7e-2; // TODO: check!!!
   const auto scint_edge_c = (middle_bin + 0.5) * 8.5e-3;
   const auto rad_to_deg = 180 / M_PI;

   auto a = cell IS_ODD? height_bottom_row: height_top_row; 
   auto d = (cell * 8.5e-3) - scint_edge_c;

   auto c = sqrt(a * a + d * d);

   const auto delta_theta = atan(r / c);
   return delta_theta * rad_to_deg;
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
      return this->current_entry++ < this->n_entries;
   }
   return false;
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
   return abs((int) a - (int) b);
}

inline unsigned char abs_diff(int a, int b) {
   return abs(a - b);
}


std::vector<std::vector<unsigned int>> get_sequences(checklist_64 hits) {
   std::vector<std::vector<unsigned int>> res = {};
   std::vector<unsigned int> current_chunk = {};

   forr (i, 0, 48) {
      bool immediate = hits.check_static(i) and hits.check_static(i + 1);
      bool skipped = hits.check_static(i) and hits.check_static(i + 2);
      if (immediate) {
         current_chunk.push_back(i);
      }
      else if (skipped) {
         current_chunk.push_back(i);
         i++;
      }

      else if (current_chunk.size() > 0){
         current_chunk.push_back(i);
         res.push_back(current_chunk);
         current_chunk = {};
      }
   }

   return res;
}


TH1D analysis::basic_angle_distrib() {
   const auto middle_bin = SCINT_CENTERED_OVER;
   // auto local_cell_edges_to_angles = [](Double_t x) -> Double_t {
   //    const auto scint_pos = 0.5 * middle_bin - 0.5;
   //    const auto height_diff = 12.5e-2;
   //    auto d = (x - scint_pos) * 17.e-3;
   //    auto theta = atan(d / height_diff);
   //    const auto rad_to_deg = 180 / M_PI;
   //    return theta * rad_to_deg;
   // };

   auto ce2a = [](Double_t x) -> Double_t {return cell_edges_to_angles(2 * x, left);};
   Double_t bin_edges[25];
   std::iota(bin_edges, bin_edges + 25, 0);
   std::transform(bin_edges, bin_edges + 25, bin_edges, ce2a);
   // std::transform(bin_edges, bin_edges + 25, bin_edges, local_cell_edges_to_angles);

   // TH1D angle_distribution = TH1D("basic_angle_distribution", "Winkelverteilung der Kosmischen Strahlung", 24, bin_edges);
   TH1D block_starts = TH1D("block_starts", "Winkelverteilung der Kosmischen Strahlung", WIRE_BINS);

   reset_entry_count();
   while(get_next_entry()) {
      // check for sequential cells that are hit
      checklist_64 wires_hit = checklist_64();
      forr (i, 0, nhits_le) {
         if (filter_exclude(i)) continue;
         else wires_hit.check(wire_le[i]);
      }

      auto sequences = get_sequences(wires_hit);
      for (auto seq: sequences) {
         int inner_most_wire = 100;
         for (auto wire: seq) {
            if (wire % 2) continue;
            if (abs_diff(wire, middle_bin) < abs_diff(inner_most_wire, middle_bin))
               inner_most_wire = wire;
         }

         block_starts.Fill(inner_most_wire);
      }
   }

   TH1D angle_distribution = TH1D("basic_angle_distribution", "Winkelverteilung der Kosmischen Strahlung", 24, bin_edges);
   forr (wire, 1, 49) {
      auto count = block_starts.GetBinContent(wire);
      int direction = wire < middle_bin? 1: -1;
      int target_bin = wire % 2? wire / 2 + direction: wire / 2;
      // if (wire % 2) continue;
      angle_distribution.AddBinContent(target_bin, count);
   }

   return angle_distribution;
}

class PathReconstruction {
   public:
      PathReconstruction(checklist_64 hit_wires, Double_t measured_dists[48]) {         
         forr (i, 0, 48) {
            Double_t angle_lb = std::min(cell_edges_to_angles(i, left), cell_edges_to_angles(i, right));
            Double_t angle_ub = std::max(cell_edges_to_angles(i, left), cell_edges_to_angles(i, right));
            angle_intervals[i] = {angle_lb, angle_ub};
            base_theta[i] = cell_edges_to_angles(i, center);
         }

         // copy just to make this technically thread safe. who tf know when root decides to
         // multithread implicitely at this point
         // (also costs like... nothing)
         forr (i, 0, 48) this->dists[i] = measured_dists[i];
         this->hits = hit_wires;
      };

      Double_t get_angle(wire_chunk_t);

   private:
      const static auto middle_bin = SCINT_CENTERED_OVER;
      static std::tuple<Double_t, Double_t> angle_intervals[48];
      static Double_t base_theta[48];

      Double_t dists[48];
      checklist_64 hits;
};


Double_t PathReconstruction::get_angle(wire_chunk_t wires) {
   auto sequences = get_sequences(hits);

   auto inner_most_even = 50;
   for (auto wire : wires) {
      if (abs_diff(wire, middle_bin) < abs_diff(inner_most_even, middle_bin) and wire IS_EVEN) {
         inner_most_even = wire;
      }
   }

   // if the particle hit inner_most_even, it mus fall within this angle interval
   auto [angle_from, angle_to] = angle_intervals[inner_most_even];

   for (auto wire : wires) {
      ;
   }

   return 0.;
}


std::tuple<Double_t, Double_t> get_angle(std::vector<unsigned int> wires) {
   const auto middle_bin = 20;
   const auto height_top_row = 12.5e-2;
   const auto height_bottom_row = height_top_row + 1.7e-2; // TODO: check!!!
   const auto scint_edge_l = middle_bin * 8.5e-3;
   const auto scint_edge_r = (middle_bin + 1) * 8.5e-3;
   const auto rad_to_deg = 180 / M_PI;

   std::tuple<Double_t, Double_t> angle_intervals[48];
   const auto cell_edges_to_angles = [=](unsigned int n, Double_t scint_pos) -> Double_t {
      if (n % 2) {
         auto d = (n * 8.5e-3) - scint_pos;
         auto theta = atan(d / height_bottom_row);
         return theta * rad_to_deg;
      } else {
         auto d = (n * 8.5e-3) - scint_pos;
         auto theta = atan(d / height_top_row);
         return theta * rad_to_deg;
      }
   };

   forr (i, 0, 48) {
      Double_t angle_lb = std::min(cell_edges_to_angles(i, scint_edge_l), cell_edges_to_angles(i, scint_edge_r));
      Double_t angle_ub = std::max(cell_edges_to_angles(i, scint_edge_l), cell_edges_to_angles(i, scint_edge_r));
      angle_intervals[i] = {angle_lb, angle_ub};
      // printf("cell %d: %lf - %lf\n", i, angle_lb, angle_ub);
   }

   auto inner_most_even = 50;
   for (auto wire : wires) {
      if (abs_diff(wire, middle_bin) < abs_diff(inner_most_even, middle_bin) and wire IS_EVEN) {
         inner_most_even = wire;
      }
   }

   // we always assume the innermost top level wire to be an actual particle
   // that passed through the scint
   // so if anything afterwards doesn't fit
   //chuck it !!!! 
   auto [valid_from, valid_to] = angle_intervals[inner_most_even];

   Double_t lb = valid_from;
   Double_t ub = valid_to;
   for (auto wire : wires) {
      auto [wire_from, wire_to] = angle_intervals[wire];

      // check if angle interval doesn't overlap with the valid region
      if (wire_from > valid_to) continue;
      if (wire_to < valid_from) continue;
      if (wire_from > lb && wire_from < ub) lb = wire_from;
      if (wire_to < ub && wire_to > lb) ub = wire_to;
   }
   return {lb, ub};
}


TH1D analysis::precise_angle_distribution() {
   TH1D angle_distribution = TH1D("precise_angle_distribution", "Winkelverteilung der Kosmischen Strahlung", 45, -60., 70.);

   reset_entry_count();
   Double_t dists[48];
   while(get_next_entry()) {
      // check for sequential cells that are hit
      checklist_64 wires_hit = checklist_64();
      forr (i, 0, nhits_le) {
         if (filter_exclude(i)) continue;
         else {
            wires_hit.check(wire_le[i]);
            dists[wire_le[i]] = this->dt_lut[time_le[i]];
         }
      }

      auto sequences = get_sequences(wires_hit);
      for (auto seq: sequences) {
         auto [angle_lb, angle_ub] = get_angle(seq);
         auto avg_angle = 0.5 * (angle_lb + angle_ub);
         angle_distribution.Fill(avg_angle);
         // printf("block from: %u %u\n", block_start, block_end);
         // printf("event_time: %lf\n", eventTime - 1778000000);
         // printf("angle_interval in (%u): %lf %lf\n", event, angle_lb, angle_ub);
      }

   }
   return angle_distribution;
}

void analysis::print_all_events() {
   reset_entry_count();
   while (get_next_entry()) {
      checklist_64 wires_hit = checklist_64();
      printf("event (%u) at: %lf\n", event, eventTime);
      printf("event contains %u hits\n", nhits_le);
      printf("wire hits are: \n\t");
      forr (i, 0, nhits_le) {
         printf("%u, ", wire_le[i]);
         if (filter_exclude(i)) continue;
         else wires_hit.check(wire_le[i]);
      }
      // printf("\n\t");
      // forr (i, 0, n_valid) {
      //    printf("%u, ", valid[i]);
      // }
      printf("\n");

      auto sequences = get_sequences(wires_hit);
      if (sequences.size() > 0) {
         printf("sequences: ");
         forr(i, 0, 48) printf("%d", wires_hit.check_static(i));
         printf("\n");
         for (auto& s : sequences) {
            for (auto& n : s) {
               printf("%u ", n);
            }
            printf("\n");
         }
         printf("\n");
      }
   }
}


Double_t analysis::get_runtime() {
   Double_t min_t = 1e11;
   Double_t max_t = 0;
   reset_entry_count();
   while(get_next_entry()) {
      if (eventTime < min_t) {
         min_t = eventTime;
      }

      if (eventTime > max_t) {
         max_t = eventTime;
      }
   }

   return max_t - min_t;
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
                  // printf("wires: %u %u    le times: %u %u    tot times: %u %u   n+-2 hits: %d %d\n", wire_le[hit], wire_le[i], time_le[hit], time_le[i], tot[hit], tot[i], hit_wires.check_static(wire_le[hit] - 2), hit_wires.check_static(wire_le[hit] + 2));
               // break;
            }
         }
      }

   }
   return dists;
}


void plot(TFitResult& p, TCanvas* C) {
   C->cd();
   p.Draw();
   C->Update();
}


template<typename Plotable>
void plot(Plotable& p, TCanvas* C) {
   C->cd();
   p.SetStats(0);
   p.Draw();
   C->Update();
}


void plot_set(std::vector<dataset*>& files, TCanvas* canv, global_object_store* gob, std::vector<TString> names, int j = 0) {
   canv->cd();
   EColor colors[] = {kBlue, kPink, kRed, kOrange, kAzure, kCyan, kMagenta, kTeal};
   int i = 0;
   for (dataset* ds : files) {
      ds->ana->filter_enabled = false;
      Double_t runtime = ds->ana->get_runtime();
      auto hist = gob->hold(ds->ana->dt_hist(names[i]));
      hist->Scale(1 / runtime);
      hist->SetTitle("");
      hist->SetStats(false);
      hist->SetXTitle("Driftzeit / ns");
      hist->SetYTitle("Rate / s");

      hist->SetLineColor(colors[i]);
      if (i == 0) hist->DrawCopy("HIST", "");
      else hist->DrawCopy("HIST same", "");

      i++;
   }
   canv->BuildLegend();
   canv->Update();
}


void plot_parameter_search() {
   TCanvas* voltage_canvas = new TCanvas("voltage_sweep", "Verschiedene Beschleunigungsspannungen");
   std::vector<dataset*> voltage_data = std::vector<dataset*> ({
      new dataset("../data/B103/run_260513_145435.root"),
      new dataset("../data/B103/run_260513_145148.root"),
      new dataset("../data/B103/run_260513_144851.root"),
      new dataset("../data/B103/run_260513_150659.root"),
      new dataset("../data/B103/run_260513_150905.root"),
      new dataset("../data/B103/run_260513_153904.root")
   });

   std::vector<TString> names = {
      "2.599kV",
      "2.499kV",
      "2.403kV",
      "2.451kV",
      "2.482kV",
      "2.801kV"
   };

   global_object_store* gob = new global_object_store();
   plot_set(voltage_data, voltage_canvas, gob, names);

   TCanvas* discriminator_canvas = new TCanvas("discriminator_sweep", "Verschiedene Diskriminatoreinstellungen");
   std::vector<dataset*> discriminator_data = std::vector<dataset*> ({
      new dataset("../data/B103/run_260513_160116.root"),
      new dataset("../data/B103/run_260513_155621.root"),
      new dataset("../data/B103/run_260513_155230.root"),
      new dataset("../data/B103/run_260513_154351.root"),
      new dataset("../data/B103/run_260513_153904.root"),
   });

   std::vector<TString> disc_names = {
      "0x40",
      "0x68",
      "0x58",
      "0x28",
      "0x20"
   };
   plot_set(discriminator_data, discriminator_canvas, gob, disc_names, 1);
}


int main(int argc, char** argv) {
   // weniger statistik mit guten parametern ..161632.root
   // viel statistik mit guten parametern ...161826.root

   TROOT root("app","app");
   Int_t dargc=1;
   char** dargv = &argv[0];
   TRint app = TRint("app", &dargc, dargv);
   std::vector<TCanvas*> canvas_vec = std::vector<TCanvas*>();
   forr(i, 0, 8) canvas_vec.push_back(new TCanvas(("c" + std::to_string(i)).c_str(), "c", 800, 600));


   dataset* data = new dataset(argv[1]);
   analysis* ana = data->ana;
   ana->wire_lut = make_wire_lut();

   dataset* calib_dataset = new dataset("../data/B103/run_260513_161826.root");
   calib_dataset->ana->max_le_time = 500;

   // auto dt_rel = ana->dt_hist();
   auto dt_rel = calib_dataset->ana->dt_hist();
   auto odb = make_odb(dt_rel);
   ana->dt_lut = odb;
   auto dt_wire_plot = ana->dt_wire_hist();
   auto wire_correlation = ana->wire_correlation();
   auto dt_tot = ana->dt_tot_relation();
   // TODO: maybe this needs to be done before filtering??
   // auto sum_vs_diff = ana->dist_plot(odb, 16, 24);
   auto sum_vs_diff = ana->dist_plot(odb, 0, 50);
   // auto sum_vs_diff = ana->dist_plot(odb, 18, 22);

   TF1 angle_dist_func = TF1("angle_dist_func", "[0] * cos((x - [2]) * pi / 180)^[1]", -90, 90);
   angle_dist_func.SetParameter(0, 700);
   angle_dist_func.SetParameter(1, 2);

   auto basic_angles = ana->basic_angle_distrib();
   TFitResultPtr fit = basic_angles.Fit(&angle_dist_func, "S");

   auto precise_angles = ana->precise_angle_distribution();

   plot(dt_rel, canvas_vec[0]);
   plot(dt_wire_plot, canvas_vec[1]);
   plot(wire_correlation, canvas_vec[2]);
   plot(dt_tot, canvas_vec[3]);
   plot(odb, canvas_vec[4]);
   plot(basic_angles, canvas_vec[5]);
   // plot(*fit, canvas_vec[5]);
   plot(sum_vs_diff, canvas_vec[6]);
   plot(precise_angles, canvas_vec[7]);

   // plot_parameter_search();

   // ana->print_all_events();

   app.Run(kTRUE);
}
