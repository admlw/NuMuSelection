#ifndef HISTOGRAMHANDLER_H
#define HISTOGRAMHANDLER_H

// cpp
#include <iostream>
#include <bitset>
#include <string>
#include <vector>

// ROOT
#include "TF1.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TLegend.h"
#include "TStyle.h"

// local
#include "DataTypes.h"
#include "Configuration.h"

namespace numusel{

  class HistogramHandler{

    private:
      numusel::Configuration _config;

    public:

      std::vector<std::string> histoNames_2D = {
        "h_dedx_resrg_mucand",
        "h_dedx_resrg_leadingpcand",
        "h_dedx_resrg_nonleadingpcand"
      };

      std::vector<std::string> histoLabels_2D = {
        ";Muon Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        ";Leading Proton Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        ";Non-Leading Proton Candidate;Residual Range (cm);dE/dx (MeV/cm)"
      };
     
      std::vector<std::vector<double>> histoBins_2D = {
        {50, 0, 20, 50, 0, 20},
        {50, 0, 20, 50, 0, 20},
        {50, 0, 20, 50, 0, 20}
      };

      std::vector<std::string> histoNames = {
        "total_purity",
        "nTracks",
        "nShowers",
        "nSelectedPfParticles",
        "trackLength",
        "vertex_x",
        "vertex_y",
        "vertex_z",
        "track_startx",
        "track_endx",
        "track_starty",
        "track_endy",
        "track_startz",
        "track_endz",
        "track_theta",
        "track_costheta",
        "track_phi",
        "log(lmipoverp)",
        "track_mcs_muassmp_fwd",
        "track_mcs_muassmp_bwd",
        "track_mcs_muassmp_energy_fwd",
        "track_mcs_muassmp_energy_bwd",
        "track_range_mom_muassumption",
        "track_range_mom_passumption",
        "track_range_energy_muassumption",
        "track_range_energy_passumption",
        "muoncand_track_length",
        "muoncand_track_theta",
        "muoncand_track_costheta",
        "muoncand_track_phi",
        "muoncand_track_mcs_fwd",
        "muoncand_track_mcs_bwd",
        "muoncand_track_mcs_energy_fwd",
        "muoncand_track_mcs_energy_bwd",
        "muoncand_track_range_momentum_muassmp",
        "muoncand_track_range_energy_muassmp",
        "muoncand_track_contained_length",
        "muoncand_track_contained_theta",
        "muoncand_track_contained_costheta",
        "muoncand_track_contained_phi",
        "muoncand_track_contained_mcs_fwd",
        "muoncand_track_contained_mcs_bwd",
        "muoncand_track_contained_mcs_energy_fwd",
        "muoncand_track_contained_mcs_energy_bwd",
        "muoncand_track_contained_range_momentum_muassmp",
        "muoncand_track_contained_range_energy_muassmp",
        "muoncand_track_uncontained_length",
        "muoncand_track_uncontained_theta",
        "muoncand_track_uncontained_costheta",
        "muoncand_track_uncontained_phi",
        "muoncand_track_uncontained_mcs_fwd",
        "muoncand_track_uncontained_mcs_bwd",
        "muoncand_track_uncontained_mcs_energy_fwd",
        "muoncand_track_uncontained_mcs_energy_bwd",
        "muoncand_track_uncontained_range_momentum_muassmp",
        "muoncand_track_uncontained_range_energy_muassmp",
        "protoncand_track_length",
        "protoncand_track_theta",
        "protoncand_track_costheta",
        "protoncand_track_phi",
        "protoncand_track_mcs_fwd",
        "protoncand_track_mcs_bwd",
        "protoncand_track_range_mom",
        "protoncand_track_range_energy",
        "protoncand_leading_track_length",
        "protoncand_leading_track_theta",
        "protoncand_leading_track_costheta",
        "protoncand_leading_track_phi",
        "protoncand_leading_track_mcs_fwd",
        "protoncand_leading_track_mcs_bwd",
        "protoncand_leading_track_range_mom",
        "protoncand_leading_track_range_energy",
        "protoncand_nonleading_track_length",
        "protoncand_nonleading_track_theta",
        "protoncand_nonleading_track_costheta",
        "protoncand_nonleading_track_phi",
        "protoncand_nonleading_track_mcs_fwd",
        "protoncand_nonleading_track_mcs_bwd",
        "protoncand_nonleading_track_range_mom",
        "protoncand_nonleading_track_range_energy"
      };

      std::vector<std::string> histoLabels = {
        ";Total purity;",
        ";Number of tracks;",
        ";Number of showers;",
        ";Number of PFParticles;",
        ";Track length (cm);",
        ";Vertex x position (cm);",
        ";Vertex y position (cm);",
        ";Vertex z position (cm);",
        ";Track start x (cm);",
        ";Track end x (cm);",
        ";Track start y (cm);",
        ";Track end y (cm);",
        ";Track start z (cm);",
        ";Track end z (cm);",
        ";Track theta (rad);",
        ";Track cos(theta);",
        ";Track phi;",
        ";Log(L_{MIP}/L_{p});",
        ";Track MCS (forward, muon assumption);",
        ";Track MCS (backward, muon assumption);",
        ";Track MCS Energy (forward, muon assumption);",
        ";Track MCS Energy (backward, muon assumption);",
        ";Track Range Momentum (muon assumption);",
        ";Track Range Momentum (proton assumption);",
        ";Track Range Energy (muon assumption);",
        ";Track Range Energy (proton assumption);",
        ";Muon candidate track length (cm);",
        ";Muon candidate track theta (rad);",
        ";Muon candidate track cos(theta);",
        ";Muon candidate track phi;",
        ";Muon candidate track MCS (forward);",
        ";Muon candidate track MCS (backward);",
        ";Muon candidate track MCS Energy (forward);",
        ";Muon candidate track MCS Energy (backward);",
        ";Muon candidate Momentum (range, muon assumption);",
        ";Muon candidate Energy (range, muon assumption);",
        ";Contained Muon candidate track length (cm);",
        ";Contained Muon candidate track theta (rad);",
        ";Contained Muon candidate track cos(theta);",
        ";Contained Muon candidate track phi;",
        ";Contained Muon candidate track MCS (forward);",
        ";Contained Muon candidate track MCS (backward);",
        ";Contained Muon candidate track MCS Energy (forward);",
        ";Contained Muon candidate track MCS Energy (backward);",
        ";Contained Muon candidate Momentum (range, muon assumption);",
        ";Contained Muon candidate Energy (range, muon assumption);",
        ";Uncontained Muon candidate track length (cm);",
        ";Uncontained Muon candidate track theta (rad);",
        ";Uncontained Muon candidate track cos(theta);",
        ";Uncontained Muon candidate track phi;",
        ";Uncontained Muon candidate track MCS (forward);",
        ";Uncontained Muon candidate track MCS (backward);",
        ";Uncontained Muon candidate track MCS Energy (forward);",
        ";Uncontained Muon candidate track MCS Energy (backward);",
        ";Uncontained Muon candidate Momentum (range, muon assumption);",
        ";Uncontained Muon candidate Energy (range, muon assumption);",
        ";Proton candidate track length (cm);",
        ";Proton candidate track theta (rad);",
        ";Proton candidate track cos(theta);",
        ";Proton candidate track phi;",
        ";Proton candidate track MCS (forward);",
        ";Proton candidate track MCS (backward);",
        ";Proton candidate track momentum (range);",
        ";Proton candidate track energy (range);",
        ";Leading Proton candidate track length (cm);",
        ";Leading Proton candidate track theta (rad);",
        ";Leading Proton candidate track cos(theta);",
        ";Leading Proton candidate track phi;",
        ";Leading Proton candidate track MCS (forward);",
        ";Leading Proton candidate track MCS (backward;",
        ";Leading Proton candidate track momentum (range);",
        ";Leading Proton candidate track energy (range);",
        ";Non-Leading Proton candidate track length (cm);",
        ";Non-Leading Proton candidate track theta (rad);",
        ";Non-Leading Proton candidate track cos(theta);",
        ";Non-Leading Proton candidate track phi;",
        ";Non-Leading Proton candidate track MCS (forward);",
        ";Non-Leading Proton candidate track MCS (backward);",
        ";Non-Leading Proton candidate track momentum (range);",
        ";Non-Leading Proton candidate track energy (range);"
      };

      std::vector<std::vector<double>> histoBins = {
        {1, 0, 1},            // total purity 
        {10, 0, 10},          // number of tracks
        {10, 0, 10},          // number of showers
        {10, 0, 10},          // number of pfparticles
        {50, 0, 700},         // track length
        {50, 0, 256},         // vertex x
        {50, -116.5, 116.5},  // vertex y
        {50, 0, 1040},        // vertex z
        {50, 0, 256},         // track start x
        {50, 0, 256},         // track end x
        {50, -116.5, 116.5},  // track start y
        {50, -116.5, 116.5},  // track end y
        {50, 0, 1040},        // track start z
        {50, 0, 1040},        // track end z
        {25, 0, 3.15},        // track theta
        {25, -1, 1},          // track cos(theta)
        {50, -3.15, 3.15},    // track phi
        {50, -10, 10},        // log(lmip/lp)
        {50, 0, 3},           // track mcs fwd
        {50, 0, 3},           // track mcs bwd
        {50, 0, 0.000015},    // track mcs energy fwd
        {50, 0, 0.000015},    // track mcs energy bwd
        {50, 0, 3},           // track range momentum muon
        {50, 0, 3},           // track range momentum proton
        {50, 0, 3},           // track range energy muon
        {50, 0, 3},           // track range energy proton
        {25, 0, 700},         // muon cand length
        {25, 0, 3.15},        // muon cand theta
        {25, -1, 1},          // muon cand cos(theta)
        {25, -3.15, 3.15},    // muon cand phi
        {25, 0, 3},           // muon candidate mcs fwd
        {25, 0, 3},           // muon candidate mcs bwd
        {25, 0, 2},           // muon candidate mc energy forward
        {25, 0, 2},           // muon candidate mc energy backward
        {25, 0, 3},           // muon candidate momentum from range
        {25, 0, 2},           // muon candidate energy from range
        {25, 0, 700},         // contained muon cand length
        {25, 0, 3.15},        // contained muon cand theta
        {25, -1, 1},          // contained muon cand cos(theta)
        {25, -3.15, 3.15},    // contained muon cand phi
        {25, 0, 3},           // contained muon candidate mcs fwd
        {25, 0, 3},           // contained muon candidate mcs bwd
        {25, 0, 2},           // contained muon candidate mc energy forward
        {25, 0, 2},           // contained muon candidate mc energy backward
        {25, 0, 3},           // contained muon candidate momentum from range
        {25, 0, 2},           // contained muon candidate energy from range
        {25, 0, 700},         // uncontained muon cand length
        {25, 0, 3.15},        // uncontained muon cand theta
        {25, -1, 1},          // uncontained muon cand cos(theta)
        {25, -3.15, 3.15},    // uncontained muon cand phi
        {25, 0, 3},           // uncontained muon candidate mcs fwd
        {25, 0, 3},           // uncontained muon candidate mcs bwd
        {25, 0, 2},           // uncontained muon candidate mc energy forward
        {25, 0, 2},           // uncontained muon candidate mc energy backward
        {25, 0, 3},           // uncontained muon candidate momentum from range
        {25, 0, 2},           // uncontained muon candidate energy from range
        {25, 0, 300},         // proton cand length
        {25, 0, 3.15},        // proton cand theta
        {25, -1, 1},          // proton cand cos(theta)
        {25, -3.15, 3.15},    // proton cand phi
        {25, 0, 3},           // proton candidate mcs fwd
        {25, 0, 3},           // proton candidate mcs bwd
        {25, 0, 3.0},         // proton candidate range momentum
        {25, 0, 0.7},         // proton candidate range energy
        {25, 0, 300},         // leading proton cand length
        {25, 0, 3.15},        // leading proton cand theta
        {25, -1, 1},          // leading proton cand cos(theta)
        {25, -3.15, 3.15},    // leading proton cand phi
        {25, 0, 3},           // leading proton candidate mcs fwd
        {25, 0, 3},           // leading proton candidate mcs bwd
        {25, 0, 3.0},         // leading proton candidate range momentum
        {25, 0, 0.7},         // leading proton candidate range energy
        {25, 0, 300},         // non-leading proton cand length
        {25, 0, 3.15},        // non-leading proton cand theta
        {25, -1, 1},          // non-leading proton cand cos(theta)
        {25, -3.15, 3.15},    // non-leading proton cand phi
        {25, 0, 3},           // non-leading proton candidate mcs fwd
        {25, 0, 3},           // non-leading proton candidate mcs bwd
        {25, 0, 3.0},         // non-leading proton candidate range momentum
        {25, 0, 0.7}          // non-leading proton candidate range energy
      };

      /*
       * These three vectors fully define all of the efficiency/purity histograms to produce
       */
      std::vector<std::string> effpurNames = {
        "totel_eff",
        "true_enu",
        "true_mu_p"
      };

      std::vector<std::string> effpurLabels = {
        ";Total Efficiency;",
        ";E_{#nu}^{true} (Gev);",
        ";p_{#mu} (Gev);"
      };

      std::vector< std::vector<double> > effpurBins = {
        {1, 0, 1},
        {25, 0, 3},
        {25, 0, 2.5}
      };

      /**
       * Fill 2D histograms
       */
      void Fill2DHist(TH2D* h2d, std::vector<double> variable_x, std::vector<double> variable_y);

      /**
       * Fill MC histograms
       */
      void FillHistMC(hists_1d* h1d, std::vector<double> variable, std::bitset<8> eventCat);

      /**
       * Fill On-beam histograms
       */
      void FillHistOnBeam(hists_1d* h1d, std::vector<double> variable);

      /**
       * Fill Off-beam histograms
       */
      void FillHistOffBeam(hists_1d* h1d, std::vector<double> variable);

      /**
       * Initialise vector of hists_1ds
       */
      void InitialiseHistoVec(std::vector< std::vector<hists_1d*> >* plots_to_make, int n_plots);

      /**
       * Initialise vector of TH2s
       */
      void InitialiseHistoVec(std::vector< std::vector<TH2D*> >* plots_to_make, int n_plots);

      /**
       * Initialise vector of eff_1ds
       */
      void InitialiseHistoVec(std::vector< std::vector<eff_1d*> >* plots_to_make, int n_plots);

      /**
       * Style histograms
       */
      void StyleHistograms(hists_1d* hists);

      /**
       * Scale histograms
       */
      void ScaleHistograms(std::vector< std::vector<hists_1d*> > hists);

      /**
       * Make stacked data/mc histograms and save to png/pdf file
       */
      void MakeStackedHistogramsAndSave(std::vector< std::vector<hists_1d*> > hists);

      /**
       * Make efficiency histograms and save to png/pdf files
       */
      void MakeEfficiencyHistogramsAndSave(std::vector< std::vector<eff_1d*> > effplots);

  };

}

#endif
