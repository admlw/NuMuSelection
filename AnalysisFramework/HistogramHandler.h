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
        "muoncand_track_mcs_bwd"
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
        ";Muon candidate track MCS (backward;",
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
        {25, 0, 3}            // muon candidate mcs bwd
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
