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
        "track_mcs_fwd",
        "track_mcs_bwd",
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
        ";Track MCS (forward);",
        ";Track MCS (backward);",
        ";Muon candidate track length (cm);",
        ";Muon candidate track theta (rad);",
        ";Muon candidate track cos(theta);",
        ";Muon candidate track phi;",
        ";Muon candidate track MCS (forward);",
        ";Muon candidate track MCS (backward;",
      };

      std::vector<std::vector<double>> histoBins = {
        {1, 0, 1},
        {10, 0, 10},
        {10, 0, 10},
        {50, 0, 700},
        {50, 0, 256},
        {50, -116.5, 116.5},
        {50, 0, 1040},
        {50, 0, 256},
        {50, 0, 256},
        {50, -116.5, 116.5},
        {50, -116.5, 116.5},
        {50, 0, 1040},
        {50, 0, 1040},
        {25, 0, 3.15},
        {25, -1, 1},
        {50, -3.15, 3.15},
        {50, -10, 10},
        {50, 0, 3},
        {50, 0, 3},
        {25, 0, 700},
        {25, 0, 3.15},
        {25, -1, 1},
        {25, -3.15, 3.15},
        {25, 0, 3},
        {25, 0, 3}
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