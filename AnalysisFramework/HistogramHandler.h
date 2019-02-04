#ifndef HISTOGRAMHANDLER_H
#define HISTOGRAMHANDLER_H

// cpp
#include <iostream>
#include <bitset>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

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
#include "TPaveText.h"
#include "TFile.h"

// local
#include "DataTypes.h"
#include "Configuration.h"

namespace numusel{

  class HistogramHandler{

    private:
      numusel::Configuration _config;

    public:

      std::vector<std::string> histoNames_2D = {
        "h_track_theta_phi",
        "h_dedx_resrg_mucand",
        "h_dedx_resrg_mucand_contained",
        "h_dedx_resrg_mucand_uncontained",
        "h_dedx_resrg_pcand",
        "h_dedx_resrg_leadingpcand",
        "h_dedx_resrg_nonleadingpcand",
        "h_trueenu_recoenu"
      };

      std::vector<std::string> histoLabels_2D = {
        "Track Theta; Track Phi",
        "Muon Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        "Contained Muon Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        "Uncontained Muon Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        "Proton Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        "Leading Proton Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        "Non-Leading Proton Candidate;Residual Range (cm);dE/dx (MeV/cm)",
        ";true E_{#nu};Reco E_{#nu}"
      };
     
      std::vector<std::vector<double>> histoBins_2D = {
        {50, 0, 3.15, 50, -3.15, 3.15},
        {50, 0, 30, 50, 0, 20},
        {50, 0, 30, 50, 0, 20},
        {50, 0, 30, 50, 0, 20},
        {50, 0, 30, 50, 0, 20},
        {50, 0, 30, 50, 0, 20},
        {50, 0, 30, 50, 0, 20},
        {50, 0, 3, 50, 0, 3}
      };

      std::vector<std::string> histoNames = {
        "total_purity",
        "nTracks",
        "nShowers",
        "nSelectedPfParticles",
        "vertex_x",
        "vertex_y",
        "vertex_z",
        "trackLength",
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
        "track_dep_energy_yplane",
        "track_dep_energy_yplane_minus_range_energy_muassmp",
        "track_dep_energy_yplane_minus_range_energy_passmp",
        "track_residualrms",
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
        "muoncand_track_dep_energy_yplane",
        "muoncand_track_dep_energy_yplane_minus_range_energy_muassmp",
        "muoncand_track_dep_energy_yplane_minus_range_energy_passmp",
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
        "muoncand_track_contained_dep_energy_yplane",
        "muoncand_track_contained_dep_energy_yplane_minus_range_energy_muassmp",
        "muoncand_track_contained_dep_energy_yplane_minus_range_energy_passmp",
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
        "muoncand_track_uncontained_dep_energy_yplane",
        "muoncand_track_uncontained_dep_energy_yplane_minus_range_energy_muassmp",
        "muoncand_track_uncontained_dep_energy_yplane_minus_range_energy_passmp",
        "protoncand_track_length",
        "protoncand_track_theta",
        "protoncand_track_costheta",
        "protoncand_track_phi",
        "protoncand_track_mcs_fwd",
        "protoncand_track_mcs_bwd",
        "protoncand_track_range_mom",
        "protoncand_track_range_energy",
        "protoncand_track_dep_energy_yplane",
        "protoncand_track_dep_energy_yplane_minus_range_energy_muassmp",
        "protoncand_track_dep_energy_yplane_minus_range_energy_passmp",
        "protoncand_leading_track_length",
        "protoncand_leading_track_theta",
        "protoncand_leading_track_costheta",
        "protoncand_leading_track_phi",
        "protoncand_leading_track_mcs_fwd",
        "protoncand_leading_track_mcs_bwd",
        "protoncand_leading_track_range_mom",
        "protoncand_leading_track_range_energy",
        "protoncand_leading_track_dep_energy_yplane",
        "protoncand_leading_track_dep_energy_yplane_minus_range_energy_muassmp",
        "protoncand_leading_track_dep_energy_yplane_minus_range_energy_passmp",
        "protoncand_nonleading_track_length",
        "protoncand_nonleading_track_theta",
        "protoncand_nonleading_track_costheta",
        "protoncand_nonleading_track_phi",
        "protoncand_nonleading_track_mcs_fwd",
        "protoncand_nonleading_track_mcs_bwd",
        "protoncand_nonleading_track_range_mom",
        "protoncand_nonleading_track_range_energy",
        "protoncand_nonleading_track_dep_energy_yplane",
        "protoncand_nonleading_track_dep_energy_yplane_minus_range_energy_muassmp",
        "protoncand_nonleading_track_dep_energy_yplane_minus_range_energy_passmp",
        "reconstructed_neutrino_energy_uncalib",
        "reconstructed_neutrino_energy_calib"
      };

      std::vector<std::string> histoLabels = {
        ";Total purity;",
        ";Number of tracks;",
        ";Number of showers;",
        ";Number of PFParticles;",
        ";Vertex x position (cm);",
        ";Vertex y position (cm);",
        ";Vertex z position (cm);",
        ";Track length (cm);",
        ";Track start x (cm);",
        ";Track end x (cm);",
        ";Track start y (cm);",
        ";Track end y (cm);",
        ";Track start z (cm);",
        ";Track end z (cm);",
        ";Track theta (rad);",
        ";Track cos(#theta);",
        ";Track #phi;",
        ";Log(L_{MIP}/L_{p});",
        ";Track MCS (GeV, forward, muon assumption);",
        ";Track MCS (GeV, backward, muon assumption);",
        ";Track MCS energy (GeV, forward, muon assumption);",
        ";Track MCS energy (GeV, backward, muon assumption);",
        ";Track range momentum (GeV, muon assumption);",
        ";Track range momentum (GeV, proton assumption);",
        ";Track range energy (GeV, muon assumption);",
        ";Track range energy (GeV, proton assumption);",
        ";Track deposited energy (GeV, Y plane);",
        ";Track E_{dep}^{Y}-E_{range}^{#mu};",
        ";Track E_{dep}^{Y}-E_{range}^{p};",
        ";Track residual (cm);",
        ";Muon candidate track length (cm);",
        ";Muon candidate track #theta (rad);",
        ";Muon candidate track cos(#theta);",
        ";Muon candidate track #phi;",
        ";Muon candidate track MCS (GeV, forward);",
        ";Muon candidate track MCS (GeV, backward);",
        ";Muon candidate track MCS energy (GeV, forward);",
        ";Muon candidate track MCS energy (GeV, backward);",
        ";Muon candidate momentum (GeV, range, muon assumption);",
        ";Muon candidate energy (GeV, range, muon assumption);",
        ";Muon candidate deposited energy (GeV, Y plane);",
        ";Muon candidate E_{dep}^{Y}-E_{range}^{#mu};",
        ";Muon candidate E_{dep}^{Y}-E_{range}^{p};",
        ";Contained muon candidate track length (cm);",
        ";Contained muon candidate track #theta (rad);",
        ";Contained muon candidate track cos(#theta);",
        ";Contained muon candidate track #phi;",
        ";Contained muon candidate track MCS (GeV, forward);",
        ";Contained muon candidate track MCS (GeV, backward);",
        ";Contained muon candidate track MCS energy (GeV, forward);",
        ";Contained muon candidate track MCS energy (GeV, backward);",
        ";Contained muon candidate momentum (GeV, range, muon assumption);",
        ";Contained muon candidate energy (GeV, range, muon assumption);",
        ";Contained muon candidate deposited energy (GeV, Y plane);",
        ";Contained muon candidate E_{dep}^{Y}-E_{range}^{#mu};",
        ";Contained muon candidate E_{dep}^{Y}-E_{range}^{p};",
        ";Uncontained muon candidate track length (cm);",
        ";Uncontained muon candidate track #theta (rad);",
        ";Uncontained muon candidate track cos(#theta);",
        ";Uncontained muon candidate track #phi;",
        ";Uncontained muon candidate track MCS (GeV, forward);",
        ";Uncontained muon candidate track MCS (GeV, backward);",
        ";Uncontained muon candidate track MCS energy (forward);",
        ";Uncontained muon candidate track MCS energy (GeV, backward);",
        ";Uncontained muon candidate momentum (GeV, range, muon assumption);",
        ";Uncontained muon candidate energy (GeV, range, muon assumption);",
        ";Uncontained muon candidate deposited energy (GeV, Y plane);",
        ";Uncontained muon candidate E_{dep}^{Y}-E_{range}^{#mu};",
        ";Uncontained muon candidate E_{dep}^{Y}-E_{range}^{p};",
        ";Proton candidate track length (cm);",
        ";Proton candidate track #theta (rad);",
        ";Proton candidate track cos(#theta);",
        ";Proton candidate track #phi;",
        ";Proton candidate track MCS (GeV, forward);",
        ";Proton candidate track MCS (GeV, backward);",
        ";Proton candidate track momentum (GeV, range);",
        ";Proton candidate track energy (GeV, range);",
        ";Proton candidate deposited energy (GeV, Y plane);",
        ";Proton candidate E_{dep}^{Y}-E_{range}^{#mu};",
        ";Proton candidate E_{dep}^{Y}-E_{range}^{p};",
        ";Leading proton candidate track length (cm);",
        ";Leading proton candidate track #theta (rad);",
        ";Leading proton candidate track cos(#theta);",
        ";Leading proton candidate track #phi;",
        ";Leading proton candidate track MCS (GeV, forward);",
        ";Leading proton candidate track MCS (GeV, backward);",
        ";Leading proton candidate track momentum (GeV, range);",
        ";Leading proton candidate track energy (GeV, range);",
        ";Leading proton candidate deposited energy (GeV, Y plane);",
        ";Leading proton candidate E_{dep}^{Y}-E_{range}^{#mu};",
        ";Leading proton candidate E_{dep}^{Y}-E_{range}^{p};",
        ";Non-leading proton candidate track length (cm);",
        ";Non-leading proton candidate track #theta (rad);",
        ";Non-leading proton candidate track cos(#theta);",
        ";Non-leading proton candidate track #phi;",
        ";Non-leading proton candidate track MCS (GeV, forward);",
        ";Non-leading proton candidate track MCS (GeV, backward);",
        ";Non-leading proton candidate track momentum (GeV, range);",
        ";Non-leading proton candidate track energy (GeV, range);",
        ";Non-leading proton candidate deposited energy (GeV, Y plane);",
        ";Non-leading proton candidate E_{dep}^{Y}-E_{range}^{#mu};",
        ";Non-leading proton candidate E_{dep}^{Y}-E_{range}^{p};",
        ";Total deposited energy (GeV);",
        ";Total deposited energy (GeV);"
      };

      std::vector<std::vector<double>> histoBins = {
        {1, 0, 1},            // total purity 
        {10, 0, 10},          // number of tracks
        {10, 0, 10},          // number of showers
        {10, 0, 10},          // number of pfparticles
        {40, 0, 260},         // vertex x
        {40, -116.5, 116.5},  // vertex y
        {40, 0, 1040},        // vertex z
        {50, 0, 700},         // track length
        {40, 0, 260},       // track start x
        {40, 0, 260},       // track end x
        {40, -116.5, 116.5},  // track start y
        {40, -116.5, 116.5},  // track end y
        {40, 0, 1040},        // track start z
        {40, 0, 1040},        // track end z
        {20, 0, 3.15},        // track theta
        {20, -1, 1},          // track cos(theta)
        {20, -3.15, 3.15},    // track phi
        {50, -10, 10},        // log(lmip/lp)
        {50, 0, 3},           // track mcs fwd
        {50, 0, 3},           // track mcs bwd
        {50, 0, 1.5},           // track mcs energy fwd
        {50, 0, 1.5},           // track mcs energy bwd
        {50, 0, 3},           // track range momentum muon
        {50, 0, 3},           // track range momentum proton
        {50, 0, 1.5},           // track range energy muon
        {50, 0, 1.5},           // track range energy proton
        {50, 0, 1},           // track deposited energy y plane
        {50, -0.4, 0.4},      // track dep energy yplane minus range energy muassmp
        {50, -0.4, 0.4},      // track dep energy yplane minus range energy passmp
        {50, 0, 5},           // track_residualrms
        {25, 0, 700},         // muon cand length
        {20, 0, 3.15},        // muon cand theta
        {20, -1, 1},          // muon cand cos(theta)
        {20, -3.15, 3.15},    // muon cand phi
        {25, 0, 3},           // muon candidate mcs fwd
        {25, 0, 3},           // muon candidate mcs bwd
        {25, 0, 1.5},           // muon candidate mc energy forward
        {25, 0, 1.5},           // muon candidate mc energy backward
        {25, 0, 3},           // muon candidate momentum from range
        {25, 0, 1.5},           // muon candidate energy from range
        {50, 0, 1},           // muon candidate deposited energy y plane
        {50, -0.4, 0.4},      // muon candidate dep energy yplane minus range energy muassmp
        {50, -0.4, 0.4},      // muon candidate dep energy yplane minus range energy passmp
        {20, 0, 700},         // contained muon cand length
        {10, 0, 3.15},        // contained muon cand theta
        {10, -1, 1},          // contained muon cand cos(theta)
        {10, -3.15, 3.15},    // contained muon cand phi
        {20, 0, 3},           // contained muon candidate mcs fwd
        {20, 0, 3},           // contained muon candidate mcs bwd
        {20, 0, 1.5},           // contained muon candidate mc energy forward
        {20, 0, 1.5},           // contained muon candidate mc energy backward
        {20, 0, 3},           // contained muon candidate momentum from range
        {20, 0, 1.5},           // contained muon candidate energy from range
        {50, 0, 1},           // contained muon candidate deposited energy y plane
        {20, -0.4, 0.4},      // contained muon candidate dep energy yplane minus range energy muassmp
        {20, -0.4, 0.4},      // contained muon candidate dep energy yplane minus range energy passmp
        {25, 0, 700},         // uncontained muon cand length
        {20, 0, 3.15},        // uncontained muon cand theta
        {20, -1, 1},          // uncontained muon cand cos(theta)
        {20, -3.15, 3.15},    // uncontained muon cand phi
        {25, 0, 3},           // uncontained muon candidate mcs fwd
        {25, 0, 3},           // uncontained muon candidate mcs bwd
        {25, 0, 1.5},           // uncontained muon candidate mc energy forward
        {25, 0, 1.5},           // uncontained muon candidate mc energy backward
        {25, 0, 3},           // uncontained muon candidate momentum from range
        {25, 0, 1.5},           // uncontained muon candidate energy from range
        {50, 0, 1},           // uncontained muon candidate deposited energy y plane
        {50, -0.4, 0.4},      // uncontained muon candidate dep energy yplane minus range energy muassmp
        {50, -0.4, 0.4},      // uncontained muon candidate dep energy yplane minus range energy passmp
        {25, 0, 150},         // proton cand length
        {20, 0, 3.15},        // proton cand theta
        {20, -1, 1},          // proton cand cos(theta)
        {20, -3.15, 3.15},    // proton cand phi
        {25, 0, 3},           // proton candidate mcs fwd
        {25, 0, 3},           // proton candidate mcs bwd
        {25, 0, 3.0},         // proton candidate range momentum
        {25, 0, 0.7},         // proton candidate range energy
        {50, 0, 1},           // proton candidate deposited energy y plane
        {50, -0.4, 0.4},      // proton candidate dep energy yplane minus range energy muassmp
        {50, -0.4, 0.4},      // proton candidate dep energy yplane minus range energy passmp
        {25, 0, 150},         // leading proton cand length
        {20, 0, 3.15},        // leading proton cand theta
        {20, -1, 1},          // leading proton cand cos(theta)
        {20, -3.15, 3.15},    // leading proton cand phi
        {25, 0, 3},           // leading proton candidate mcs fwd
        {25, 0, 3},           // leading proton candidate mcs bwd
        {25, 0, 3.0},         // leading proton candidate range momentum
        {25, 0, 0.7},         // leading proton candidate range energy
        {50, 0, 1},           // leading proton candidate deposited energy y plane
        {50, -0.4, 0.4},      // leading proton candidate dep energy yplane minus range energy muassmp
        {50, -0.4, 0.4},      // leading proton candidate dep energy yplane minus range energy passmp
        {25, 0, 150},         // non-leading proton cand length
        {10, 0, 3.15},        // non-leading proton cand theta
        {10, -1, 1},          // non-leading proton cand cos(theta)
        {10, -3.15, 3.15},    // non-leading proton cand phi
        {25, 0, 3},           // non-leading proton candidate mcs fwd
        {25, 0, 3},           // non-leading proton candidate mcs bwd
        {25, 0, 3.0},         // non-leading proton candidate range momentum
        {25, 0, 0.7},         // non-leading proton candidate range energy
        {50, 0, 1},           // non-leading proton candidate deposited energy y plane
        {50, -0.4, 0.4},      // non-leading proton candidate dep energy yplane minus range energy muassmp
        {50, -0.4, 0.4},      // non-leading proton candidate dep energy yplane minus range energy passmp
        {25, 0, 3},           // reconstructed neutrino energy
        {25, 0, 3}            // calibrated reconstructed neutrino energy                  
      };

      std::vector<bool> histoMakeTrackPlot = {
        false,                // total purity 
        false,                // number of tracks
        false,                // number of showers
        false,                // number of pfparticles
        false,                // vertex x
        false,                // vertex y
        false,                // vertex z
        true,                 // track length
        true,                 // track start x
        true,                 // track end x
        true,                 // track start y
        true,                 // track end y
        true,                 // track start z
        true,                 // track end z
        true,                 // track theta
        true,                 // track cos(theta)
        true,                 // track phi
        true,                 // log(lmip/lp)
        true,                 // track mcs fwd
        true,                 // track mcs bwd
        true,                 // track mcs energy fwd
        true,                 // track mcs energy bwd
        true,                 // track range momentum muon
        true,                 // track range momentum proton
        true,                 // track range energy muon
        true,                 // track range energy proton
        true,                 // track dep energy yplane
        true,                 // track dep energy yplane minus range energy muassmp
        true,                 // track dep energy yplane minus range energy passmp
        true,                 // track residualrms
        false,                // muon cand length
        false,                // muon cand theta
        false,                // muon cand cos(theta)
        false,                // muon cand phi
        false,                // muon candidate mcs fwd
        false,                // muon candidate mcs bwd
        false,                // muon candidate mc energy forward
        false,                // muon candidate mc energy backward
        false,                // muon candidate momentum from range
        false,                // muon candidate energy from range
        false,                // muon candidate dep energy yplane
        false,                // muon candidate dep energy yplane minus range energy muassmp
        false,                // muon candidate dep energy yplane minus range energy passmp
        false,                // contained muon cand length
        false,                // contained muon cand theta
        false,                // contained muon cand cos(theta)
        false,                // contained muon cand phi
        false,                // contained muon candidate mcs fwd
        false,                // contained muon candidate mcs bwd
        false,                // contained muon candidate mc energy forward
        false,                // contained muon candidate mc energy backward
        false,                // contained muon candidate momentum from range
        false,                // contained muon candidate energy from range
        false,                // contained muon candidate dep energy yplane
        false,                // contained muon candidate dep energy yplane minus range energy muassmp
        false,                // contained muon candidate dep energy yplane minus range energy passmp
        false,                // uncontained muon cand length
        false,                // uncontained muon cand theta
        false,                // uncontained muon cand cos(theta)
        false,                // uncontained muon cand phi
        false,                // uncontained muon candidate mcs fwd
        false,                // uncontained muon candidate mcs bwd
        false,                // uncontained muon candidate mc energy forward
        false,                // uncontained muon candidate mc energy backward
        false,                // uncontained muon candidate momentum from range
        false,                // uncontained muon candidate energy from range
        false,                // uncontained muon candidate dep energy yplane
        false,                // uncontained muon candidate dep energy yplane minus range energy muassmp
        false,                // uncontained muon candidate dep energy yplane minus range energy passmp
        false,                // proton cand length
        false,                // proton cand theta
        false,                // proton cand cos(theta)
        false,                // proton cand phi
        false,                // proton candidate mcs fwd
        false,                // proton candidate mcs bwd
        false,                // proton candidate range momentum
        false,                // proton candidate range energy
        false,                // proton candidate dep energy yplane
        false,                // proton candidate dep energy yplane minus range energy muassmp
        false,                // proton candidate dep energy yplane minus range energy passmp
        false,                // leading proton cand length
        false,                // leading proton cand theta
        false,                // leading proton cand cos(theta)
        false,                // leading proton cand phi
        false,                // leading proton candidate mcs fwd
        false,                // leading proton candidate mcs bwd
        false,                // leading proton candidate range momentum
        false,                // leading proton candidate range energy
        false,                // leading proton candidate dep energy yplane
        false,                // leading proton candidate dep energy yplane minus range energy muassmp
        false,                // leading proton candidate dep energy yplane minus range energy passmp
        false,                // non-leading proton cand length
        false,                // non-leading proton cand theta
        false,                // non-leading proton cand cos(theta)
        false,                // non-leading proton cand phi
        false,                // non-leading proton candidate mcs fwd
        false,                // non-leading proton candidate mcs bwd
        false,                // non-leading proton candidate range momentum
        false,                // non-leading proton candidate range energy
        false,                // non-leading proton candidate dep energy yplane
        false,                // non-leading proton candidate dep energy yplane minus range energy muassmp
        false,                // non-leading proton candidate dep energy yplane minus range energy passmp
        false,                // reconstructed neutrino energy
        false,                // calibrated reconstructed neutrino energy                  
      };

      /*
       * These three vectors fully define all of the efficiency/purity histograms to produce
       */
      std::vector<std::string> effNames = {
        "totel_eff",
        "true_enu",
        "true_mu_p"
      };

      std::vector<std::string> effLabels = {
        ";Total Efficiency;",
        ";E_{#nu}^{true} (Gev);",
        ";p_{#mu} (Gev);"
      };

      std::vector< std::vector<double> > effBins = {
        {1, 0, 1},
        {30, 0, 3},
        {25, 0, 2.5}
      };

      /**
       * Calculate chi2
       */
      std::pair<double,int> CalcChi2(TH1D* E, TH1D* O);

      /**
       * Fill 2D histograms
       */
      void Fill2DHistMC(hists_2d* h2d, std::vector<std::pair<double, double>> variable);

      /**
       * Fill 2D histograms
       */
      void Fill2DHistOnbeam(hists_2d* h2d, std::vector<std::pair<double, double>> variable);

      /**
       * Fill 2D histograms
       */
      void Fill2DHistOffbeam(hists_2d* h2d, std::vector<std::pair<double, double>> variable);

      /**
       * Fill MC histograms
       */
      void FillHistMC(hists_1d* h1d, std::vector<double> variable, std::bitset<8> eventCat);

      /**
       * Fill Dirt histograms
       */
      void FillHistDirt(hists_1d* h1d, std::vector<double> variable);

      /**
       * Fill On-beam histograms
       */
      void FillHistOnBeam(hists_1d* h1d, std::vector<double> variable);

      /**
       * Fill Off-beam histograms
       */
      void FillHistOffBeam(hists_1d* h1d, std::vector<double> variable);

      /**
       * Fill MC trackhistos
       */
      void FillTrackHistMC(trackhists_1d* h1d, std::vector<double> variable, std::vector<double> pid, std::vector<double> cut);

      /**
       * Fill MC dirt trackhistos
       */
      void FillTrackHistDirtMC(trackhists_1d* h1d, std::vector<double> variable, std::vector<double> pid, std::vector<double> cut, float weight);

      /**
       * Fill On-beam trackhistos
       */
      void FillTrackHistOnBeam(trackhists_1d* h1d, std::vector<double> variable, std::vector<double> cut);

      /**
       * Fill Off-beam trackhistos
       */
      void FillTrackHistOffBeam(trackhists_1d* h1d, std::vector<double> variable, std::vector<double> cut);

      /**
       * Initialise vector of hists_1ds
       */
      void InitialiseHistoVec(std::vector< std::vector<hists_1d*> >* plots_to_make, int n_plots);

      /**
       * Initialise vecotr of trackhists_1ds
       */
      void InitialiseTrackHistoVec(std::vector< std::vector<trackhists_1d*> >* trackhists_plots_to_make, int ntrackplots);

      /**
       * Initialise vector of TH2s
       */
      void InitialiseHistoVec(std::vector< std::vector<hists_2d*> >* plots_to_make, int n_plots);

      /**
       * Initialise vector of eff_1ds
       */
      void InitialiseHistoVec(std::vector< std::vector<eff_1d*> >* plots_to_make, int n_plots);

      /**
       * Style histograms
       */
      void StyleHistograms(hists_1d* hists);

      /**
       * Style track histograms
       */
      void StyleTrackHistograms(trackhists_1d* hists);

      /**
       * Scale histograms
       */
      void ScaleHistograms(std::vector< std::vector<hists_1d*> > hists);

      /**
       * Scale track histograms
       */
      void ScaleTrackHistograms(std::vector< std::vector<trackhists_1d*> > hists);

      /**
       * Make stacked data/mc histograms and save to png/pdf file
       */
      void MakeStackedHistogramsAndSave(std::vector< std::vector<hists_1d*> > hists);

      /**
       * Make stacked data/mc track histograms and save to png/pdf file
       */
      void MakeStackedTrackHistogramsAndSave(std::vector< std::vector<trackhists_1d*> > hists);

      /**
       * Make 2d histograms and save to png/pdf file
       */
      void Make2DHistogramsAndSave(std::vector< std::vector<hists_2d*> > hists);

      /**
       * Make efficiency histograms and save to png/pdf files
       */
      void MakeEfficiencyHistogramsAndSave(std::vector< std::vector<eff_1d*> > effplots);

  };

}

#endif
