#ifndef DATATYPES_H
#define DATATYPES_H

// cpp
#include <vector>
#include <string>

// ROOT
#include "TH1.h"
#include "TH2D.h"

/**
 * enum for variable type
 */
enum kVarType {
  HISTOGRAM_1D,
  HISTOGRAM_2D,
  EFFICIENCY,
  PDG,
  TRACKCUTVAR,
};

/**
 * holds all variables to be pulled from the root file for each event
 */
struct var_list {

  bool isUBXSecSelected = false;
  bool isSimulation = false;
  int nSelectedTracks = -999;
  int nSelectedShowers = -999;
  int nSelectedPfparticles = -999;
  std::vector<double>* track_length = nullptr;
  double vertex_x = -999;
  double vertex_y = -999;
  double vertex_z = -999;
  std::vector<int>* pfp_pdgCode = nullptr;
  std::vector<double>* track_startx = nullptr;
  std::vector<double>* track_endx = nullptr;
  std::vector<double>* track_starty = nullptr;
  std::vector<double>* track_endy = nullptr;
  std::vector<double>* track_startz = nullptr;
  std::vector<double>* track_endz = nullptr;
  std::vector<double>* track_theta = nullptr;
  std::vector<double>* track_costheta = nullptr;
  std::vector<double>* track_phi = nullptr;
  std::vector<double>* bragg_fwd_p = nullptr;
  std::vector<double>* bragg_bwd_p = nullptr;
  std::vector<double>* noBragg_fwd_mip = nullptr;
  std::vector<double>* track_mcs_muassmp_fwd = nullptr;
  std::vector<double>* track_mcs_muassmp_bwd = nullptr;
  std::vector<double>* track_mcs_muassmp_fwd_loglikelihood = nullptr;
  std::vector<double>* track_mcs_muassmp_bwd_loglikelihood = nullptr;
  std::vector<double>* track_mcs_muassmp_energy_fwd = nullptr;
  std::vector<double>* track_mcs_muassmp_energy_bwd = nullptr;
  std::vector<double>* track_range_mom_muassumption = nullptr;
  std::vector<double>* track_range_mom_passumption = nullptr;
  std::vector<double>* track_range_energy_muassumption = nullptr;
  std::vector<double>* track_range_energy_passumption = nullptr;
  std::vector<double>* track_residualrms = nullptr;
  std::vector<bool>* track_isContained = nullptr;
  std::vector<std::vector<double>>* track_dedxperhit_smeared = nullptr;
  std::vector<std::vector<double>>* track_resrangeperhit = nullptr;

  bool isBeamNeutrino = false;
  bool isCosmic = false;
  bool isMixed = false;
  bool isInFV = false;
  std::vector<double>* true_genie_starte = nullptr;
  std::vector<double>* true_genie_startp = nullptr;
  std::vector<double>* true_genie_pdg = nullptr;
  int true_nu_ccnc = -999;
  std::vector<int>* true_match_pdg = nullptr;
  std::vector<double>* true_match_starte = nullptr;
  std::vector<double>* true_mcp_pdg = nullptr;
  std::vector<std::string>* true_mcp_process = nullptr;
  std::vector<double>* true_mcp_starte = nullptr;
  std::vector<double>* true_mcp_startp = nullptr;

};

struct hists_2d {
  TH2D* h_mc;
  TH2D* h_onbeam;
  TH2D* h_offbeam;

  hists_2d(std::string name, std::string title, double nbinsx, double binlowx, double binhighx, double nbinsy, double binlowy, double binhighy){
    h_mc = new TH2D(std::string(name+"").c_str(), title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy);
    h_onbeam = new TH2D(std::string(name+"onbeam").c_str(), title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy);
    h_offbeam = new TH2D(std::string(name+"offbeam").c_str(), title.c_str(), nbinsx, binlowx, binhighx, nbinsy, binlowy, binhighy);
  }
};

/**
 * holds all histograms for data-simulation comparisons
 */
struct hists_1d {
  TH1D* h_mccosmic;
  TH1D* h_mcmixed;
  TH1D* h_mcoofv;
  TH1D* h_mcnc;
  TH1D* h_mcnuenuebar;
  TH1D* h_mcnumubar;
  TH1D* h_mcnumuccother;
  TH1D* h_mcnumucc0pinp;
  TH1D* h_onbeam;
  TH1D* h_offbeam;

  hists_1d(std::string name, std::string title, double nbinsx, double binlowx, double binhighx){
    h_mccosmic = new TH1D(std::string(name+"cosmic").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mcmixed = new TH1D(std::string(name+"mixed").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mcoofv = new TH1D(std::string(name+"oofv").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mcnc = new TH1D(std::string(name+"nc").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mcnuenuebar = new TH1D(std::string(name+"nuenuebar").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mcnumubar = new TH1D(std::string(name+"numubar").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mcnumuccother = new TH1D(std::string(name+"numuccother").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mcnumucc0pinp = new TH1D(std::string(name+"numucc0pinp").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_onbeam = new TH1D(std::string(name+"onbeam").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_offbeam = new TH1D(std::string(name+"offbeam").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
  }
};

struct trackhists_1d {
  TH1D* h_muon;
  TH1D* h_proton;
  TH1D* h_pion;
  TH1D* h_kaon;
  TH1D* h_electron;
  TH1D* h_other;
  TH1D* h_onbeam;
  TH1D* h_offbeam;

  trackhists_1d(std::string name, std::string title, double nbinsx, double binlowx, double binhighx){
    h_muon = new TH1D(std::string(name+"muon").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_proton = new TH1D(std::string(name+"proton").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_pion = new TH1D(std::string(name+"pion").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_kaon = new TH1D(std::string(name+"kaon").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_electron = new TH1D(std::string(name+"electron").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_other = new TH1D(std::string(name+"other").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_onbeam = new TH1D(std::string(name+"onbeam").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_offbeam = new TH1D(std::string(name+"offbeam").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
  }
};

struct eff_1d{
  TH1D* h_num;
  TH1D* h_denom;

  eff_1d(std::string name, std::string title, double nbinsx, double binlowx, double binhighx){
    h_num = new TH1D(std::string(name+"_effnum").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_denom = new TH1D(std::string(name+"_effdenom").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
  }

};

struct pur_1d{
  TH1D* h_num;
  TH1D* h_denom;

  pur_1d(std::string name, std::string title, double nbinsx, double binlowx, double binhighx){
    h_num = new TH1D(std::string(name+"_purnum").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_denom = new TH1D(std::string(name+"_purdenom").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
  }

};

#endif