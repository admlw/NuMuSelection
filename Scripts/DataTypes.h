#ifndef DATATYPES_H
#define DATATYPES_H

// cpp
#include <vector>
#include <string>

// ROOT
#include "TH1.h"

// define structs to use
struct var_list {

  bool isUBXSecSelected = false;
  int nSelectedTracks = -999;
  int nSelectedShowers = -999;
  std::vector<double>* bragg_fwd_p = nullptr;
  std::vector<double>* bragg_bwd_p = nullptr;
  std::vector<double>* noBragg_fwd_mip = nullptr;
  std::vector<double>* track_length = nullptr;
  double vertex_x = -999;
  double vertex_y = -999;
  double vertex_z = -999;

  bool isBeamNeutrino = false;
  bool isCosmic = false;
  bool isMixed = false;
  bool isInFV = false;
  int true_nu_ccnc = -999;
  std::vector<double>* true_match_pdg = nullptr;
  std::vector<double>* true_mcp_pdg = nullptr;
  std::vector<std::string>* true_mcp_process = nullptr;

};

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

#endif
