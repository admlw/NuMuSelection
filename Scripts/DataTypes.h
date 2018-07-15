#ifndef DATATYPES_H
#define DATATYPES_H

// cpp
#include <vector>
#include <string>

// ROOT
#include "TH1.h"

// define structs to use
struct var_list {
  int nSelectedTracks = -999;
  int nSelectedShowers = -999;
  std::vector<double>* track_length = nullptr;
  double vertex_x = -999;
  double vertex_y = -999;
  double vertex_z = -999;

  bool isBeamNeutrino = false;
  bool isCosmic = false;
  bool isMixed = false;
  bool isInFV = false;
  int CCNC = 0;
  std::vector<double>* true_match_pdg = nullptr;
  std::vector<double>* true_mcp_pdg = nullptr;
  std::vector<std::string>* true_mcp_process = nullptr;

};

struct hists_1d {
  TH1D* h_cosmic;
  TH1D* h_mixed;
  TH1D* h_oofv;
  TH1D* h_nc;
  TH1D* h_nuenuebar;
  TH1D* h_numubar;
  TH1D* h_numuccother;
  TH1D* h_numucc0pinp;

  hists_1d(std::string name, std::string title, double nbinsx, double binlowx, double binhighx){
    h_cosmic = new TH1D(std::string(name+"cosmic").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_mixed = new TH1D(std::string(name+"mixed").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_oofv = new TH1D(std::string(name+"oofv").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_nc = new TH1D(std::string(name+"nc").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_nuenuebar = new TH1D(std::string(name+"nuenuebar").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_numubar = new TH1D(std::string(name+"numubar").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_numuccother = new TH1D(std::string(name+"numuccother").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
    h_numucc0pinp = new TH1D(std::string(name+"numucc0pinp").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
  }
};

#endif
