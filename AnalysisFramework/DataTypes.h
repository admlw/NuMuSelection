#ifndef DATATYPES_H
#define DATATYPES_H

// cpp
#include <vector>
#include <bitset>
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

  // ---
  // reconstructed
  // ---

  // --- one number per event
  int run = -999;
  int subrun = -999;
  int event = -999;
  bool isUBXSecSelected = false;
  bool isSimulation = false;
  int nSelectedTracks = -999;
  int nSelectedShowers = -999;
  int nSelectedPfparticles = -999;
  double vertex_x = -999;
  double vertex_y = -999;
  double vertex_z = -999;

  // --- pfp info
  std::vector<int>* pfp_pdgCode = nullptr;

  // --- track info
  std::vector<double>* track_length = nullptr;
  std::vector<double>* track_startx = nullptr;
  std::vector<double>* track_endx = nullptr;
  std::vector<double>* track_starty = nullptr;
  std::vector<double>* track_endy = nullptr;
  std::vector<double>* track_startz = nullptr;
  std::vector<double>* track_endz = nullptr;
  std::vector<double>* track_theta = nullptr;
  std::vector<double>* track_costheta = nullptr;
  std::vector<double>* track_phi = nullptr;
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
  std::vector<bool>* track_isCollectionPID = nullptr;
  std::vector<std::vector<double>>* track_dedxperhit_smeared = nullptr;
  std::vector<std::vector<double>>* track_resrangeperhit = nullptr;
  std::vector<double>* bragg_fwd_p = nullptr;
  std::vector<double>* bragg_bwd_p = nullptr;
  std::vector<double>* noBragg_fwd_mip = nullptr;

  // --- track hit variables
  std::vector<int>* track_hit_nhits_uplane = nullptr;
  std::vector<int>* track_hit_nhits_vplane = nullptr;
  std::vector<int>* track_hit_nhits_yplane = nullptr;
  std::vector<double>* track_hit_median_peak_amplitude_uplane = nullptr;
  std::vector<double>* track_hit_median_peak_amplitude_vplane = nullptr;
  std::vector<double>* track_hit_median_peak_amplitude_yplane = nullptr;
  std::vector<double>* track_hit_median_integral_uplane = nullptr;
  std::vector<double>* track_hit_median_integral_vplane = nullptr;
  std::vector<double>* track_hit_median_integral_yplane = nullptr;
  std::vector<double>* track_hit_median_multiplicity_uplane = nullptr;
  std::vector<double>* track_hit_median_multiplicity_vplane = nullptr;
  std::vector<double>* track_hit_median_multiplicity_yplane = nullptr;

  // --- track calorimetry objects
  std::vector<int>* track_ncaloobj_uplane = nullptr;
  std::vector<int>* track_ncaloobj_vplane = nullptr;
  std::vector<int>* track_ncaloobj_yplane = nullptr;

  // ---
  // trueth
  // ---
  
  bool isBeamNeutrino = false;
  bool isCosmic = false;
  bool isMixed = false;
  bool isInFV = false;

  // --- genie information
  std::vector<double>* true_genie_starte = nullptr;
  std::vector<double>* true_genie_startp = nullptr;
  std::vector<double>* true_genie_pdg = nullptr;
  int true_nu_ccnc = -999;

  // --- g4 information
  std::vector<double>* true_mcp_pdg = nullptr;
  std::vector<int>* true_mcp_trackid = nullptr;
  std::vector<std::string>* true_mcp_process = nullptr;
  std::vector<double>* true_mcp_starte = nullptr;
  std::vector<double>* true_mcp_startp = nullptr;

  // --- truth match to g4 information
  std::vector<int>* true_match_pdg = nullptr;
  std::vector<double>* true_match_starte = nullptr;
  std::vector<int>* true_match_motherid = nullptr;
  std::vector<int>* true_match_trackid = nullptr;
  std::vector<std::string>* true_match_process = nullptr;
  std::vector<double>* true_match_purity = nullptr;
  std::vector<double>* true_match_completeness = nullptr;

  // ---
  // other
  // ---
  double reconstructedNeutrinoEnergy = -999;
  double reconstructedNeutrinoEnergyCalib = -999;
  std::vector<bool> eventCat;
  std::vector<int> * track_isprotoncand =  nullptr;
  std::vector<int> * track_ismuoncand = nullptr;


};

/**
 * holds all variables to be pulled from the event weight root file
 */
struct ew_list {
  static const int MaxArraySize = 50;

  int run = -999;
  int subrun = -999;
  int event = -999;
  int MCFlux_evtno = -999;
  double MCFlux_NuPosX = -999;
  double MCFlux_NuPosY = -999;
  double MCFlux_NuPosZ = -999;
  double MCFlux_NuMomX = -999;
  double MCFlux_NuMomY = -999;
  double MCFlux_NuMomZ = -999;
  double MCFlux_NuMomE = -999;
  double MCFlux_genx = -999;
  double MCFlux_geny = -999;
  double MCFlux_genz = -999;
  int MCFlux_ntype = -999;
  int MCFlux_ptype = -999;
  double MCFlux_nimpwt = -999;
  double MCFlux_dk2gen = -999;
  double MCFlux_nenergyn = -999;
  double MCFlux_tpx = -999;
  double MCFlux_tpy = -999;
  double MCFlux_tpz = -999;
  int MCFlux_tptype = -999;
  double MCFlux_vx = -999;
  double MCFlux_vy = -999;
  double MCFlux_vz = -999;

  int MCTruth_NParticles = -999;
  int MCTruth_particles_TrackId[MaxArraySize];
  int MCTruth_particles_PdgCode[MaxArraySize];
  int MCTruth_particles_Mother[MaxArraySize];
  int MCTruth_particles_StatusCode[MaxArraySize];
  int MCTruth_particles_NumberDaughters[MaxArraySize];
  int MCTruth_particles_Daughters[MaxArraySize][100];
  double MCTruth_particles_Gvx[MaxArraySize];
  double MCTruth_particles_Gvy[MaxArraySize];
  double MCTruth_particles_Gvz[MaxArraySize];
  double MCTruth_particles_Gvt[MaxArraySize];
  double MCTruth_particles_px0[MaxArraySize];
  double MCTruth_particles_py0[MaxArraySize];
  double MCTruth_particles_pz0[MaxArraySize];
  double MCTruth_particles_e0[MaxArraySize];
  int MCTruth_particles_Rescatter[MaxArraySize];
  double MCTruth_particles_polx[MaxArraySize];
  double MCTruth_particles_poly[MaxArraySize];
  double MCTruth_particles_polz[MaxArraySize];
  int MCTruth_neutrino_CCNC = -999;
  int MCTruth_neutrino_mode = -999;
  int MCTruth_neutrino_interactionType = -999;
  int MCTruth_neutrino_target = -999;
  int MCTruth_neutrino_nucleon = -999;
  int MCTruth_neutrino_quark = -999;
  double MCTruth_neutrino_W = -999;
  double MCTruth_neutrino_X = -999;
  double MCTruth_neutrino_Y = -999;
  double MCTruth_neutrino_QSqr = -999;

  bool GTruth_IsSeaQuark;
  int GTruth_tgtPDG = -999;
  double GTruth_weight = -999;
  double GTruth_probability = -999;
  double GTruth_Xsec = -999;
  double GTruth_DiffXsec = -999;
  double GTruth_vertexX = -999;
  double GTruth_vertexY = -999;
  double GTruth_vertexZ = -999;
  double GTruth_vertexT = -999;
  int GTruth_Gscatter = -999;
  int GTruth_Gint = -999;
  int GTruth_ResNum = -999;
  int GTruth_NumPiPlus = -999;
  int GTruth_NumPi0 = -999;
  int GTruth_NumPiMinus = -999;
  int GTruth_NumProton = -999;
  int GTruth_NumNeutron = -999;
  bool GTruth_IsCharm = -999;
  double GTruth_gX = -999;
  double GTruth_gY = -999;
  double GTruth_gT = -999;
  double GTruth_gW = -999;
  double GTruth_gQ2 = -999;
  double GTruth_gq2 = -999;
  int GTruth_ProbePDG = -999;
  double GTruth_ProbeP4x = -999;
  double GTruth_ProbeP4y = -999;
  double GTruth_ProbeP4z = -999;
  double GTruth_ProbeP4E = -999;
  double GTruth_HitNucP4x = -999;
  double GTruth_HitNucP4y = -999;
  double GTruth_HitNucP4z = -999;
  double GTruth_HitNucP4E = -999;
  double GTruth_FShadSystP4x = -999;
  double GTruth_FShadSystP4y = -999;
  double GTruth_FShadSystP4z = -999;
  double GTruth_FShadSystP4E = -999;

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
  TH1D* h_dirt;
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
    h_dirt = new TH1D(std::string(name+"dirt").c_str(), title.c_str(), nbinsx, binlowx, binhighx);
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
