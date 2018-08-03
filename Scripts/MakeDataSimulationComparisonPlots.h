// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPad.h"

// local
#include "Configuration.h"
#include "EventCategoriser.h"
#include "DataTypes.h"
#include "AnalysisCuts.h"
#include "SelectionMaker.h"
#include "EfficiencyPurity.h"
#include "HistogramHandler.h"

// cpp
#include <iostream>
#include <bitset>
#include <string>

std::vector<std::vector<hists_1d*>> plots_to_make;
std::vector<std::vector<eff_1d*>> eff_to_make;
std::vector<std::vector<pur_1d*>> pur_to_make;

void setTreeVars(TTree* tree, var_list* varstoset, bool b_isSimulation){

  tree->SetBranchAddress("isUBXSecSelected" , &(varstoset->isUBXSecSelected));
  tree->SetBranchAddress("isSimulation"     , &(varstoset->isSimulation));
  tree->SetBranchAddress("nSelectedTracks"  , &(varstoset->nSelectedTracks));
  tree->SetBranchAddress("nSelectedShowers" , &(varstoset->nSelectedShowers));
  tree->SetBranchAddress("track_length"     , &(varstoset->track_length));
  tree->SetBranchAddress("vertex_x"         , &(varstoset->vertex_x));
  tree->SetBranchAddress("vertex_y"         , &(varstoset->vertex_y));
  tree->SetBranchAddress("vertex_z"         , &(varstoset->vertex_z));
  tree->SetBranchAddress("track_startx"    , &(varstoset->track_startx));
  tree->SetBranchAddress("track_endx"      , &(varstoset->track_endx));
  tree->SetBranchAddress("track_starty"    , &(varstoset->track_starty));
  tree->SetBranchAddress("track_endy"      , &(varstoset->track_endy));
  tree->SetBranchAddress("track_startz"    , &(varstoset->track_startz));
  tree->SetBranchAddress("track_endz"      , &(varstoset->track_endz));
  tree->SetBranchAddress("track_theta"      , &(varstoset->track_theta));
  tree->SetBranchAddress("track_costheta"   , &(varstoset->track_costheta));
  tree->SetBranchAddress("track_phi"        , &(varstoset->track_phi));
  tree->SetBranchAddress("bragg_fwd_p"      , &(varstoset->bragg_fwd_p));
  tree->SetBranchAddress("bragg_bwd_p"      , &(varstoset->bragg_bwd_p));
  tree->SetBranchAddress("noBragg_fwd_mip"  , &(varstoset->noBragg_fwd_mip));
  tree->SetBranchAddress("track_mcs_fwd"    , &(varstoset->track_mcs_fwd));
  tree->SetBranchAddress("track_mcs_bwd"    , &(varstoset->track_mcs_bwd));

  if (b_isSimulation){
    tree->SetBranchAddress("isBeamNeutrino"   , &(varstoset->isBeamNeutrino));
    tree->SetBranchAddress("isCosmic"         , &(varstoset->isCosmic));
    tree->SetBranchAddress("isMixed"          , &(varstoset->isMixed));
    tree->SetBranchAddress("isInFV"           , &(varstoset->isInFV));
    tree->SetBranchAddress("true_genie_starte", &(varstoset->true_genie_starte));
    tree->SetBranchAddress("true_genie_startp", &(varstoset->true_genie_startp));
    tree->SetBranchAddress("true_genie_pdg"   , &(varstoset->true_genie_pdg));
    tree->SetBranchAddress("true_nu_ccnc"     , &(varstoset->true_nu_ccnc));
    tree->SetBranchAddress("true_match_pdg"   , &(varstoset->true_match_pdg));
    tree->SetBranchAddress("true_mcp_pdg"     , &(varstoset->true_mcp_pdg));
    tree->SetBranchAddress("true_mcp_process" , &(varstoset->true_mcp_process));
    tree->SetBranchAddress("true_mcp_starte"  , &(varstoset->true_mcp_starte));
    tree->SetBranchAddress("true_mcp_startp"  , &(varstoset->true_mcp_startp));
  }

};

