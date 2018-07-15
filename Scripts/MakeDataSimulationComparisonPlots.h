// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"

// local
#include "EventCategoriser.h"
#include "DataTypes.h"

// cpp
#include <iostream>

// vars
std::string s_onbeam = "/uboone/data/users/alister1/numuSelection/18-07-09-NuMuSelection/selinfo_onbeam.root";
std::string s_offbeam = "/uboone/data/users/alister1/numuSelection/18-07-09-NuMuSelection/selinfo_offbeam.root";
std::string s_simulation = "/uboone/data/users/alister1/numuSelection/18-07-09-NuMuSelection/selinfo_bnbcos.root";
double offbeamscaling = 0.7; // not good
double simscaling = 1.8; // not good

void setTreeVars(TTree* tree, var_list* varstoset, bool isSimulation){

  tree->SetBranchAddress("nSelectedTracks"  , &(varstoset->nSelectedTracks));
  tree->SetBranchAddress("nSelectedShowers" , &(varstoset->nSelectedShowers));
  tree->SetBranchAddress("track_length"     , &(varstoset->track_length));
  tree->SetBranchAddress("vertex_x"         , &(varstoset->vertex_x));
  tree->SetBranchAddress("vertex_y"         , &(varstoset->vertex_x));
  tree->SetBranchAddress("vertex_z"         , &(varstoset->vertex_x));
  
  if (isSimulation){
    tree->SetBranchAddress("isBeamNeutrino"   , &(varstoset->isBeamNeutrino));
    tree->SetBranchAddress("isCosmic"         , &(varstoset->isCosmic));
    tree->SetBranchAddress("isMixed"          , &(varstoset->isMixed));
    tree->SetBranchAddress("isInFV"           , &(varstoset->isInFV));
    tree->SetBranchAddress("true_match_pdg"   , &(varstoset->true_match_pdg));
    tree->SetBranchAddress("true_mcp_pdg"     , &(varstoset->true_mcp_pdg));
    tree->SetBranchAddress("true_mcp_process" , &(varstoset->true_mcp_process));
  }

}

