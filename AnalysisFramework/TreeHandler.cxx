#include "TreeHandler.h"

namespace numusel{

  void TreeHandler::SetTreeVars(TTree* tree, var_list* varstoset, bool b_isSimulation){

    tree->SetBranchStatus("*",0);
    
    tree->SetBranchStatus("isUBXSecSelected"     , 1);
    tree->SetBranchStatus("isSimulation"         , 1);
    tree->SetBranchStatus("nSelectedTracks"      , 1);
    tree->SetBranchStatus("nSelectedShowers"     , 1);
    tree->SetBranchStatus("nSelectedPfparticles" , 1);
    tree->SetBranchStatus("track_length"         , 1);
    tree->SetBranchStatus("vertex_x"             , 1);
    tree->SetBranchStatus("vertex_y"             , 1);
    tree->SetBranchStatus("vertex_z"             , 1);
    tree->SetBranchStatus("track_startx"         , 1);
    tree->SetBranchStatus("track_endx"           , 1);
    tree->SetBranchStatus("track_starty"         , 1);
    tree->SetBranchStatus("track_endy"           , 1);
    tree->SetBranchStatus("track_startz"         , 1);
    tree->SetBranchStatus("track_endz"           , 1);
    tree->SetBranchStatus("track_theta"          , 1);
    tree->SetBranchStatus("track_costheta"       , 1);
    tree->SetBranchStatus("track_phi"            , 1);
    tree->SetBranchStatus("bragg_fwd_p"          , 1);
    tree->SetBranchStatus("bragg_bwd_p"          , 1);
    tree->SetBranchStatus("noBragg_fwd_mip"      , 1);
    tree->SetBranchStatus("track_mcs_fwd"        , 1);
    tree->SetBranchStatus("track_mcs_bwd"        , 1);

    tree->SetBranchAddress("isUBXSecSelected"     , &(varstoset->isUBXSecSelected));
    tree->SetBranchAddress("isSimulation"         , &(varstoset->isSimulation));
    tree->SetBranchAddress("nSelectedTracks"      , &(varstoset->nSelectedTracks));
    tree->SetBranchAddress("nSelectedShowers"     , &(varstoset->nSelectedShowers));
    tree->SetBranchAddress("nSelectedPfparticles" , &(varstoset->nSelectedPfparticles));
    tree->SetBranchAddress("track_length"         , &(varstoset->track_length));
    tree->SetBranchAddress("vertex_x"             , &(varstoset->vertex_x));
    tree->SetBranchAddress("vertex_y"             , &(varstoset->vertex_y));
    tree->SetBranchAddress("vertex_z"             , &(varstoset->vertex_z));
    tree->SetBranchAddress("track_startx"         , &(varstoset->track_startx));
    tree->SetBranchAddress("track_endx"           , &(varstoset->track_endx));
    tree->SetBranchAddress("track_starty"         , &(varstoset->track_starty));
    tree->SetBranchAddress("track_endy"           , &(varstoset->track_endy));
    tree->SetBranchAddress("track_startz"         , &(varstoset->track_startz));
    tree->SetBranchAddress("track_endz"           , &(varstoset->track_endz));
    tree->SetBranchAddress("track_theta"          , &(varstoset->track_theta));
    tree->SetBranchAddress("track_costheta"       , &(varstoset->track_costheta));
    tree->SetBranchAddress("track_phi"            , &(varstoset->track_phi));
    tree->SetBranchAddress("bragg_fwd_p"          , &(varstoset->bragg_fwd_p));
    tree->SetBranchAddress("bragg_bwd_p"          , &(varstoset->bragg_bwd_p));
    tree->SetBranchAddress("noBragg_fwd_mip"      , &(varstoset->noBragg_fwd_mip));
    tree->SetBranchAddress("track_mcs_fwd"        , &(varstoset->track_mcs_fwd));
    tree->SetBranchAddress("track_mcs_bwd"        , &(varstoset->track_mcs_bwd));


    if (b_isSimulation){

      tree->SetBranchStatus("isBeamNeutrino"   , 1);
      tree->SetBranchStatus("isCosmic"         , 1);
      tree->SetBranchStatus("isMixed"          , 1);
      tree->SetBranchStatus("isInFV"           , 1);
      tree->SetBranchStatus("true_genie_starte", 1);
      tree->SetBranchStatus("true_genie_startp", 1);
      tree->SetBranchStatus("true_genie_pdg"   , 1);
      tree->SetBranchStatus("true_nu_ccnc"     , 1);
      tree->SetBranchStatus("true_match_pdg"   , 1);
      tree->SetBranchStatus("true_match_starte", 1);
      tree->SetBranchStatus("true_mcp_pdg"     , 1);
      tree->SetBranchStatus("true_mcp_process" , 1);
      tree->SetBranchStatus("true_mcp_starte"  , 1);
      tree->SetBranchStatus("true_mcp_startp"  , 1);

      tree->SetBranchAddress("isBeamNeutrino"   , &(varstoset->isBeamNeutrino));
      tree->SetBranchAddress("isCosmic"         , &(varstoset->isCosmic));
      tree->SetBranchAddress("isMixed"          , &(varstoset->isMixed));
      tree->SetBranchAddress("isInFV"           , &(varstoset->isInFV));
      tree->SetBranchAddress("true_genie_starte", &(varstoset->true_genie_starte));
      tree->SetBranchAddress("true_genie_startp", &(varstoset->true_genie_startp));
      tree->SetBranchAddress("true_genie_pdg"   , &(varstoset->true_genie_pdg));
      tree->SetBranchAddress("true_nu_ccnc"     , &(varstoset->true_nu_ccnc));
      tree->SetBranchAddress("true_match_pdg"   , &(varstoset->true_match_pdg));
      tree->SetBranchAddress("true_match_starte", &(varstoset->true_match_starte));
      tree->SetBranchAddress("true_mcp_pdg"     , &(varstoset->true_mcp_pdg));
      tree->SetBranchAddress("true_mcp_process" , &(varstoset->true_mcp_process));
      tree->SetBranchAddress("true_mcp_starte"  , &(varstoset->true_mcp_starte));
      tree->SetBranchAddress("true_mcp_startp"  , &(varstoset->true_mcp_startp));

    }

  };

}
