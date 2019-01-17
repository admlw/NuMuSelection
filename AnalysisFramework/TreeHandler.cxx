#include "TreeHandler.h"

namespace numusel{
  
  void TreeHandler::PrepareTreeForSearching(TTree* tree){

    tree->SetBranchStatus("*", 0);

    tree->SetBranchStatus("run", 1);
    tree->SetBranchStatus("subrun", 1);
    tree->SetBranchStatus("event", 1);

  }

  void TreeHandler::PrepareTreeForWriting(TTree* tree){

    tree->SetBranchStatus("*", 1);

  }

  void TreeHandler::SetEWTreeVars(TTree* tree, ew_list* varstoset){
/*
    tree->SetBranchStatus("*", 0);

    tree->SetBranchStatus("run", 1);
    tree->SetBranchStatus("subrun", 1);
    tree->SetBranchStatus("event", 1);
    tree->SetBranchStatus("MCFlux_evtno", 1);
    tree->SetBranchStatus("MCFlux_NuPosX", 1);
    tree->SetBranchStatus("MCFlux_NuPosY", 1);
    tree->SetBranchStatus("MCFlux_NuPosZ", 1);
    tree->SetBranchStatus("MCFlux_NuMomX", 1);
    tree->SetBranchStatus("MCFlux_NuMomY", 1);
    tree->SetBranchStatus("MCFlux_NuMomZ", 1);
    tree->SetBranchStatus("MCFlux_NuMomE", 1);
    tree->SetBranchStatus("MCFlux_genx", 1);
    tree->SetBranchStatus("MCFlux_geny", 1);
    tree->SetBranchStatus("MCFlux_genz", 1);
    tree->SetBranchStatus("MCFlux_ntype", 1);
    tree->SetBranchStatus("MCFlux_ptype", 1);
    tree->SetBranchStatus("MCFlux_nimpwt", 1);
    tree->SetBranchStatus("MCFlux_dk2gen", 1);
    tree->SetBranchStatus("MCFlux_nenergyn", 1);
    tree->SetBranchStatus("MCFlux_tpx", 1);
    tree->SetBranchStatus("MCFlux_tpy", 1);
    tree->SetBranchStatus("MCFlux_tpz", 1);
    tree->SetBranchStatus("MCFlux_tptype", 1);
    tree->SetBranchStatus("MCFlux_vx", 1);
    tree->SetBranchStatus("MCFlux_vy", 1);
    tree->SetBranchStatus("MCFlux_vz", 1);
    tree->SetBranchStatus("MCTruth_NParticles", 1);
    tree->SetBranchStatus("MCTruth_particles_TrackId", 1);
    tree->SetBranchStatus("MCTruth_particles_PdgCode", 1);
    tree->SetBranchStatus("MCTruth_particles_Mother", 1);
    tree->SetBranchStatus("MCTruth_particles_StatusCode", 1);
    tree->SetBranchStatus("MCTruth_particles_NumberDaughters", 1);
    tree->SetBranchStatus("MCTruth_particles_Daughters", 1);
    tree->SetBranchStatus("MCTruth_particles_Gvx", 1);
    tree->SetBranchStatus("MCTruth_particles_Gvy", 1);
    tree->SetBranchStatus("MCTruth_particles_Gvz", 1);
    tree->SetBranchStatus("MCTruth_particles_Gvt", 1);
    tree->SetBranchStatus("MCTruth_particles_px0", 1);
    tree->SetBranchStatus("MCTruth_particles_py0", 1);
    tree->SetBranchStatus("MCTruth_particles_pz0", 1);
    tree->SetBranchStatus("MCTruth_particles_e0", 1); 
    tree->SetBranchStatus("MCTruth_particles_Rescatter", 1);
    tree->SetBranchStatus("MCTruth_particles_polx", 1);
    tree->SetBranchStatus("MCTruth_particles_poly", 1);
    tree->SetBranchStatus("MCTruth_particles_polz", 1);
    tree->SetBranchStatus("MCTruth_neutrino_CCNC", 1);
    tree->SetBranchStatus("MCTruth_neutrino_mode", 1);
    tree->SetBranchStatus("MCTruth_neutrino_interactionType", 1);
    tree->SetBranchStatus("MCTruth_neutrino_target", 1);
    tree->SetBranchStatus("MCTruth_neutrino_nucleon", 1);
    tree->SetBranchStatus("MCTruth_neutrino_quark", 1);
    tree->SetBranchStatus("MCTruth_neutrino_W", 1);
    tree->SetBranchStatus("MCTruth_neutrino_X", 1);
    tree->SetBranchStatus("MCTruth_neutrino_Y", 1);
    tree->SetBranchStatus("MCTruth_neutrino_QSqr", 1);
    tree->SetBranchStatus("GTruth_IsSeaQuark", 1);
    tree->SetBranchStatus("GTruth_tgtPDG", 1);
    tree->SetBranchStatus("GTruth_weight", 1);
    tree->SetBranchStatus("GTruth_probability", 1);
    tree->SetBranchStatus("GTruth_Xsec", 1);
    tree->SetBranchStatus("GTruth_DiffXsec", 1);
    tree->SetBranchStatus("GTruth_vertexX", 1);
    tree->SetBranchStatus("GTruth_vertexY", 1);
    tree->SetBranchStatus("GTruth_vertexZ", 1);
    tree->SetBranchStatus("GTruth_vertexT", 1);
    tree->SetBranchStatus("GTruth_Gscatter", 1);
    tree->SetBranchStatus("GTruth_Gint", 1);
    tree->SetBranchStatus("GTruth_ResNum", 1);
    tree->SetBranchStatus("GTruth_NumPiPlus", 1);
    tree->SetBranchStatus("GTruth_NumPi0", 1);
    tree->SetBranchStatus("GTruth_NumPiMinus", 1);
    tree->SetBranchStatus("GTruth_NumProton", 1);
    tree->SetBranchStatus("GTruth_NumNeutron", 1);
    tree->SetBranchStatus("GTruth_IsCharm", 1);
    tree->SetBranchStatus("GTruth_gX", 1);
    tree->SetBranchStatus("GTruth_gY", 1);
    tree->SetBranchStatus("GTruth_gT", 1);
    tree->SetBranchStatus("GTruth_gW", 1);
    tree->SetBranchStatus("GTruth_gQ2", 1);
    tree->SetBranchStatus("GTruth_gq2", 1);
    tree->SetBranchStatus("GTruth_ProbePDG", 1);
    tree->SetBranchStatus("GTruth_ProbeP4x", 1);
    tree->SetBranchStatus("GTruth_ProbeP4y", 1);
    tree->SetBranchStatus("GTruth_ProbeP4z", 1);
    tree->SetBranchStatus("GTruth_ProbeP4E", 1);
    tree->SetBranchStatus("GTruth_HitNucP4x", 1);
    tree->SetBranchStatus("GTruth_HitNucP4y", 1);
    tree->SetBranchStatus("GTruth_HitNucP4z", 1);
    tree->SetBranchStatus("GTruth_HitNucP4E", 1);
    tree->SetBranchStatus("GTruth_FShadSystP4x", 1);
    tree->SetBranchStatus("GTruth_FShadSystP4y", 1);
    tree->SetBranchStatus("GTruth_FShadSystP4z", 1);
    tree->SetBranchStatus("GTruth_FShadSystP4E", 1);
*/
    tree->SetBranchAddress("run", &(varstoset->run));
    tree->SetBranchAddress("subrun", &(varstoset->subrun));
    tree->SetBranchAddress("event", &(varstoset->event));
    tree->SetBranchAddress("MCFlux_evtno" , &(varstoset->MCFlux_evtno));
    tree->SetBranchAddress("MCFlux_NuPosX", &(varstoset->MCFlux_NuPosX));
    tree->SetBranchAddress("MCFlux_NuPosY", &(varstoset->MCFlux_NuPosY));
    tree->SetBranchAddress("MCFlux_NuPosZ", &(varstoset->MCFlux_NuPosZ));
    tree->SetBranchAddress("MCFlux_NuMomX", &(varstoset->MCFlux_NuMomX));
    tree->SetBranchAddress("MCFlux_NuMomY", &(varstoset->MCFlux_NuMomY));
    tree->SetBranchAddress("MCFlux_NuMomZ", &(varstoset->MCFlux_NuMomZ));
    tree->SetBranchAddress("MCFlux_NuMomE", &(varstoset->MCFlux_NuMomE));
    tree->SetBranchAddress("MCFlux_genx", &(varstoset->MCFlux_genx));
    tree->SetBranchAddress("MCFlux_geny", &(varstoset->MCFlux_geny));
    tree->SetBranchAddress("MCFlux_genz", &(varstoset->MCFlux_genz));
    tree->SetBranchAddress("MCFlux_ntype", &(varstoset->MCFlux_ntype));
    tree->SetBranchAddress("MCFlux_ptype", &(varstoset->MCFlux_ptype));
    tree->SetBranchAddress("MCFlux_nimpwt", &(varstoset->MCFlux_nimpwt));
    tree->SetBranchAddress("MCFlux_dk2gen", &(varstoset->MCFlux_dk2gen));
    tree->SetBranchAddress("MCFlux_nenergyn", &(varstoset->MCFlux_nenergyn));
    tree->SetBranchAddress("MCFlux_tpx", &(varstoset->MCFlux_tpx));
    tree->SetBranchAddress("MCFlux_tpy", &(varstoset->MCFlux_tpy));
    tree->SetBranchAddress("MCFlux_tpz", &(varstoset->MCFlux_tpz));
    tree->SetBranchAddress("MCFlux_tptype", &(varstoset->MCFlux_tptype));
    tree->SetBranchAddress("MCFlux_vx", &(varstoset->MCFlux_vx));
    tree->SetBranchAddress("MCFlux_vy", &(varstoset->MCFlux_vy));
    tree->SetBranchAddress("MCFlux_vz", &(varstoset->MCFlux_vz));
    tree->SetBranchAddress("MCTruth_NParticles", &(varstoset->MCTruth_NParticles));
    tree->SetBranchAddress("MCTruth_particles_TrackId", &(varstoset->MCTruth_particles_TrackId));
    tree->SetBranchAddress("MCTruth_particles_PdgCode", &(varstoset->MCTruth_particles_PdgCode));
    tree->SetBranchAddress("MCTruth_particles_Mother", &(varstoset->MCTruth_particles_Mother));
    tree->SetBranchAddress("MCTruth_particles_StatusCode", &(varstoset->MCTruth_particles_StatusCode));
    tree->SetBranchAddress("MCTruth_particles_NumberDaughters", &(varstoset->MCTruth_particles_NumberDaughters));
    tree->SetBranchAddress("MCTruth_particles_Daughters", &(varstoset->MCTruth_particles_Daughters));
    tree->SetBranchAddress("MCTruth_particles_Gvx", &(varstoset->MCTruth_particles_Gvx));
    tree->SetBranchAddress("MCTruth_particles_Gvy", &(varstoset->MCTruth_particles_Gvy));
    tree->SetBranchAddress("MCTruth_particles_Gvz", &(varstoset->MCTruth_particles_Gvz));
    tree->SetBranchAddress("MCTruth_particles_Gvt", &(varstoset->MCTruth_particles_Gvt));
    tree->SetBranchAddress("MCTruth_particles_px0", &(varstoset->MCTruth_particles_px0));
    tree->SetBranchAddress("MCTruth_particles_py0", &(varstoset->MCTruth_particles_py0));
    tree->SetBranchAddress("MCTruth_particles_pz0", &(varstoset->MCTruth_particles_pz0));
    tree->SetBranchAddress("MCTruth_particles_e0", &(varstoset->MCTruth_particles_e0));
    tree->SetBranchAddress("MCTruth_particles_Rescatter", &(varstoset->MCTruth_particles_Rescatter));
    tree->SetBranchAddress("MCTruth_particles_polx", &(varstoset->MCTruth_particles_polx));
    tree->SetBranchAddress("MCTruth_particles_poly", &(varstoset->MCTruth_particles_poly));
    tree->SetBranchAddress("MCTruth_particles_polz", &(varstoset->MCTruth_particles_polz));
    tree->SetBranchAddress("MCTruth_neutrino_CCNC", &(varstoset->MCTruth_neutrino_CCNC));
    tree->SetBranchAddress("MCTruth_neutrino_mode", &(varstoset->MCTruth_neutrino_mode));
    tree->SetBranchAddress("MCTruth_neutrino_interactionType", &(varstoset->MCTruth_neutrino_interactionType));
    tree->SetBranchAddress("MCTruth_neutrino_target", &(varstoset->MCTruth_neutrino_target));
    tree->SetBranchAddress("MCTruth_neutrino_nucleon", &(varstoset->MCTruth_neutrino_nucleon));
    tree->SetBranchAddress("MCTruth_neutrino_quark", &(varstoset->MCTruth_neutrino_quark));
    tree->SetBranchAddress("MCTruth_neutrino_W", &(varstoset->MCTruth_neutrino_W));
    tree->SetBranchAddress("MCTruth_neutrino_X", &(varstoset->MCTruth_neutrino_X));
    tree->SetBranchAddress("MCTruth_neutrino_Y", &(varstoset->MCTruth_neutrino_Y));
    tree->SetBranchAddress("MCTruth_neutrino_QSqr", &(varstoset->MCTruth_neutrino_QSqr));
    tree->SetBranchAddress("GTruth_IsSeaQuark", &(varstoset->GTruth_IsSeaQuark));
    tree->SetBranchAddress("GTruth_tgtPDG", &(varstoset->GTruth_tgtPDG));
    tree->SetBranchAddress("GTruth_weight", &(varstoset->GTruth_weight));
    tree->SetBranchAddress("GTruth_probability", &(varstoset->GTruth_probability));
    tree->SetBranchAddress("GTruth_Xsec", &(varstoset->GTruth_Xsec));
    tree->SetBranchAddress("GTruth_DiffXsec", &(varstoset->GTruth_DiffXsec));
    tree->SetBranchAddress("GTruth_vertexX", &(varstoset->GTruth_vertexX));
    tree->SetBranchAddress("GTruth_vertexY", &(varstoset->GTruth_vertexY));
    tree->SetBranchAddress("GTruth_vertexZ", &(varstoset->GTruth_vertexZ));
    tree->SetBranchAddress("GTruth_vertexT", &(varstoset->GTruth_vertexT));
    tree->SetBranchAddress("GTruth_Gscatter", &(varstoset->GTruth_Gscatter));
    tree->SetBranchAddress("GTruth_Gint", &(varstoset->GTruth_Gint));
    tree->SetBranchAddress("GTruth_ResNum", &(varstoset->GTruth_ResNum));
    tree->SetBranchAddress("GTruth_NumPiPlus", &(varstoset->GTruth_NumPiPlus));
    tree->SetBranchAddress("GTruth_NumPi0", &(varstoset->GTruth_NumPi0));
    tree->SetBranchAddress("GTruth_NumPiMinus", &(varstoset->GTruth_NumPiMinus));
    tree->SetBranchAddress("GTruth_NumProton", &(varstoset->GTruth_NumProton));
    tree->SetBranchAddress("GTruth_NumNeutron", &(varstoset->GTruth_NumNeutron));
    tree->SetBranchAddress("GTruth_IsCharm", &(varstoset->GTruth_IsCharm));
    tree->SetBranchAddress("GTruth_gX", &(varstoset->GTruth_gX));
    tree->SetBranchAddress("GTruth_gY", &(varstoset->GTruth_gY));
    tree->SetBranchAddress("GTruth_gT", &(varstoset->GTruth_gT));
    tree->SetBranchAddress("GTruth_gW", &(varstoset->GTruth_gW));
    tree->SetBranchAddress("GTruth_gQ2", &(varstoset->GTruth_gQ2));
    tree->SetBranchAddress("GTruth_gq2", &(varstoset->GTruth_gq2));
    tree->SetBranchAddress("GTruth_ProbePDG", &(varstoset->GTruth_ProbePDG));
    tree->SetBranchAddress("GTruth_ProbeP4x", &(varstoset->GTruth_ProbeP4x));
    tree->SetBranchAddress("GTruth_ProbeP4y", &(varstoset->GTruth_ProbeP4y));
    tree->SetBranchAddress("GTruth_ProbeP4z", &(varstoset->GTruth_ProbeP4z));
    tree->SetBranchAddress("GTruth_ProbeP4E", &(varstoset->GTruth_ProbeP4E));
    tree->SetBranchAddress("GTruth_HitNucP4x", &(varstoset->GTruth_HitNucP4x));
    tree->SetBranchAddress("GTruth_HitNucP4y", &(varstoset->GTruth_HitNucP4y));
    tree->SetBranchAddress("GTruth_HitNucP4z", &(varstoset->GTruth_HitNucP4z));
    tree->SetBranchAddress("GTruth_HitNucP4E", &(varstoset->GTruth_HitNucP4E));
    tree->SetBranchAddress("GTruth_FShadSystP4x", &(varstoset->GTruth_FShadSystP4x));
    tree->SetBranchAddress("GTruth_FShadSystP4y", &(varstoset->GTruth_FShadSystP4y));
    tree->SetBranchAddress("GTruth_FShadSystP4z", &(varstoset->GTruth_FShadSystP4z));
    tree->SetBranchAddress("GTruth_FShadSystP4E", &(varstoset->GTruth_FShadSystP4E));

  };

  void TreeHandler::SetTreeVars(TTree* tree, var_list* varstoset, bool b_isSimulation){

    tree->SetBranchStatus("*",0);
   
    tree->SetBranchStatus("run", 1);
    tree->SetBranchStatus("subrun", 1);
    tree->SetBranchStatus("event", 1);

    tree->SetBranchStatus("isUBXSecSelected"     , 1);
    tree->SetBranchStatus("isSimulation"         , 1);
    tree->SetBranchStatus("nSelectedTracks"      , 1);
    tree->SetBranchStatus("nSelectedShowers"     , 1);
    tree->SetBranchStatus("nSelectedPfparticles" , 1);
    tree->SetBranchStatus("vertex_x"             , 1);
    tree->SetBranchStatus("vertex_y"             , 1);
    tree->SetBranchStatus("vertex_z"             , 1);
    tree->SetBranchStatus("pfp_pdgCode"          , 1);
    tree->SetBranchStatus("track_length"         , 1);
    tree->SetBranchStatus("track_startx"         , 1);
    tree->SetBranchStatus("track_endx"           , 1);
    tree->SetBranchStatus("track_starty"         , 1);
    tree->SetBranchStatus("track_endy"           , 1);
    tree->SetBranchStatus("track_startz"         , 1);
    tree->SetBranchStatus("track_endz"           , 1);
    tree->SetBranchStatus("track_theta"          , 1);
    tree->SetBranchStatus("track_costheta"       , 1);
    tree->SetBranchStatus("track_isContained"    , 1);
    tree->SetBranchStatus("track_isCollectionPID", 1);
    tree->SetBranchStatus("track_phi"            , 1);
    tree->SetBranchStatus("bragg_fwd_p"          , 1);
    tree->SetBranchStatus("bragg_bwd_p"          , 1);
    tree->SetBranchStatus("noBragg_fwd_mip"      , 1);
    tree->SetBranchStatus("track_mcs_muassmp_fwd_loglikelihood", 1);
    tree->SetBranchStatus("track_mcs_muassmp_bwd_loglikelihood", 1);
    tree->SetBranchStatus("track_mcs_muassmp_fwd"        , 1);
    tree->SetBranchStatus("track_mcs_muassmp_bwd"        , 1);
    tree->SetBranchStatus("track_mcs_muassmp_energy_fwd", 1);
    tree->SetBranchStatus("track_mcs_muassmp_energy_bwd", 1);
    tree->SetBranchStatus("track_range_mom_muassumption", 1);
    tree->SetBranchStatus("track_range_mom_passumption", 1);
    tree->SetBranchStatus("track_range_energy_muassumption", 1);
    tree->SetBranchStatus("track_range_energy_passumption", 1);
    tree->SetBranchStatus("track_dedxperhit_smeared", 1);
    tree->SetBranchStatus("track_resrangeperhit", 1);
    tree->SetBranchStatus("track_residualrms", 1);
    tree->SetBranchStatus("track_hit_nhits_uplane", 1);
    tree->SetBranchStatus("track_hit_nhits_vplane", 1);
    tree->SetBranchStatus("track_hit_nhits_yplane", 1);
    tree->SetBranchStatus("track_hit_median_peak_amplitude_uplane", 1);
    tree->SetBranchStatus("track_hit_median_peak_amplitude_vplane", 1);
    tree->SetBranchStatus("track_hit_median_peak_amplitude_yplane", 1);
    tree->SetBranchStatus("track_hit_median_integral_uplane", 1);
    tree->SetBranchStatus("track_hit_median_integral_vplane", 1);
    tree->SetBranchStatus("track_hit_median_integral_yplane", 1);
    tree->SetBranchStatus("track_hit_median_multiplicity_uplane", 1);
    tree->SetBranchStatus("track_hit_median_multiplicity_vplane", 1);
    tree->SetBranchStatus("track_hit_median_multiplicity_yplane", 1);
    tree->SetBranchStatus("track_ncaloobj_uplane", 1);
    tree->SetBranchStatus("track_ncaloobj_vplane", 1);
    tree->SetBranchStatus("track_ncaloobj_yplane", 1);

    tree->SetBranchAddress("run",    &(varstoset->run));
    tree->SetBranchAddress("subrun", &(varstoset->subrun));
    tree->SetBranchAddress("event",  &(varstoset->event));
    tree->SetBranchAddress("isUBXSecSelected"     , &(varstoset->isUBXSecSelected));
    tree->SetBranchAddress("isSimulation"         , &(varstoset->isSimulation));
    tree->SetBranchAddress("nSelectedTracks"      , &(varstoset->nSelectedTracks));
    tree->SetBranchAddress("nSelectedShowers"     , &(varstoset->nSelectedShowers));
    tree->SetBranchAddress("nSelectedPfparticles" , &(varstoset->nSelectedPfparticles));
    tree->SetBranchAddress("vertex_x"             , &(varstoset->vertex_x));
    tree->SetBranchAddress("vertex_y"             , &(varstoset->vertex_y));
    tree->SetBranchAddress("vertex_z"             , &(varstoset->vertex_z));
    tree->SetBranchAddress("pfp_pdgCode"          , &(varstoset->pfp_pdgCode));
    tree->SetBranchAddress("track_length"         , &(varstoset->track_length));
    tree->SetBranchAddress("track_startx"         , &(varstoset->track_startx));
    tree->SetBranchAddress("track_endx"           , &(varstoset->track_endx));
    tree->SetBranchAddress("track_starty"         , &(varstoset->track_starty));
    tree->SetBranchAddress("track_endy"           , &(varstoset->track_endy));
    tree->SetBranchAddress("track_startz"         , &(varstoset->track_startz));
    tree->SetBranchAddress("track_endz"           , &(varstoset->track_endz));
    tree->SetBranchAddress("track_theta"          , &(varstoset->track_theta));
    tree->SetBranchAddress("track_costheta"       , &(varstoset->track_costheta));
    tree->SetBranchAddress("track_phi"            , &(varstoset->track_phi));
    tree->SetBranchAddress("track_isContained"    , &(varstoset->track_isContained));
    tree->SetBranchAddress("track_isCollectionPID", &(varstoset->track_isCollectionPID));
    tree->SetBranchAddress("bragg_fwd_p"          , &(varstoset->bragg_fwd_p));
    tree->SetBranchAddress("bragg_bwd_p"          , &(varstoset->bragg_bwd_p));
    tree->SetBranchAddress("noBragg_fwd_mip"      , &(varstoset->noBragg_fwd_mip));
    tree->SetBranchAddress("track_mcs_muassmp_fwd"        , &(varstoset->track_mcs_muassmp_fwd));
    tree->SetBranchAddress("track_mcs_muassmp_bwd"        , &(varstoset->track_mcs_muassmp_bwd));
    tree->SetBranchAddress("track_mcs_muassmp_fwd_loglikelihood", &(varstoset->track_mcs_muassmp_fwd_loglikelihood));
    tree->SetBranchAddress("track_mcs_muassmp_bwd_loglikelihood", &(varstoset->track_mcs_muassmp_bwd_loglikelihood));
    tree->SetBranchAddress("track_mcs_muassmp_energy_fwd", &(varstoset->track_mcs_muassmp_energy_fwd));
    tree->SetBranchAddress("track_mcs_muassmp_energy_bwd", &(varstoset->track_mcs_muassmp_energy_bwd));
    tree->SetBranchAddress("track_range_mom_muassumption", &(varstoset->track_range_mom_muassumption));
    tree->SetBranchAddress("track_range_mom_passumption", &(varstoset->track_range_mom_passumption));
    tree->SetBranchAddress("track_range_energy_muassumption", &(varstoset->track_range_energy_muassumption));
    tree->SetBranchAddress("track_range_energy_passumption", &(varstoset->track_range_energy_passumption));
    tree->SetBranchAddress("track_dedxperhit_smeared", &(varstoset->track_dedxperhit_smeared));
    tree->SetBranchAddress("track_resrangeperhit", &(varstoset->track_resrangeperhit));
    tree->SetBranchAddress("track_residualrms", &(varstoset->track_residualrms));
    tree->SetBranchAddress("track_hit_nhits_uplane", &(varstoset->track_hit_nhits_uplane));
    tree->SetBranchAddress("track_hit_nhits_vplane", &(varstoset->track_hit_nhits_vplane));
    tree->SetBranchAddress("track_hit_nhits_yplane", &(varstoset->track_hit_nhits_yplane));
    tree->SetBranchAddress("track_hit_median_peak_amplitude_uplane", &(varstoset->track_hit_median_peak_amplitude_uplane));
    tree->SetBranchAddress("track_hit_median_peak_amplitude_vplane", &(varstoset->track_hit_median_peak_amplitude_vplane));
    tree->SetBranchAddress("track_hit_median_peak_amplitude_yplane", &(varstoset->track_hit_median_peak_amplitude_yplane));
    tree->SetBranchAddress("track_hit_median_integral_uplane"      , &(varstoset->track_hit_median_integral_uplane));
    tree->SetBranchAddress("track_hit_median_integral_vplane"      , &(varstoset->track_hit_median_integral_vplane));
    tree->SetBranchAddress("track_hit_median_integral_yplane"      , &(varstoset->track_hit_median_integral_yplane));
    tree->SetBranchAddress("track_hit_median_multiplicity_uplane"  , &(varstoset->track_hit_median_multiplicity_uplane));
    tree->SetBranchAddress("track_hit_median_multiplicity_vplane"  , &(varstoset->track_hit_median_multiplicity_vplane));
    tree->SetBranchAddress("track_hit_median_multiplicity_yplane"  , &(varstoset->track_hit_median_multiplicity_yplane));
    tree->SetBranchAddress("track_ncaloobj_uplane", &(varstoset->track_ncaloobj_uplane));
    tree->SetBranchAddress("track_ncaloobj_vplane", &(varstoset->track_ncaloobj_vplane));
    tree->SetBranchAddress("track_ncaloobj_yplane", &(varstoset->track_ncaloobj_yplane));

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
      tree->SetBranchStatus("true_match_motherid", 1);
      tree->SetBranchStatus("true_match_trackid", 1);
      tree->SetBranchStatus("true_match_process", 1);
      tree->SetBranchStatus("true_match_purity", 1);
      tree->SetBranchStatus("true_match_completeness", 1);
      tree->SetBranchStatus("true_mcp_pdg"     , 1);
      tree->SetBranchStatus("true_mcp_process" , 1);
      tree->SetBranchStatus("true_mcp_starte"  , 1);
      tree->SetBranchStatus("true_mcp_startp"  , 1);
      tree->SetBranchStatus("true_mcp_trackid", 1);

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
      tree->SetBranchAddress("true_match_trackid", &(varstoset->true_match_trackid));
      tree->SetBranchAddress("true_match_motherid", &(varstoset->true_match_motherid));
      tree->SetBranchAddress("true_match_process", &(varstoset->true_match_process));
      tree->SetBranchAddress("true_match_purity", &(varstoset->true_match_purity));
      tree->SetBranchAddress("true_match_completeness", &(varstoset->true_match_completeness));
      tree->SetBranchAddress("true_mcp_pdg"     , &(varstoset->true_mcp_pdg));
      tree->SetBranchAddress("true_mcp_process" , &(varstoset->true_mcp_process));
      tree->SetBranchAddress("true_mcp_starte"  , &(varstoset->true_mcp_starte));
      tree->SetBranchAddress("true_mcp_startp"  , &(varstoset->true_mcp_startp));
      tree->SetBranchAddress("true_mcp_trackid", &(varstoset->true_mcp_trackid));

      tree->Branch("event_cat", "std::vector<bool>", &(varstoset->eventCat));
    }

    tree->Branch("reconstructed_neutrino_energy", &(varstoset->reconstructedNeutrinoEnergy));
    tree->Branch("track_isprotoncand", &(varstoset->track_isprotoncand));
    tree->Branch("track_ismuoncand", &(varstoset->track_ismuoncand));

  };

  int TreeHandler::FindEntryFromEvent(TTree* ewin, ew_list* ewvars, int run, int subrun, int event, int startentry){

    int entry = -1;
    
    for (size_t i = startentry/2; i < ewin->GetEntries(); i++){

      ewin->GetEntry(i);

      if (ewvars->run == run && ewvars->subrun == subrun && ewvars->event == event){

        entry = i;
        break;
      
      }

    }

    return entry;

   };

}
