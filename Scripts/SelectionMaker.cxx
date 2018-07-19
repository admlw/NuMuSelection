#include "SelectionMaker.h"

namespace numusel{

  std::vector<std::vector<std::vector<double>>> SelectionMaker::GetPlottingVariables(var_list* vars, bool isEffPur, float cutval){

    thisMatrix.clear();
    m_stage0.resize(0);
    m_stage1.resize(0);
    m_stage2.resize(0);
    m_stage3.resize(0);

    numusel::AnalysisCuts anacuts; 
    numusel::SelectionMaker selmaker;

    if (isEffPur == false)
      selmaker.PushBackVVectors(s_stage0, vars, false);
    if (isEffPur == true)
      selmaker.PushBackEPVectors(s_stage0, vars);

    // passes first cut
    if (anacuts.isPassNTrackSelection(vars)){

      if (isEffPur == false)
        selmaker.PushBackVVectors(s_stage1, vars, false);
      if (isEffPur == true)
        selmaker.PushBackEPVectors(s_stage1, vars);

      // passes second cut
      if (anacuts.isPassNShowerSelection(vars)){

        if (isEffPur == false)
          selmaker.PushBackVVectors(s_stage2, vars, false);
        if (isEffPur == true)
          selmaker.PushBackEPVectors(s_stage2, vars);

        // passes third cut
        if (anacuts.isPassParticleIDSelection(vars).first, cutval){

          if (isEffPur == false)
            selmaker.PushBackVVectors(s_stage3, vars, true);
          if (isEffPur == true)
            selmaker.PushBackEPVectors(s_stage3, vars);

        }

      }

    }

    thisMatrix.push_back(m_stage0);
    thisMatrix.push_back(m_stage1);
    thisMatrix.push_back(m_stage2);
    thisMatrix.push_back(m_stage3); 

    return thisMatrix;
  };

  std::vector<std::vector<std::vector<double>>> SelectionMaker::GetPlottingVariables(var_list* vars, bool isEffPur){

      numusel::AnalysisCuts anacuts; 
      thisMatrix = GetPlottingVariables(vars, isEffPur, anacuts.pid_cutvalue);

      return thisMatrix;

  }

  void SelectionMaker::PushBackEPVectors(std::vector<std::vector<double>>* m_stagex, var_list* vars){

    // for genie mcparicles, the first entry is always the neutrino, and the second is the outgoing lepton

    m_stagex->push_back(std::vector<double>({(double)vars->true_genie_starte->at(0)}));
    if (vars->true_genie_startp->size() >= 5)
      m_stagex->push_back(std::vector<double>({(double)vars->true_genie_startp->at(4)}));
    else m_stagex->push_back(std::vector<double>({-1.0}));

  }

  void SelectionMaker::PushBackVVectors(std::vector<std::vector<double>>* m_stagex, var_list* vars, bool isHasPID){

    numusel::AnalysisCuts anacuts;

    m_stagex->push_back(std::vector<double>({(double)vars->nSelectedTracks}));
    m_stagex->push_back(std::vector<double>({(double)vars->nSelectedShowers}));
    m_stagex->push_back(*vars->track_length);
    m_stagex->push_back(std::vector<double>({(double)vars->vertex_x}));
    m_stagex->push_back(std::vector<double>({(double)vars->vertex_y}));
    m_stagex->push_back(std::vector<double>({(double)vars->vertex_z}));
    m_stagex->push_back(*vars->track_startx);
    m_stagex->push_back(*vars->track_endx);
    m_stagex->push_back(*vars->track_starty);
    m_stagex->push_back(*vars->track_endy);
    m_stagex->push_back(*vars->track_startz);
    m_stagex->push_back(*vars->track_endz);
    m_stagex->push_back(*vars->track_theta);
    m_stagex->push_back(*vars->track_costheta);
    m_stagex->push_back(*vars->track_phi);

    std::vector<double> pid;
    for (int i = 0; i < vars->noBragg_fwd_mip->size(); i++){
      double lmip = vars->noBragg_fwd_mip->at(i);
      double lp = std::max(vars->bragg_fwd_p->at(i), vars->bragg_bwd_p->at(i));
      pid.push_back(std::log(lmip/lp));
    }
    m_stagex->push_back(pid);

    m_stagex->push_back(*vars->track_mcs_fwd);
    m_stagex->push_back(*vars->track_mcs_bwd);

    if (isHasPID == true){

      std::vector<double> candMuonLength;
      std::vector<double> candMuonTheta;
      std::vector<double> candMuonCosTheta;
      std::vector<double> candMuonPhi;
      std::vector<double> candMuonMCSFwd;
      std::vector<double> candMuonMCSBwd;

      for (int i = 0; i < pid.size(); i++){

        // get muon candidate
        if (pid.at(i) > anacuts.pid_cutvalue){

          candMuonLength.push_back(vars->track_length->at(i));
          candMuonTheta.push_back(vars->track_theta->at(i));
          candMuonCosTheta.push_back(vars->track_costheta->at(i));
          candMuonPhi.push_back(vars->track_phi->at(i));
          candMuonMCSFwd.push_back(vars->track_mcs_fwd->at(i));
          candMuonMCSBwd.push_back(vars->track_mcs_bwd->at(i));

        }

      }

      m_stagex->push_back(candMuonLength);
      m_stagex->push_back(candMuonTheta);
      m_stagex->push_back(candMuonCosTheta);
      m_stagex->push_back(candMuonPhi);
      m_stagex->push_back(candMuonMCSFwd);
      m_stagex->push_back(candMuonMCSBwd);

    }

  }

}
