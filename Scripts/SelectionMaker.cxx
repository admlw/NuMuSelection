#include "SelectionMaker.h"

namespace numusel{

  std::vector<std::vector<std::vector<double>>> SelectionMaker::GetPlottingVariables(var_list* vars){

    thisMatrix.clear();
    m_stage0.resize(0);
    m_stage1.resize(0);
    m_stage2.resize(0);
    m_stage3.resize(0);

    numusel::AnalysisCuts anacuts; 
    numusel::SelectionMaker selmaker;

    selmaker.PushBackVectors(s_stage0, vars);

    // passes first cut
    if (anacuts.isPassNTrackSelection(vars)){
  
      selmaker.PushBackVectors(s_stage1, vars);

      // passes second cut
      if (anacuts.isPassNShowerSelection(vars)){

        selmaker.PushBackVectors(s_stage2, vars);

        if (anacuts.isPassParticleIDSelection(vars).first){

          selmaker.PushBackVectors(s_stage3, vars);

        }

      }

    }

    thisMatrix.push_back(m_stage0);
    thisMatrix.push_back(m_stage1);
    thisMatrix.push_back(m_stage2);
    thisMatrix.push_back(m_stage3); 

    return thisMatrix;
  };

  void SelectionMaker::PushBackVectors(std::vector<std::vector<double>>* m_stagex, var_list* vars){

    m_stagex->push_back(std::vector<double>({(double)vars->nSelectedTracks}));
    m_stagex->push_back(std::vector<double>({(double)vars->nSelectedShowers}));
    m_stagex->push_back(*vars->track_length);

    std::vector<double> pid;
    for (int i = 0; i < vars->noBragg_fwd_mip->size(); i++){
      double lmip = vars->noBragg_fwd_mip->at(i);
      double lp = std::max(vars->bragg_fwd_p->at(i), vars->bragg_bwd_p->at(i));
      pid.push_back(std::log(lmip/lp));
    }
    m_stagex->push_back(pid);

  }

}
