#include "SelectionMaker.h"

namespace numusel{

  std::vector<std::vector<double>> SelectionMaker::GetPlottingVariables(var_list* vars){

    thisMatrix.clear();
    m_stage0.clear();
    m_stage1.clear();
    m_stage2.clear();
    m_stage3.clear();

    numusel::AnalysisCuts anacuts; 

    m_stage0.push_back(vars->nSelectedTracks);
    m_stage0.push_back(vars->nSelectedShowers);

    // passes first cut
    if (anacuts.isPassNTrackSelection(vars)){

      m_stage1.push_back(vars->nSelectedTracks);
      m_stage1.push_back(vars->nSelectedShowers);

      // passes second cut
      if (anacuts.isPassNShowerSelection(vars)){

        m_stage2.push_back(vars->nSelectedTracks);
        m_stage2.push_back(vars->nSelectedShowers);

        if (anacuts.isPassParticleIDSelection(vars).first){

          m_stage3.push_back(vars->nSelectedTracks);
          m_stage3.push_back(vars->nSelectedShowers);

        }

      }

    }

    thisMatrix.push_back(m_stage0);
    thisMatrix.push_back(m_stage1);
    thisMatrix.push_back(m_stage2);
    thisMatrix.push_back(m_stage3); 

    return thisMatrix;
  };

}
