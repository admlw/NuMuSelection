#include "AnalysisCuts.h"

namespace numusel{

  bool AnalysisCuts::isPassNTrackSelection(var_list* vars){

    if (vars->nSelectedTracks >= n_tracks_min_cut_val && 
        vars->nSelectedTracks <= n_tracks_max_cut_val)
      return true;
    else return false;

  };

  bool AnalysisCuts::isPassNShowerSelection(var_list* vars){

    if (vars->nSelectedShowers >= n_showers_min_cut_val &&
        vars->nSelectedShowers <= n_showers_max_cut_val)
       return true;
    else return false;

  };

  std::pair<bool, std::vector<int>> AnalysisCuts::isPassParticleIDSelection(var_list* vars){

    std::pair<bool, std::vector<int>> thisPair;

    int n_muon_cand = 0;
    std::vector<int> pidPdgCodes;

    // calculate PID variables of interest
    for (int i = 0; i < vars->nSelectedTracks; i++){
      double proton_likelihood = std::max(vars->bragg_fwd_p->at(i), vars->bragg_bwd_p->at(i));
      double mip_likelihood = vars->noBragg_fwd_mip->at(i);

      double loglmipoverp = std::log(mip_likelihood/proton_likelihood);
  
      if (loglmipoverp > pid_cutvalue){
        n_muon_cand++;
        pidPdgCodes.push_back(13); 
      }
      else pidPdgCodes.push_back(2212);

    }

    if (n_muon_cand == 1){

      thisPair.first = true;
      thisPair.second = pidPdgCodes;

    }
    else{

      thisPair.first = false;
      thisPair.second = pidPdgCodes;

    }

    return thisPair;

  };

}
