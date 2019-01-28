#include "AnalysisCuts.h"

namespace numusel{

  bool AnalysisCuts::isPassNPfparticleSelection(var_list* vars){

    if (vars->nSelectedPfparticles >= n_pfParticles_min_cut_val &&
        vars->nSelectedPfparticles <= n_pfParticles_max_cut_val)
      return true;
    else return false;

  };

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

  std::pair<bool, std::vector<int>> AnalysisCuts::isPassParticleIDSelection(var_list* vars, float cutval){

    numusel::Configuration _config;

    // do any tracks have fewer than 3 hits on the collection plane?
    /*int n_lowyplanehittracks = 0;
    for (int j = 0; j < vars->track_hit_nhits_yplane->size(); j++){
        if (vars->track_hit_nhits_yplane->at(j) < min_no_calo_hits) 
            if (vars->track_hit_nhits_uplane->at(j) < 2 || vars->track_hit_nhits_vplane->at(j) < 2)
                n_lowyplanehittracks++;
    }
*/

    std::pair<bool, std::vector<int>> thisPair;

    // first perform quality checks
    if (((vars->nSelectedPfparticles != vars->bragg_fwd_p->size()) && vars->isUBXSecSelected)){
      thisPair.first = false;
      thisPair.second = {0};

      return thisPair;

    }

    int n_muon_cand = 0;
    bool muon_quality = true;
    bool protons_contained = true;
    bool protons_quality = true;
    std::vector<int> pidPdgCodes;

    // calculate PID variables of interest
    for (int i = 0; i < vars->bragg_fwd_p->size(); i++){

      // flags for turning on PID for showers/PID for tracks
      if (_config.DoPIDForShowers == false && vars->pfp_pdgCode->at(i) == 11){
        continue;
      }

      if (_config.DoPIDForTracks == false && vars->pfp_pdgCode->at(i) == 13){
        continue;
      }

      double proton_likelihood = std::max(vars->bragg_fwd_p->at(i), vars->bragg_bwd_p->at(i));
      double mip_likelihood = vars->noBragg_fwd_mip->at(i);

      double loglmipoverp = std::log(mip_likelihood/proton_likelihood);

      // muon candidate
      if (loglmipoverp > cutval && loglmipoverp != 0){
         if (vars->track_dep_energy_yplane->at(i)-vars->track_range_energy_muassumption->at(i) < -0.15 
                || vars->track_dep_energy_yplane->at(i)-vars->track_range_energy_muassumption->at(i) > 0.1)
        muon_quality = false;

        n_muon_cand++;
        pidPdgCodes.push_back(13); 
      }
      // proton candidates
      else {
        pidPdgCodes.push_back(2212);

        if (vars->track_isContained->at(i) == 0)
          protons_contained = false;

        if (vars->track_dep_energy_yplane->at(i)-vars->track_range_energy_passumption->at(i) < -0.15 
                || vars->track_dep_energy_yplane->at(i)-vars->track_range_energy_passumption->at(i) > 0.1)
        protons_quality = false;
      }

    }

    if (n_muon_cand == 1 && muon_quality == true && protons_contained == true && protons_quality == true){

      thisPair.first = true;
      thisPair.second = pidPdgCodes;

    }
    else{

      thisPair.first = false;
      thisPair.second = {0};

    }

    return thisPair;

    };

    std::pair<bool, std::vector<int>> AnalysisCuts::isPassParticleIDSelection(var_list* vars){

      std::pair<bool, std::vector<int>> thisPair = isPassParticleIDSelection(vars, pid_cutvalue);
      return thisPair;

    };

  }
