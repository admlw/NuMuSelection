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

    int n_nohits = 0;
    for (int j = 0; j < vars->track_dedxperhit_smeared->size(); j++){
      if (vars->track_dedxperhit_smeared->at(j).size() == 0) n_nohits++; 
    }

    std::pair<bool, std::vector<int>> thisPair;

    // first perform quality checks
    if (((vars->nSelectedPfparticles != vars->bragg_fwd_p->size()) && vars->isUBXSecSelected)
        || n_nohits !=0){
      //|| (vars->nSelectedPfparticles != vars->bragg_fwd_p->size())){
      /*
         if ((vars->nSelectedPfparticles != vars->bragg_fwd_p->size())){

         std::cout << "----" << std::endl;
         std::cout << "npfps: " << vars->nSelectedPfparticles << std::endl;
         std::cout << "ntracks: " << vars->nSelectedTracks << std::endl;
         std::cout << "nshowers: " << vars->nSelectedShowers << std::endl;
         std::cout << "npids: " << vars->bragg_fwd_p->size() << std::endl;
         for (int j = 0; j < vars->track_dedxperhit_smeared->size(); j++){
         std::cout << "track has " << vars->track_dedxperhit_smeared->at(j).size() << " dedx values" << std::endl; 
         std::cout << "track has length " << vars->track_length->at(j) << std::endl;
         std::cout << "track pid info: " << vars->bragg_fwd_p->at(j) << std::endl;
         }  

         }
         */
      thisPair.first = false;
      thisPair.second = {0};

      return thisPair;

    }

    int n_muon_cand = 0;
    bool protons_contained = true;
    bool protons_quality = true;
    std::vector<int> pidPdgCodes;
    int n_dedx_lt2=0;
    int n_dedx_gt2=0;

    // calculate PID variables of interest
    for (int i = 0; i < vars->bragg_fwd_p->size(); i++){

      if (_config.DoPIDForShowers == false && vars->pfp_pdgCode->at(i) == 11){
        continue;
      }

      if (_config.DoPIDForTracks == false && vars->pfp_pdgCode->at(i) == 13){
        continue;
      }

      double proton_likelihood = std::max(vars->bragg_fwd_p->at(i), vars->bragg_bwd_p->at(i));
      double mip_likelihood = vars->noBragg_fwd_mip->at(i);

      double loglmipoverp;

      if (mip_likelihood == -999)
        loglmipoverp = -999;
      else
        loglmipoverp = std::log(mip_likelihood/proton_likelihood);

      if (loglmipoverp > cutval){
        n_muon_cand++;
        pidPdgCodes.push_back(13); 
      }
      else {
        pidPdgCodes.push_back(2212);
        if (vars->track_isContained->at(i) == 0)
          protons_contained = false;

        for (int j = 0; j < vars->track_dedxperhit_smeared->at(i).size(); j++){

          if (vars->track_dedxperhit_smeared->at(i).at(j) <= 2)
            n_dedx_lt2++;
          else
            n_dedx_gt2++;

        }

        if ((float)n_dedx_lt2 / (n_dedx_lt2+n_dedx_gt2) > 0.25 && vars->track_isCollectionPID->at(i) == true)
          protons_quality = true;

      }

    }


    if (n_muon_cand == 1 && protons_contained == true &&  protons_quality == true){

      thisPair.first = true;
      thisPair.second = pidPdgCodes;

    }
    else{

      thisPair.first = false;
      thisPair.second = pidPdgCodes;

    }

    return thisPair;

    };

    std::pair<bool, std::vector<int>> AnalysisCuts::isPassParticleIDSelection(var_list* vars){

      std::pair<bool, std::vector<int>> thisPair = isPassParticleIDSelection(vars, pid_cutvalue);
      return thisPair;

    };

  }
