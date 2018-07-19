#include "MakeDataSimulationComparisonPlots.h"

int main(){

  // pull out TTrees from provided TFiles
  TTree* t_onbeam = (TTree*)(new TFile(s_onbeam.c_str(), "read"))->Get("numuselection/analysis_tree");
  TTree* t_offbeam = (TTree*)(new TFile(s_offbeam.c_str(), "read"))->Get("numuselection/analysis_tree");
  TTree* t_simulation = (TTree*)(new TFile(s_simulation.c_str(), "read"))->Get("numuselection/analysis_tree");

  // initialise variables and trees
  var_list onbeam_vars_tmp;
  var_list offbeam_vars_tmp;
  var_list simulation_vars_tmp;

  var_list* onbeam_vars = &onbeam_vars_tmp;
  var_list* offbeam_vars = &offbeam_vars_tmp;
  var_list* simulation_vars = &simulation_vars_tmp;

  setTreeVars(t_onbeam, &onbeam_vars_tmp, false);
  setTreeVars(t_offbeam, &offbeam_vars_tmp, false);
  setTreeVars(t_simulation, &simulation_vars_tmp, true);

  // initialise classes
  numusel::EventCategoriser evcat;
  numusel::SelectionMaker selmaker;

  // find the number of plots we're going to make
  int n_plots = (int)histoNames.size();
  plots_to_make = std::vector<std::vector<hists_1d*> >(n_stages, std::vector<hists_1d*>(n_plots));

  for (int i_st = 0; i_st < n_stages; i_st++){
    for (int i_pl = 0; i_pl < n_plots; i_pl++){
  
      plots_to_make.at(i_st).at(i_pl) = new hists_1d(
          std::string("h_"+histoNames.at(i_pl)+"_stage"+std::to_string(i_st)),
          histoLabels.at(i_pl),
          histoBins.at(i_pl).at(0),
          histoBins.at(i_pl).at(1),
          histoBins.at(i_pl).at(2)
          );
    }
  }

  // loop simulation
  std::cout << "[MDSCP] Beginning simulation loop..." << std::endl;
  for (int i = 0; i < t_simulation->GetEntries(); i++){

    t_simulation->GetEntry(i);

    if (simulation_vars->isUBXSecSelected){

      if (simulation_vars->nSelectedTracks != simulation_vars->bragg_fwd_p->size()) continue;

      // get bitset
      std::bitset<8> eventCat = evcat.CategoriseEvent(simulation_vars);

      std::vector<std::vector<std::vector<double>>> plottingVariables = selmaker.GetPlottingVariables(simulation_vars);

      for (size_t i_st = 0; i_st < plottingVariables.size(); i_st++){ 
        for (size_t i_pl = 0; i_pl < plottingVariables.at(i_st).size(); i_pl++){

          FillHistMC(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), eventCat);

        }
      }
    }
  }

  // loop onbeam
  std::cout << "[MDSCP] Beginning onbeam loop..." << std::endl;
  for (int i = 0; i < t_onbeam->GetEntries(); i++){

    t_onbeam->GetEntry(i);

    if (onbeam_vars->isUBXSecSelected){

      if (onbeam_vars->nSelectedTracks != onbeam_vars->bragg_fwd_p->size()) continue;

      std::vector< std::vector<std::vector<double>> > plottingVariables = selmaker.GetPlottingVariables(onbeam_vars);

      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          FillHistOnBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));
        }
      }
    }
  }

  // loop offbeam
  std::cout << "[MDSCP] Beginning offbeam loop..." << std::endl;
  for (int i = 0; i < t_offbeam->GetEntries(); i++){

    t_offbeam->GetEntry(i);

    if (offbeam_vars->isUBXSecSelected){

      if (offbeam_vars->nSelectedTracks != offbeam_vars->bragg_fwd_p->size()) continue;

      std::vector< std::vector<std::vector<double>> > plottingVariables = selmaker.GetPlottingVariables(offbeam_vars);
      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          FillHistOffBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));

        }
      }
    }
  }

  MakeStackedHistogramsAndSave(plots_to_make);

  return 0;

}