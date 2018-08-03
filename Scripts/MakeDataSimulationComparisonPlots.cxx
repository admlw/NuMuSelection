#include "MakeDataSimulationComparisonPlots.h"

int main(){

  // initialise classes
  numusel::Configuration    _config;
  numusel::EventCategoriser _evcat;
  numusel::SelectionMaker   _selmaker;
  numusel::EfficiencyPurity _effpur;
  numusel::HistogramHandler _histhandler;

  // pull out TTrees from provided TFiles
  TTree* t_onbeam = (TTree*)(new TFile(_config.s_onbeam.c_str(), "read"))->Get("numuselection/analysis_tree");
  TTree* t_offbeam = (TTree*)(new TFile(_config.s_offbeam.c_str(), "read"))->Get("numuselection/analysis_tree");
  TTree* t_simulation = (TTree*)(new TFile(_config.s_simulation.c_str(), "read"))->Get("numuselection/analysis_tree");

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

  // find the number of plots we're going to make
  int n_plots = (int)_histhandler.histoNames.size();
  plots_to_make = std::vector<std::vector<hists_1d*> >(_config.n_stages, std::vector<hists_1d*>(n_plots));

  for (int i_st = 0; i_st < _config.n_stages; i_st++){
    for (int i_pl = 0; i_pl < n_plots; i_pl++){

      plots_to_make.at(i_st).at(i_pl) = new hists_1d(
          std::string("h_"+_histhandler.histoNames.at(i_pl)+"_stage"+std::to_string(i_st)),
          _histhandler.histoLabels.at(i_pl),
          _histhandler.histoBins.at(i_pl).at(0),
          _histhandler.histoBins.at(i_pl).at(1),
          _histhandler.histoBins.at(i_pl).at(2)
          );
    }
  }

  // find the number of eff/pur plots we're going to make
  int n_effpur = (int)_histhandler.effpurNames.size();
  eff_to_make = std::vector<std::vector<eff_1d*> >(_config.n_stages, std::vector<eff_1d*>(n_effpur));

  for (int i_st = 0; i_st < _config.n_stages; i_st++){
    for (int i_pl = 0; i_pl < n_effpur; i_pl++){

      eff_to_make.at(i_st).at(i_pl) = new eff_1d(
          std::string("h_"+_histhandler.effpurNames.at(i_pl)+"_stage"+std::to_string(i_st)),
          _histhandler.effpurLabels.at(i_pl),
          _histhandler.effpurBins.at(i_pl).at(0),
          _histhandler.effpurBins.at(i_pl).at(1),
          _histhandler.effpurBins.at(i_pl).at(2)
          );
    }
  }


  // loop simulation
  std::cout << "[MDSCP] Beginning simulation loop..." << std::endl;
  for (int i = 0; i < t_simulation->GetEntries(); i++){

    t_simulation->GetEntry(i);

    // get bitset
    std::bitset<8> eventCat = _evcat.CategoriseEvent(simulation_vars);


    // protect against tracks which don't have PID in the collection plane
    if ( (simulation_vars->nSelectedTracks != simulation_vars->bragg_fwd_p->size()) && simulation_vars->isUBXSecSelected) continue;

    // get eff/pur plots
    // n.b. this is some duplication of code whcih could be improved but just
    // do things the stupid way for now
    std::vector<std::vector<std::vector<double>>> effVariables = _selmaker.GetPlottingVariables(simulation_vars, true); 

    for (size_t i_st = 0; i_st < effVariables.size(); i_st++){
      for (size_t i_pl = 0; i_pl < effVariables.at(i_st).size(); i_pl++){

        _effpur.FillEfficiencyNumerator(eff_to_make.at(i_st).at(i_pl), effVariables.at(i_st).at(i_pl), eventCat, simulation_vars);

        _effpur.FillEfficiencyDenominator(eff_to_make.at(i_st).at(i_pl), effVariables.at(0).at(i_pl), eventCat, simulation_vars);


      }
    }

    // get plots for data/MC comparisons
    if (simulation_vars->isUBXSecSelected){
      std::vector<std::vector<std::vector<double>>> plottingVariables = _selmaker.GetPlottingVariables(simulation_vars, false);

      for (size_t i_st = 0; i_st < plottingVariables.size(); i_st++){ 
        for (size_t i_pl = 0; i_pl < plottingVariables.at(i_st).size(); i_pl++){

          _histhandler.FillHistMC(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), eventCat);

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

      std::vector< std::vector<std::vector<double>> > plottingVariables = _selmaker.GetPlottingVariables(onbeam_vars, false);

      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          _histhandler.FillHistOnBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));
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

      std::vector< std::vector<std::vector<double>> > plottingVariables = _selmaker.GetPlottingVariables(offbeam_vars, false);
      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          _histhandler.FillHistOffBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));

        }
      }
    }
  }

  _histhandler.MakeStackedHistogramsAndSave(plots_to_make);
  _histhandler.MakeEfficiencyHistogramsAndSave(eff_to_make);

  return 0;

}
