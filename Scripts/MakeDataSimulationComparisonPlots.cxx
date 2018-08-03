#include "MakeDataSimulationComparisonPlots.h"

int main(){

  // initialise classes
  numusel::Configuration    _config;
  numusel::EventCategoriser _evcat;
  numusel::SelectionMaker   _selmaker;
  numusel::EfficiencyPurity _effpur;
  numusel::HistogramHandler _histoHandler;
  numusel::TreeHandler      _treeHandler;

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

  _treeHandler.SetTreeVars(t_onbeam, &onbeam_vars_tmp, false);
  _treeHandler.SetTreeVars(t_offbeam, &offbeam_vars_tmp, false);
  _treeHandler.SetTreeVars(t_simulation, &simulation_vars_tmp, true);

  // find the number of plots we're going to make
  int n_plots = (int)_histoHandler.histoNames.size();

  plots_to_make = std::vector<std::vector<hists_1d*> >(_config.n_stages, std::vector<hists_1d*>(n_plots));

  for (int i_st = 0; i_st < _config.n_stages; i_st++){
    for (int i_pl = 0; i_pl < n_plots; i_pl++){

      plots_to_make.at(i_st).at(i_pl) = new hists_1d(
          std::string("h_"+_histoHandler.histoNames.at(i_pl)+"_stage"+std::to_string(i_st)),
          _histoHandler.histoLabels.at(i_pl),
          _histoHandler.histoBins.at(i_pl).at(0),
          _histoHandler.histoBins.at(i_pl).at(1),
          _histoHandler.histoBins.at(i_pl).at(2)
          );
    }
  }

  // find the number of eff/pur plots we're going to make
  int n_effpur = (int)_histoHandler.effpurNames.size();
  eff_to_make = std::vector<std::vector<eff_1d*> >(_config.n_stages, std::vector<eff_1d*>(n_effpur));

  for (int i_st = 0; i_st < _config.n_stages; i_st++){
    for (int i_pl = 0; i_pl < n_effpur; i_pl++){

      eff_to_make.at(i_st).at(i_pl) = new eff_1d(
          std::string("h_"+_histoHandler.effpurNames.at(i_pl)+"_stage"+std::to_string(i_st)),
          _histoHandler.effpurLabels.at(i_pl),
          _histoHandler.effpurBins.at(i_pl).at(0),
          _histoHandler.effpurBins.at(i_pl).at(1),
          _histoHandler.effpurBins.at(i_pl).at(2)
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

    std::vector<std::vector<std::vector<double>>> effVariables = _selmaker.GetPlottingVariables(simulation_vars, true); 


    for (size_t i_st = 0; i_st < effVariables.size(); i_st++){

      for (size_t i_pl = 0; i_pl < effVariables.at(i_st).size(); i_pl++){

        _effpur.FillEfficiencyNumerator(eff_to_make.at(i_st).at(i_pl), effVariables.at(i_st).at(i_pl), eventCat, simulation_vars);

        _effpur.FillEfficiencyDenominator(eff_to_make.at(i_st).at(i_pl), effVariables.at(i_st).at(i_pl), eventCat, simulation_vars);


      }
    }

    // get plots for data/MC comparisons
    if (simulation_vars->isUBXSecSelected){
      std::vector<std::vector<std::vector<double>>> plottingVariables = _selmaker.GetPlottingVariables(simulation_vars, false);

      for (size_t i_st = 0; i_st < plottingVariables.size(); i_st++){ 
        for (size_t i_pl = 0; i_pl < plottingVariables.at(i_st).size(); i_pl++){

          _histoHandler.FillHistMC(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), eventCat);

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

          _histoHandler.FillHistOnBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));
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

          _histoHandler.FillHistOffBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));

        }
      }
    }
  }

  _histoHandler.ScaleHistograms(plots_to_make);
  _histoHandler.MakeStackedHistogramsAndSave(plots_to_make);
  _histoHandler.MakeEfficiencyHistogramsAndSave(eff_to_make);

  std::cout << "efficiency: " << eff_to_make.at(3).at(0)->h_num->GetBinContent(1)/eff_to_make.at(0).at(0)->h_denom->GetBinContent(1) << std::endl;;

  TH1D *h_tot = (TH1D*)plots_to_make.at(3).at(0)->h_offbeam->Clone("h_tot");
  h_tot->Add(plots_to_make.at(3).at(0)->h_mccosmic);
  h_tot->Add(plots_to_make.at(3).at(0)->h_mcmixed);
  h_tot->Add(plots_to_make.at(3).at(0)->h_mcoofv);
  h_tot->Add(plots_to_make.at(3).at(0)->h_mcnc);
  h_tot->Add(plots_to_make.at(3).at(0)->h_mcnuenuebar);
  h_tot->Add(plots_to_make.at(3).at(0)->h_mcnumubar);
  h_tot->Add(plots_to_make.at(3).at(0)->h_mcnumuccother);
  h_tot->Add(plots_to_make.at(3).at(0)->h_mcnumucc0pinp);

  std::cout << "purity:" << plots_to_make.at(3).at(0)->h_mcnumucc0pinp->GetBinContent(1)/h_tot->GetBinContent(1) << std::endl;

  return 0;

}
