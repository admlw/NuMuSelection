/**
 * \author Adam Lister
 *
 * \email a.lister1@lancaster.ac.uk
 *
 * \description this framework is designed to work with the output of the 
 *              NuMuSelection module contained in the Module directory
 *              in this repository
 *              it pulls out many different data/mc comparisons at
 *              each stage of the selection
 *              the stages are currently set to be
 *                0. pure UBXSec CC-inclusive
 *                1. Topology cut: N tracks
 *                2. Topology cut: N showers
 *                3. ParticleID cut
 */

#include "Main.h"

int main(){

  // pull out TTrees from provided TFiles
  // input trees are set in Configuration.h
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

  // initialise data/mc comparison plots to make
  int n_plots = (int)_histoHandler.histoNames.size();
  plots_to_make = std::vector<std::vector<hists_1d*> >(_config.n_stages, std::vector<hists_1d*>(n_plots));
  _histoHandler.InitialiseHistoVec(&plots_to_make, n_plots);

  // initialise efficiency plots to make
  int n_effpur = (int)_histoHandler.effpurNames.size();
  eff_to_make = std::vector<std::vector<eff_1d*> >(_config.n_stages, std::vector<eff_1d*>(n_effpur));
  _histoHandler.InitialiseHistoVec(&eff_to_make, n_effpur);

  //------------------------------------
  // loop simulation
  //------------------------------------
 
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

  //------------------------------------
  // loop onbeam
  //------------------------------------
 
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

  //------------------------------------
  // loop offbeam
  //------------------------------------

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

  //------------------------------------
  // scale histograms, draw, and save
  //------------------------------------

  _histoHandler.ScaleHistograms(plots_to_make);
  _histoHandler.MakeStackedHistogramsAndSave(plots_to_make);
  _histoHandler.MakeEfficiencyHistogramsAndSave(eff_to_make);

  //------------------------------------
  // print efficiency and purity info
  //------------------------------------

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
