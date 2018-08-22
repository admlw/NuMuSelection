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
 *                1. Topology cut: N Pfparticles
 *                2. Topology cut: N tracks
 *                3. Topology cut: N showers
 *                4. ParticleID cut
 */

#include "DoProducePlots.h"

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

  // initialise output trees

  TFile *outfile = new TFile("selectedEvents.root", "RECREATE");
  outfile->cd();

  TTree *t_onbeam_out = (TTree*)t_onbeam->CloneTree(0);
  TTree *t_offbeam_out = (TTree*)t_offbeam->CloneTree(0);
  TTree *t_simulation_out = (TTree*)t_simulation->CloneTree(0);

  //_treeHandler.SetTreeVars(t_onbeam_out, &onbeam_vars_tmp, false);
  //_treeHandler.SetTreeVars(t_offbeam_out, &offbeam_vars_tmp, false);
  //_treeHandler.SetTreeVars(t_simulation_out, &simulation_vars_tmp, true);

  t_onbeam_out->SetName("onbeam");
  t_offbeam_out->SetName("offbeam");
  t_simulation_out->SetName("simulation");

  // debugging
  std::cout << "histoNames size: " << _histoHandler.histoNames.size() << std::endl;
  std::cout << "histoLabels size: " << _histoHandler.histoLabels.size() << std::endl;
  std::cout << "histoBins size: " << _histoHandler.histoBins.size() << std::endl;
  std::cout << "histoNames_2D size: " << _histoHandler.histoNames_2D.size() << std::endl;
  std::cout << "histoLabels_2D size: " << _histoHandler.histoLabels_2D.size() << std::endl;
  std::cout << "histoBins_2D size: " << _histoHandler.histoBins_2D.size() << std::endl;

  // initialise data/mc comparison plots to make
  int n_plots = (int)_histoHandler.histoNames.size();
  plots_to_make = std::vector<std::vector<hists_1d*> >(_config.n_stages, std::vector<hists_1d*>(n_plots));
  _histoHandler.InitialiseHistoVec(&plots_to_make, n_plots);

  // 2d plots 
  int n_plots2D = (int)_histoHandler.histoNames_2D.size();
  plots_to_make_2D = std::vector<std::vector<hists_2d*> >(_config.n_stages, std::vector<hists_2d*>(n_plots2D)); 
  _histoHandler.InitialiseHistoVec(&plots_to_make_2D, n_plots2D);

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
      if (simulation_vars->track_dedxperhit_smeared->size() == 0) continue;
      int n_nohits = 0;
      for (int j = 0; j < simulation_vars->track_dedxperhit_smeared->size(); j++){
        if (simulation_vars->track_dedxperhit_smeared->at(j).size() == 0) n_nohits++; 
      }
      if (n_nohits !=0) continue;

    std::vector<std::vector<std::vector<double>>> effVariables = _selmaker.GetPlottingVariables(simulation_vars, EFFICIENCY); 


    for (size_t i_st = 0; i_st < effVariables.size(); i_st++){

      for (size_t i_pl = 0; i_pl < effVariables.at(i_st).size(); i_pl++){

        _effpur.FillEfficiencyNumerator(eff_to_make.at(i_st).at(i_pl), effVariables.at(i_st).at(i_pl), eventCat, simulation_vars);

        _effpur.FillEfficiencyDenominator(eff_to_make.at(i_st).at(i_pl), effVariables.at(i_st).at(i_pl), eventCat, simulation_vars);


      }
    }

    if (simulation_vars->isUBXSecSelected){

      // get plots for data/MC comparisons
      std::vector<std::vector<std::vector<double>>> plottingVariables = _selmaker.GetPlottingVariables(simulation_vars, HISTOGRAM_1D, t_simulation, t_simulation_out, i);

      for (size_t i_st = 0; i_st < plottingVariables.size(); i_st++){ 

        for (size_t i_pl = 0; i_pl < plottingVariables.at(i_st).size(); i_pl++){

          _histoHandler.FillHistMC(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), eventCat);

        }
      }

      // get plots for 2d histograms
      std::vector<std::vector<std::vector<std::pair<double,double>>>> plottingVariables_2d = _selmaker.Get2DPlottingVariables(simulation_vars, HISTOGRAM_2D, t_simulation, t_simulation_out, i);


      for (size_t i_st = 0; i_st < plottingVariables_2d.size(); i_st++){ 
        for (size_t i_pl = 0; i_pl < plottingVariables_2d.at(i_st).size(); i_pl++){

          _histoHandler.Fill2DHistMC(plots_to_make_2D.at(i_st).at(i_pl), plottingVariables_2d.at(i_st).at(i_pl));

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
      if (onbeam_vars->track_dedxperhit_smeared->size() == 0) continue;
      int n_nohits = 0;
      for (int j = 0; j < onbeam_vars->track_dedxperhit_smeared->size(); j++){
        if (onbeam_vars->track_dedxperhit_smeared->at(j).size() == 0) n_nohits++; 
      }
      if (n_nohits !=0) continue;

      std::vector< std::vector<std::vector<double>> > plottingVariables = _selmaker.GetPlottingVariables(onbeam_vars, HISTOGRAM_1D, t_onbeam, t_onbeam_out, i);

      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          _histoHandler.FillHistOnBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));
        }
      }

      // get plots for 2d histograms
      std::vector<std::vector<std::vector<std::pair<double,double>>>> plottingVariables_2d = _selmaker.Get2DPlottingVariables(onbeam_vars, HISTOGRAM_2D, t_onbeam, t_onbeam_out, i);

      for (size_t i_st = 0; i_st < plottingVariables_2d.size(); i_st++){ 

        for (size_t i_pl = 0; i_pl < plottingVariables_2d.at(i_st).size(); i_pl++){

          _histoHandler.Fill2DHistOnbeam(plots_to_make_2D.at(i_st).at(i_pl), plottingVariables_2d.at(i_st).at(i_pl));

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
      if (offbeam_vars->track_dedxperhit_smeared->size() == 0) continue;
      int n_nohits = 0;
      for (int j = 0; j < offbeam_vars->track_dedxperhit_smeared->size(); j++){
        if (offbeam_vars->track_dedxperhit_smeared->at(j).size() == 0) n_nohits++; 
      }
      if (n_nohits !=0) continue;


      std::vector< std::vector<std::vector<double>> > plottingVariables = _selmaker.GetPlottingVariables(offbeam_vars, HISTOGRAM_1D, t_offbeam, t_offbeam_out, i);
      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          _histoHandler.FillHistOffBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));

        }
      }

      // get plots for 2d histograms
      std::vector<std::vector<std::vector<std::pair<double, double>>>> plottingVariables_2d = _selmaker.Get2DPlottingVariables(offbeam_vars, HISTOGRAM_2D, t_offbeam, t_offbeam_out, i);

      for (size_t i_st = 0; i_st < plottingVariables_2d.size(); i_st++){ 
        for (size_t i_pl = 0; i_pl < plottingVariables_2d.at(i_st).size(); i_pl++){

          _histoHandler.Fill2DHistOffbeam(plots_to_make_2D.at(i_st).at(i_pl), plottingVariables_2d.at(i_st).at(i_pl));

        }
      }


    }
  }

  //------------------------------------
  // scale histograms, draw, and save
  //------------------------------------

  _histoHandler.ScaleHistograms(plots_to_make);
  _histoHandler.MakeStackedHistogramsAndSave(plots_to_make);
  _histoHandler.Make2DHistogramsAndSave(plots_to_make_2D);
  _histoHandler.MakeEfficiencyHistogramsAndSave(eff_to_make);

  //------------------------------------
  // print efficiency and purity info
  //------------------------------------

  int finalStage = plots_to_make.size() - 1;

  std::cout << "efficiency: " << eff_to_make.at(finalStage).at(0)->h_num->GetBinContent(1)/eff_to_make.at(0).at(0)->h_denom->GetBinContent(1) << std::endl;;

  TH1D *h_tot = (TH1D*)plots_to_make.at(finalStage).at(0)->h_offbeam->Clone("h_tot");
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mccosmic);
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mcmixed);
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mcoofv);
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mcnc);
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mcnuenuebar);
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mcnumubar);
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mcnumuccother);
  h_tot->Add(plots_to_make.at(finalStage).at(0)->h_mcnumucc0pinp);

  std::cout << "purity:" << plots_to_make.at(finalStage).at(0)->h_mcnumucc0pinp->GetBinContent(1)/h_tot->GetBinContent(1) << std::endl;

 t_onbeam_out->Write();
 t_offbeam_out->Write();
 t_simulation_out->Write();

  return 0;

}
