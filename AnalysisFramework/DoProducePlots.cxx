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

  // print configuration information
  std::cout << "--- POT INFORMATION" << std::endl;
  std::cout << "----- simulation POT:          " << _config.bnbcosPOT << std::endl;
  std::cout << "----- onbeam tor860_wcut:      " << _config.onbeam_tor860_wcut << std::endl;
  std::cout << "----- onbeam E1DCNT_wcut:      " << _config.onbeam_E1DCNT_wcut << std::endl;
  std::cout << "----- offbeam EXT:             " << _config.offbeam_EXT << std::endl;
  std::cout << "-----" << std::endl;
  std::cout << "----- simulation scale factor: " << _config.simscaling << std::endl;
  std::cout << "----- offbeam scale factor:    " << _config.offbeamscaling << std::endl;
  std::cout << " " << std::endl;
  std::cout << "--- CONFIGURATION" << std::endl;
  std::cout << "----- Do PID For Tracks?       " << _config.DoPIDForTracks << std::endl;
  std::cout << "----- Do PID For Showers?      " << _config.DoPIDForShowers << std::endl;
  std::cout << "----- Make Track Plots?        " << _config.MakeTrackPlots << std::endl;
  std::cout << "----- Use Track Cut?           " << _config.UseTrackCut << std::endl;
  std::cout << " " << std::endl;
  std::cout << "--- ANALYSIS CUT VALUES" << std::endl;
  std::cout << "----- nPFParticles LOW:        " << _anacuts.n_pfParticles_min_cut_val << std::endl; 
  std::cout << "----- nPFParticles HIGH:       " << _anacuts.n_pfParticles_max_cut_val << std::endl; 
  std::cout << "----- nTracks LOW:             " << _anacuts.n_tracks_min_cut_val << std::endl; 
  std::cout << "----- nTracks HIGH:            " << _anacuts.n_tracks_max_cut_val << std::endl; 
  std::cout << "----- nShowers LOW:            " << _anacuts.n_showers_min_cut_val << std::endl; 
  std::cout << "----- nShowers HIGH:           " << _anacuts.n_showers_max_cut_val << std::endl; 
  std::cout << "----- PID cut value:           " << _anacuts.pid_cutvalue << std::endl; 

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

  TFile *outfile_sim = new TFile("selectedEvents_sim.root", "UPDATE");
  outfile_sim->cd();
  TTree *t_simulation_out = (TTree*)t_simulation->CloneTree(0);
 
  TFile *outfile_onbeam = new TFile("selectedEvents_onbeam.root", "UPDATE");
  outfile_onbeam->cd();
  TTree *t_onbeam_out = (TTree*)t_onbeam->CloneTree(0);

  TFile *outfile_offbeam = new TFile("selectedEvents_offbeam.root", "UPDATE");
  outfile_offbeam->cd();
  TTree *t_offbeam_out = (TTree*)t_offbeam->CloneTree(0);

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

  // ntrackplots
  int n_trackplots = 0;
  for (int i = 0; i < _histoHandler.histoMakeTrackPlot.size(); i++){
    if (_histoHandler.histoMakeTrackPlot.at(i) == true)
      i++;
  }
  trackplots_to_make = std::vector< std::vector<trackhists_1d*> >(_config.n_stages, std::vector<trackhists_1d*>(n_plots));
  _histoHandler.InitialiseTrackHistoVec(&trackplots_to_make, n_plots);

  // 2d plots 
  int n_plots2D = (int)_histoHandler.histoNames_2D.size();
  plots_to_make_2D = std::vector<std::vector<hists_2d*> >(_config.n_stages, std::vector<hists_2d*>(n_plots2D)); 
  _histoHandler.InitialiseHistoVec(&plots_to_make_2D, n_plots2D);

  // initialise efficiency plots to make
  int n_eff = (int)_histoHandler.effNames.size();
  eff_to_make = std::vector<std::vector<eff_1d*> >(_config.n_stages, std::vector<eff_1d*>(n_eff));
  _histoHandler.InitialiseHistoVec(&eff_to_make, n_eff);

  //------------------------------------
  // loop simulation
  //------------------------------------

  std::cout << "[MDSCP] Beginning simulation loop..." << std::endl;
  for (int i = 0; i < t_simulation->GetEntries(); i++){

    t_simulation->GetEntry(i);

    // get bitset
    std::bitset<8> eventCat = _evcat.CategoriseEvent(simulation_vars);
    
    std::vector<std::vector<std::vector<double>>> effVariables = _selmaker.GetPlottingVariables(simulation_vars, EFFICIENCY); 

    for (size_t i_st = 0; i_st < effVariables.size(); i_st++){

      for (size_t i_pl = 0; i_pl < effVariables.at(i_st).size(); i_pl++){

        _eff.FillEfficiencyNumerator(eff_to_make.at(i_st).at(i_pl), effVariables.at(i_st).at(i_pl), eventCat, simulation_vars);

        _eff.FillEfficiencyDenominator(eff_to_make.at(i_st).at(i_pl), effVariables.at(i_st).at(i_pl), eventCat, simulation_vars);


      }
    }

    if (simulation_vars->isUBXSecSelected){

      // get plots for data/MC comparisons
      std::vector<std::vector<std::vector<double>>> plottingVariables = _selmaker.GetPlottingVariables(simulation_vars, HISTOGRAM_1D, t_simulation, t_simulation_out, i);
      std::vector<std::vector<std::vector<double>>> plottingVariablesPDG = _selmaker.GetPlottingVariables(simulation_vars, PDG, t_simulation, t_simulation_out, i);
      std::vector<std::vector<std::vector<double>>> plottingVariablesTrackCut = _selmaker.GetPlottingVariables(simulation_vars, TRACKCUTVAR, t_simulation, t_simulation_out, i);

      for (size_t i_st = 0; i_st < plottingVariables.size(); i_st++){ 

        for (size_t i_pl = 0; i_pl < plottingVariables.at(i_st).size(); i_pl++){

          _histoHandler.FillHistMC(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), eventCat);

          if (_histoHandler.histoMakeTrackPlot.at(i_pl) == true) {

            _histoHandler.FillTrackHistMC(trackplots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), plottingVariablesPDG.at(i_st).at(0), plottingVariablesTrackCut.at(i_st).at(0));

          }
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

      std::vector< std::vector<std::vector<double>> > plottingVariables = _selmaker.GetPlottingVariables(onbeam_vars, HISTOGRAM_1D, t_onbeam, t_onbeam_out, i);
      std::vector<std::vector<std::vector<double>>> plottingVariablesTrackCut = _selmaker.GetPlottingVariables(onbeam_vars, TRACKCUTVAR, t_onbeam, t_onbeam_out, i);


      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          _histoHandler.FillHistOnBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));

          if (_histoHandler.histoMakeTrackPlot.at(i_pl) == true){
            _histoHandler.FillTrackHistOnBeam(trackplots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), plottingVariablesTrackCut.at(i_st).at(0));
          }
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

      std::vector< std::vector<std::vector<double>> > plottingVariables = _selmaker.GetPlottingVariables(offbeam_vars, HISTOGRAM_1D, t_offbeam, t_offbeam_out, i);
      std::vector<std::vector<std::vector<double>>> plottingVariablesTrackCut = _selmaker.GetPlottingVariables(offbeam_vars, TRACKCUTVAR, t_offbeam, t_offbeam_out, i);

      for (int i_st = 0; i_st < (int)plottingVariables.size(); i_st++){ 
        for (int i_pl = 0; i_pl < (int)plottingVariables.at(i_st).size(); i_pl++){

          _histoHandler.FillHistOffBeam(plots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl));

          if (_histoHandler.histoMakeTrackPlot.at(i_pl) == true) {
            _histoHandler.FillTrackHistOffBeam(trackplots_to_make.at(i_st).at(i_pl), plottingVariables.at(i_st).at(i_pl), plottingVariablesTrackCut.at(i_st).at(0));
          }
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

  _histoHandler.ScaleTrackHistograms(trackplots_to_make);
  _histoHandler.MakeStackedTrackHistogramsAndSave(trackplots_to_make);

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

  outfile_onbeam->cd();  
  t_onbeam_out->Write();

  outfile_offbeam->cd();  
  t_offbeam_out->Write();
  
  outfile_sim->cd();  
  t_simulation_out->Write();

  return 0;

}
