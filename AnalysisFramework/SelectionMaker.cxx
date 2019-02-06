#include "SelectionMaker.h"

namespace numusel{

  std::vector<std::vector<std::vector<std::pair<double, double>>>> SelectionMaker::Get2DPlottingVariables(var_list* vars, kVarType var_type, float cutval, TTree* infile, TTree* outfile, int entry){

    thisMatrix_2d.clear();
    m_stage0_2d.resize(0);
    m_stage1_2d.resize(0);
    m_stage2_2d.resize(0);
    m_stage3_2d.resize(0);
    m_stage4_2d.resize(0);

    numusel::AnalysisCuts _anacuts; 
    numusel::SelectionMaker _selmaker;

    switch(var_type){
      case HISTOGRAM_2D:
        _selmaker.PushBack2DVectors(s_stage0_2d, vars, false);
        break;
    }
    // passes first cut
    if (_anacuts.isPassNPfparticleSelection(vars)){

      switch(var_type){
        case HISTOGRAM_2D:
          _selmaker.PushBack2DVectors(s_stage1_2d, vars, false);
          break;
      }
      // passes second cut
      if (_anacuts.isPassNTrackSelection(vars)){

        switch(var_type){
          case HISTOGRAM_2D:
            _selmaker.PushBack2DVectors(s_stage2_2d, vars, false);
            break;
        }

        // passes third cut
        if (_anacuts.isPassNShowerSelection(vars)){
          switch(var_type){
            case HISTOGRAM_2D:
              _selmaker.PushBack2DVectors(s_stage3_2d, vars, false);
              break;
          }

          // passes fourth cut
          if (_anacuts.isPassParticleIDSelection(vars, cutval).first){
            switch(var_type){
              case HISTOGRAM_2D:
                _selmaker.PushBack2DVectors(s_stage4_2d, vars, true);
                break;
            }
          }

        }

      }

    }

    thisMatrix_2d.push_back(m_stage0_2d);
    thisMatrix_2d.push_back(m_stage1_2d);
    thisMatrix_2d.push_back(m_stage2_2d);
    thisMatrix_2d.push_back(m_stage3_2d); 
    thisMatrix_2d.push_back(m_stage4_2d); 

    return thisMatrix_2d;
  };

  std::vector<std::vector<std::vector<std::pair<double,double>>>> SelectionMaker::Get2DPlottingVariables(var_list* vars, kVarType var_type, TTree* infile,  TTree* outfile, int entry){

    numusel::AnalysisCuts _anacuts; 
    thisMatrix_2d = Get2DPlottingVariables(vars, var_type, _anacuts.pid_cutvalue, infile, outfile, entry);

    return thisMatrix_2d;

  }

  std::vector<std::vector<std::vector<double>>> SelectionMaker::GetPlottingVariables(var_list* vars, kVarType var_type, float cutval, TTree* infile, TTree* outfile, int entry, ew_list* ewvars, TTree* ewin, TTree* ewout){

    thisMatrix.clear();
    m_stage0.resize(0);
    m_stage1.resize(0);
    m_stage2.resize(0);
    m_stage3.resize(0);
    m_stage4.resize(0);

    numusel::AnalysisCuts _anacuts; 
    numusel::SelectionMaker _selmaker;
    numusel::TreeHandler _treehandler;
    numusel::Configuration _config;

    switch(var_type){
      case HISTOGRAM_1D:
        _selmaker.PushBack1DVectors(s_stage0, vars, false);
        break;
      case EFFICIENCY:
        _selmaker.PushBackEPVectors(s_stage0, vars);
        break;
      case PDG:
        _selmaker.PushBackPIDVectors(s_stage0, vars);
        break;
      case TRACKCUTVAR:
        _selmaker.PushBackTrackCutVar(s_stage0, vars);
        break;
    }
    // passes first cut
    // nPFP cut
    if (_anacuts.isPassNPfparticleSelection(vars)){
      switch(var_type){
        case HISTOGRAM_1D:
          _selmaker.PushBack1DVectors(s_stage1, vars, false);
          break;
        case EFFICIENCY:
          _selmaker.PushBackEPVectors(s_stage1, vars);
          break;
        case PDG:
          _selmaker.PushBackPIDVectors(s_stage1, vars);
          break;
        case TRACKCUTVAR:
          _selmaker.PushBackTrackCutVar(s_stage1, vars);
          break;
      }
      // passes second cut
      // nTrack cut (obsolete)
      if (_anacuts.isPassNTrackSelection(vars)){
        switch(var_type){
          case HISTOGRAM_1D:
            _selmaker.PushBack1DVectors(s_stage2, vars, false);
            break;
          case EFFICIENCY:
            _selmaker.PushBackEPVectors(s_stage2, vars);
            break;
          case PDG:
            _selmaker.PushBackPIDVectors(s_stage2, vars);
            break;
          case TRACKCUTVAR:
            _selmaker.PushBackTrackCutVar(s_stage2, vars);
            break;
        }

        // passes third cut
        // nShower cut (obslete)
        if (_anacuts.isPassNShowerSelection(vars)){
          switch(var_type){
            case HISTOGRAM_1D:
              _selmaker.PushBack1DVectors(s_stage3, vars, false);
              break;
            case EFFICIENCY:
              _selmaker.PushBackEPVectors(s_stage3, vars);
              break;
            case PDG:
              _selmaker.PushBackPIDVectors(s_stage3, vars);
              break;
            case TRACKCUTVAR:
              _selmaker.PushBackTrackCutVar(s_stage3, vars);
              break;

          }

          // passes fourth cut
          // PID cut and proton quality cuts
          if (_anacuts.isPassParticleIDSelection(vars, cutval).first){
            bool hasPid = false;
            switch(var_type){
              case HISTOGRAM_1D:
                hasPid = true;
                _selmaker.PushBack1DVectors(s_stage4, vars, hasPid);
                break;
              case EFFICIENCY:
                _selmaker.PushBackEPVectors(s_stage4, vars); 
                break;
              case PDG:
                _selmaker.PushBackPIDVectors(s_stage4, vars);
                break;
              case TRACKCUTVAR:
                _selmaker.PushBackTrackCutVar(s_stage4, vars);
                break;

            }
           
/*
            if (vars->isSimulation == 0) 
                std::cout << "event " << vars->run << "." << vars->subrun << "." << vars->event << " x,z" << vars->vertex_x << ", " << vars->vertex_z << " nPFP: " << vars->nSelectedPfparticles <<  std::endl; 
*/
            if (entry != -1 && infile != nullptr && hasPid == true){
              //infile->GetEntry(entry);
              outfile->Fill();

              int thisEntry = 0;
              if (ewin != nullptr && _config.DoEventWeightMatching == true){
                int pastentry = thisEntry;
                _treehandler.PrepareTreeForSearching(ewin);
                thisEntry = _treehandler.FindEntryFromEvent(ewin, ewvars, vars->run, vars->subrun, vars->event, pastentry);
                //_treehandler.SetEWTreeVars(ewin, ewvars);
                _treehandler.PrepareTreeForWriting(ewin);
                ewin->GetEntry(thisEntry);
                ewout->Fill();
              }

            }

          }

        }

      }

    }

    thisMatrix.push_back(m_stage0);
    thisMatrix.push_back(m_stage1);
    thisMatrix.push_back(m_stage2);
    thisMatrix.push_back(m_stage3); 
    thisMatrix.push_back(m_stage4); 

    return thisMatrix;
  };

  std::vector<std::vector<std::vector<double>>> SelectionMaker::GetPlottingVariables(var_list* vars, kVarType var_type, TTree* infile,  TTree* outfile, int entry, ew_list* ewvars, TTree* ewin, TTree* ewout){

    numusel::AnalysisCuts _anacuts; 
    thisMatrix = GetPlottingVariables(vars, var_type, _anacuts.pid_cutvalue, infile, outfile, entry, ewvars, ewin, ewout);

    return thisMatrix;

  }

  void SelectionMaker::PushBackEPVectors(std::vector<std::vector<double>>* m_stagex, var_list* vars){

    // for genie mcparicles, the first entry is always the neutrino, and the fourth is the outgoing lepton
    m_stagex->push_back(std::vector<double>({(double)0.}));
    m_stagex->push_back(std::vector<double>({(double)vars->true_genie_starte->at(0)}));
    if (vars->true_genie_startp->size() >= 5)
      m_stagex->push_back(std::vector<double>({(double)vars->true_genie_startp->at(4)}));
    else m_stagex->push_back(std::vector<double>({-1.0}));

  }

  void SelectionMaker::PushBackPIDVectors(std::vector<std::vector<double>>* m_stagex, var_list* vars){

    m_stagex->push_back(std::vector<double>(vars->true_match_pdg->begin(),vars->true_match_pdg->end()));

  }

  void SelectionMaker::PushBackTrackCutVar(std::vector<std::vector<double>>* m_stagex, var_list* vars){

    std::vector<double> iscutpassed;

    // this currently just checks whether the pfparticles were identified
    // as a shower by pandora
    for (int i = 0; i < vars->pfp_pdgCode->size(); i++){

      double val = -1;
      if (std::abs(vars->pfp_pdgCode->at(i)) == 13)
        val = 1;
      else val = 0;

      iscutpassed.push_back(val);

    }
    m_stagex->push_back(iscutpassed);

  }


  void SelectionMaker::PushBack1DVectors(std::vector<std::vector<double>>* m_stagex, var_list* vars, bool isHasPID){

    // n.b. every time you add a variable here you need to
    // add the histogram name, title, and bin ranges in HistogramHandler.h
    // vectors

    numusel::AnalysisCuts     _anacuts;
    numusel::Configuration    _config;

    TF1* proton_range_energy_func = new TF1("proton_range_energy_func", "(x-[1])/[0]", 0, 10);
    TF1* muon_range_energy_func = new TF1("muon_range_energy_func", "(x-[1])/[0]", 0, 10);
    TF1* muon_mcs_energy_func = new TF1("muon_mcs_energy_func", "(x-[1])/[0]", 0, 10);

    proton_range_energy_func->SetParameters(_config.proton_range_m, _config.proton_range_c);
    muon_range_energy_func->SetParameters(_config.muon_range_contained_m, _config.muon_range_contained_c);
    muon_mcs_energy_func->SetParameters(_config.muon_mcs_uncontained_m, _config.muon_mcs_uncontained_c);

    m_stagex->push_back(std::vector<double>({(double)0.}));
    m_stagex->push_back(std::vector<double>({(double)vars->nSelectedTracks}));
    m_stagex->push_back(std::vector<double>({(double)vars->nSelectedShowers}));
    m_stagex->push_back(std::vector<double>({(double)vars->nSelectedPfparticles}));
    m_stagex->push_back(std::vector<double>({(double)vars->vertex_x}));
    m_stagex->push_back(std::vector<double>({(double)vars->vertex_y}));
    m_stagex->push_back(std::vector<double>({(double)vars->vertex_z}));
    m_stagex->push_back(*vars->track_length);
    m_stagex->push_back(*vars->track_startx);
    m_stagex->push_back(*vars->track_endx);
    m_stagex->push_back(*vars->track_starty);
    m_stagex->push_back(*vars->track_endy);
    m_stagex->push_back(*vars->track_startz);
    m_stagex->push_back(*vars->track_endz);
    m_stagex->push_back(*vars->track_theta);
    m_stagex->push_back(*vars->track_costheta);
    m_stagex->push_back(*vars->track_phi);

    std::vector<double> pid;
    for (int i = 0; i < vars->noBragg_fwd_mip->size(); i++){
      double lmip = vars->noBragg_fwd_mip->at(i);
      double lp = std::max(vars->bragg_fwd_p->at(i), vars->bragg_bwd_p->at(i));

      if (lmip == -999)
        pid.push_back(-999);
      else 
      {
        pid.push_back(std::log(lmip/lp));
        //if (std::log(lmip/lp) > 1 && std::log(lmip/lp) < 3 && vars->pfp_pdgCode->at(i) == 11 && vars->true_match_pdg->at(i) == 2212){
        //std::cout << vars->run << "." << vars->subrun << "." << vars->event << "   log(lmip/lp): " << log(lmip/lp) << " track start z: " << vars->track_startz->at(i) << " endz: " << vars->track_endz->at(i) << " truematchpurity: " << vars->true_match_purity->at(i) << std::endl;
        //TFile *out_dedxs = new TFile("out_dedxs.root", "update");
        //TGraph* h_tmp = new TGraph(vars->track_dedxperhit_smeared->at(i).size(), &vars->track_resrangeperhit->at(i)[0], &vars->track_dedxperhit_smeared->at(i)[0]);
        //TString name = Form("run%isubrun%ievent%i", vars->run, vars->subrun, vars->event);
        //h_tmp->SetName(name);
        //if (out_dedxs->Get(name)) continue;
        //h_tmp->Write();
        //out_dedxs->Close();
        //}
      }
    }
    m_stagex->push_back(pid);

    std::vector<double> depErangeE_mu_Y = {};
    std::vector<double> depErangeE_p_Y = {};

    for (int i = 0; i < vars->track_range_energy_muassumption->size(); i++){
       depErangeE_mu_Y.push_back(vars->track_dep_energy_yplane->at(i) - vars->track_range_energy_muassumption->at(i)); 
       depErangeE_p_Y.push_back(vars->track_dep_energy_yplane->at(i) - vars->track_range_energy_passumption->at(i)); 
    }

    m_stagex->push_back(*vars->track_mcs_muassmp_fwd);
    m_stagex->push_back(*vars->track_mcs_muassmp_bwd);
    m_stagex->push_back(*vars->track_mcs_muassmp_energy_fwd);
    m_stagex->push_back(*vars->track_mcs_muassmp_energy_bwd);
    m_stagex->push_back(*vars->track_range_mom_muassumption);
    m_stagex->push_back(*vars->track_range_mom_passumption);
    m_stagex->push_back(*vars->track_range_energy_muassumption);
    m_stagex->push_back(*vars->track_range_energy_passumption);
    m_stagex->push_back(*vars->track_dep_energy_yplane);
    m_stagex->push_back(depErangeE_mu_Y);
    m_stagex->push_back(depErangeE_p_Y);
    m_stagex->push_back(*vars->track_residualrms);

    if (isHasPID == true){

      // candidate muon variables
      std::vector<double> candMuonLength;
      std::vector<double> candMuonTheta;
      std::vector<double> candMuonCosTheta;
      std::vector<double> candMuonPhi;
      std::vector<double> candMuonMCSFwd;
      std::vector<double> candMuonMCSBwd;
      std::vector<double> candMuonMCSFwdEnergy;
      std::vector<double> candMuonMCSBwdEnergy;
      std::vector<double> candMuonRangeMomentumMuassmp;
      std::vector<double> candMuonRangeEnergyMuassmp;
      std::vector<double> candMuonEdep;
      std::vector<double> candMuonEdepMinusErange_mu;
      std::vector<double> candMuonEdepMinusErange_p;
      std::vector<double> candMuonLength_contained;
      std::vector<double> candMuonTheta_contained;
      std::vector<double> candMuonCosTheta_contained;
      std::vector<double> candMuonPhi_contained;
      std::vector<double> candMuonMCSFwd_contained;
      std::vector<double> candMuonMCSBwd_contained;
      std::vector<double> candMuonMCSFwdEnergy_contained;
      std::vector<double> candMuonMCSBwdEnergy_contained;
      std::vector<double> candMuonRangeMomentumMuassmp_contained;
      std::vector<double> candMuonRangeEnergyMuassmp_contained;
      std::vector<double> candMuonEdep_contained;
      std::vector<double> candMuonEdepMinusErange_mu_contained;
      std::vector<double> candMuonEdepMinusErange_p_contained;
      std::vector<double> candMuonLength_uncontained;
      std::vector<double> candMuonTheta_uncontained;
      std::vector<double> candMuonCosTheta_uncontained;
      std::vector<double> candMuonPhi_uncontained;
      std::vector<double> candMuonMCSFwd_uncontained;
      std::vector<double> candMuonMCSBwd_uncontained;
      std::vector<double> candMuonMCSFwdEnergy_uncontained;
      std::vector<double> candMuonMCSBwdEnergy_uncontained;
      std::vector<double> candMuonRangeMomentumMuassmp_uncontained;
      std::vector<double> candMuonRangeEnergyMuassmp_uncontained;
      std::vector<double> candMuonEdep_uncontained;
      std::vector<double> candMuonEdepMinusErange_mu_uncontained;
      std::vector<double> candMuonEdepMinusErange_p_uncontained;

      // candidate proton variables
      std::vector<double> candProtonLength;
      std::vector<double> candProtonTheta;
      std::vector<double> candProtonCosTheta;
      std::vector<double> candProtonPhi;
      std::vector<double> candProtonMCSFwd;
      std::vector<double> candProtonMCSBwd;
      std::vector<double> candProtonRangeMomentumPassmp;
      std::vector<double> candProtonRangeEnergyPassmp;
      std::vector<double> candProtonEdep;
      std::vector<double> candProtonEdepMinusErange_mu;
      std::vector<double> candProtonEdepMinusErange_p;
      std::vector<double> candLeadingProtonLength;
      std::vector<double> candLeadingProtonTheta;
      std::vector<double> candLeadingProtonCosTheta;
      std::vector<double> candLeadingProtonPhi;
      std::vector<double> candLeadingProtonMCSFwd;
      std::vector<double> candLeadingProtonMCSBwd;
      std::vector<double> candLeadingProtonRangeMomentumPassmp;
      std::vector<double> candLeadingProtonRangeEnergyPassmp;
      std::vector<double> candLeadingProtonEdep;
      std::vector<double> candLeadingProtonEdepMinusErange_mu;
      std::vector<double> candLeadingProtonEdepMinusErange_p;
      std::vector<double> candNonLeadingProtonLength;
      std::vector<double> candNonLeadingProtonTheta;
      std::vector<double> candNonLeadingProtonCosTheta;
      std::vector<double> candNonLeadingProtonPhi;
      std::vector<double> candNonLeadingProtonMCSFwd;
      std::vector<double> candNonLeadingProtonMCSBwd;
      std::vector<double> candNonLeadingProtonRangeMomentumPassmp;
      std::vector<double> candNonLeadingProtonRangeEnergyPassmp;
      std::vector<double> candNonLeadingProtonEdep;
      std::vector<double> candNonLeadingProtonEdepMinusErange_mu;
      std::vector<double> candNonLeadingProtonEdepMinusErange_p;
      std::vector<double> neutrinoReconstructedEnergyUncalib = {0.};
      std::vector<double> neutrinoReconstructedEnergyCalib = {0.};

      std::vector<std::pair<int, float>> leadingProtonFinder;

      vars->track_ismuoncand->resize(0);
      vars->track_isprotoncand->resize(0);

      for (int i = 0; i < pid.size(); i++){

        // get muon candidate
        // if the pid value is 0, then this indicated there was no calorimetry information
        // in the collection plane
        if (pid.at(i) > _anacuts.pid_cutvalue){

          vars->track_ismuoncand->push_back(1);
          vars->track_isprotoncand->push_back(0);

          candMuonLength.push_back(vars->track_length->at(i));
          candMuonTheta.push_back(vars->track_theta->at(i));
          candMuonCosTheta.push_back(vars->track_costheta->at(i));
          candMuonPhi.push_back(vars->track_phi->at(i));
          candMuonMCSFwd.push_back(vars->track_mcs_muassmp_fwd->at(i));
          candMuonMCSBwd.push_back(vars->track_mcs_muassmp_bwd->at(i));
          candMuonMCSFwdEnergy.push_back(vars->track_mcs_muassmp_energy_fwd->at(i));
          candMuonMCSBwdEnergy.push_back(vars->track_mcs_muassmp_energy_bwd->at(i));
          candMuonRangeMomentumMuassmp.push_back(vars->track_range_mom_muassumption->at(i));
          candMuonRangeEnergyMuassmp.push_back(vars->track_range_energy_muassumption->at(i));
          candMuonEdep.push_back(vars->track_dep_energy_yplane->at(i));
          candMuonEdepMinusErange_mu.push_back(depErangeE_mu_Y.at(i));
          candMuonEdepMinusErange_p.push_back(depErangeE_p_Y.at(i));

          if (vars->track_isContained->at(i) == true){

            candMuonLength_contained.push_back(vars->track_length->at(i));
            candMuonTheta_contained.push_back(vars->track_theta->at(i));
            candMuonCosTheta_contained.push_back(vars->track_costheta->at(i));
            candMuonPhi_contained.push_back(vars->track_phi->at(i));
            candMuonMCSFwd_contained.push_back(vars->track_mcs_muassmp_fwd->at(i));
            candMuonMCSBwd_contained.push_back(vars->track_mcs_muassmp_bwd->at(i));
            candMuonMCSFwdEnergy_contained.push_back(vars->track_mcs_muassmp_energy_fwd->at(i));
            candMuonMCSBwdEnergy_contained.push_back(vars->track_mcs_muassmp_energy_bwd->at(i));
            candMuonRangeMomentumMuassmp_contained.push_back(vars->track_range_mom_muassumption->at(i));
            candMuonRangeEnergyMuassmp_contained.push_back(vars->track_range_energy_muassumption->at(i));
            candMuonEdep_contained.push_back(vars->track_dep_energy_yplane->at(i));
            candMuonEdepMinusErange_mu_contained.push_back(depErangeE_mu_Y.at(i));
            candMuonEdepMinusErange_p_contained.push_back(depErangeE_p_Y.at(i));

            neutrinoReconstructedEnergyUncalib.at(0)+=vars->track_range_energy_muassumption->at(i);
            neutrinoReconstructedEnergyCalib.at(0)+=muon_range_energy_func->Eval(vars->track_range_energy_muassumption->at(i));
            
          }
          else{
            candMuonLength_uncontained.push_back(vars->track_length->at(i));
            candMuonTheta_uncontained.push_back(vars->track_theta->at(i));
            candMuonCosTheta_uncontained.push_back(vars->track_costheta->at(i));
            candMuonPhi_uncontained.push_back(vars->track_phi->at(i));
            candMuonMCSFwd_uncontained.push_back(vars->track_mcs_muassmp_fwd->at(i));
            candMuonMCSBwd_uncontained.push_back(vars->track_mcs_muassmp_bwd->at(i));
            candMuonMCSFwdEnergy_uncontained.push_back(vars->track_mcs_muassmp_energy_fwd->at(i));
            candMuonMCSBwdEnergy_uncontained.push_back(vars->track_mcs_muassmp_energy_bwd->at(i));
            candMuonRangeMomentumMuassmp_uncontained.push_back(vars->track_range_mom_muassumption->at(i));
            candMuonRangeEnergyMuassmp_uncontained.push_back(vars->track_range_energy_muassumption->at(i));
            candMuonEdep_uncontained.push_back(vars->track_dep_energy_yplane->at(i));
            candMuonEdepMinusErange_mu_uncontained.push_back(depErangeE_mu_Y.at(i));
            candMuonEdepMinusErange_p_uncontained.push_back(depErangeE_p_Y.at(i));

            float mcs_energy = 0;

            if (vars->track_mcs_muassmp_fwd_loglikelihood->at(i) < vars->track_mcs_muassmp_bwd_loglikelihood->at(i)){
              mcs_energy = vars->track_mcs_muassmp_energy_fwd->at(i);
            }
            else {
              mcs_energy = vars->track_mcs_muassmp_energy_bwd->at(i);
            }
            neutrinoReconstructedEnergyUncalib.at(0)+=mcs_energy;
            neutrinoReconstructedEnergyCalib.at(0)+=muon_mcs_energy_func->Eval(mcs_energy);

          }

        }
        // else they're proton candidates
        else{

          vars->track_ismuoncand->push_back(0);
          vars->track_isprotoncand->push_back(1);

          candProtonLength.push_back(vars->track_length->at(i));
          candProtonTheta.push_back(vars->track_theta->at(i));
          candProtonCosTheta.push_back(vars->track_costheta->at(i));
          candProtonPhi.push_back(vars->track_phi->at(i));
          candProtonMCSFwd.push_back(vars->track_mcs_muassmp_fwd->at(i));
          candProtonMCSBwd.push_back(vars->track_mcs_muassmp_bwd->at(i));
          candProtonRangeMomentumPassmp.push_back(vars->track_range_mom_passumption->at(i));
          candProtonRangeEnergyPassmp.push_back(vars->track_range_energy_passumption->at(i));
          candProtonEdep.push_back(vars->track_dep_energy_yplane->at(i));
          candProtonEdepMinusErange_mu.push_back(depErangeE_mu_Y.at(i));
          candProtonEdepMinusErange_p.push_back(depErangeE_p_Y.at(i));

          std::pair<int, float> thisProtonInformation;
          thisProtonInformation.first = i;
          thisProtonInformation.second = vars->track_length->at(i);
          leadingProtonFinder.push_back(thisProtonInformation);

          neutrinoReconstructedEnergyUncalib.at(0)+=vars->track_range_energy_passumption->at(i);
          neutrinoReconstructedEnergyCalib.at(0)+=proton_range_energy_func->Eval(vars->track_range_energy_passumption->at(i));

        }

      }

      // sort thisProtonInformation based on length
      std::sort(leadingProtonFinder.begin(), leadingProtonFinder.end(), [](auto &left, auto &right){
          return left.second > right.second;
          });

      if (leadingProtonFinder.size() > 0){
        // get leading proton information
        candLeadingProtonLength.push_back(vars->track_length->at(leadingProtonFinder.at(0).first));
        candLeadingProtonTheta.push_back(vars->track_theta->at(leadingProtonFinder.at(0).first));
        candLeadingProtonCosTheta.push_back(vars->track_costheta->at(leadingProtonFinder.at(0).first));
        candLeadingProtonPhi.push_back(vars->track_phi->at(leadingProtonFinder.at(0).first));
        candLeadingProtonMCSFwd.push_back(vars->track_mcs_muassmp_fwd->at(leadingProtonFinder.at(0).first));
        candLeadingProtonMCSBwd.push_back(vars->track_mcs_muassmp_bwd->at(leadingProtonFinder.at(0).first));
        candLeadingProtonRangeMomentumPassmp.push_back(vars->track_range_mom_passumption->at(leadingProtonFinder.at(0).first));
        candLeadingProtonRangeEnergyPassmp.push_back(vars->track_range_energy_passumption->at(leadingProtonFinder.at(0).first));
        candLeadingProtonEdep.push_back(vars->track_dep_energy_yplane->at(leadingProtonFinder.at(0).first));
        candLeadingProtonEdepMinusErange_mu.push_back(depErangeE_mu_Y.at(leadingProtonFinder.at(0).first));
        candLeadingProtonEdepMinusErange_p.push_back(depErangeE_p_Y.at(leadingProtonFinder.at(0).first));
      }

      //and non-leading proton information
      for (int i = 1; i <leadingProtonFinder.size(); i++){

        candNonLeadingProtonLength.push_back(vars->track_length->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonTheta.push_back(vars->track_theta->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonCosTheta.push_back(vars->track_costheta->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonPhi.push_back(vars->track_phi->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonMCSFwd.push_back(vars->track_mcs_muassmp_fwd->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonMCSBwd.push_back(vars->track_mcs_muassmp_bwd->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonRangeMomentumPassmp.push_back(vars->track_range_mom_passumption->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonRangeEnergyPassmp.push_back(vars->track_range_energy_passumption->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonEdep.push_back(vars->track_dep_energy_yplane->at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonEdepMinusErange_mu.push_back(depErangeE_mu_Y.at(leadingProtonFinder.at(i).first));
        candNonLeadingProtonEdepMinusErange_p.push_back(depErangeE_p_Y.at(leadingProtonFinder.at(i).first));
      }

      m_stagex->push_back(candMuonLength);
      m_stagex->push_back(candMuonTheta);
      m_stagex->push_back(candMuonCosTheta);
      m_stagex->push_back(candMuonPhi);
      m_stagex->push_back(candMuonMCSFwd);
      m_stagex->push_back(candMuonMCSBwd);
      m_stagex->push_back(candMuonMCSFwdEnergy);
      m_stagex->push_back(candMuonMCSBwdEnergy);
      m_stagex->push_back(candMuonRangeMomentumMuassmp);
      m_stagex->push_back(candMuonRangeEnergyMuassmp);
      m_stagex->push_back(candMuonEdep);
      m_stagex->push_back(candMuonEdepMinusErange_mu);
      m_stagex->push_back(candMuonEdepMinusErange_p);
      m_stagex->push_back(candMuonLength_contained);
      m_stagex->push_back(candMuonTheta_contained);
      m_stagex->push_back(candMuonCosTheta_contained);
      m_stagex->push_back(candMuonPhi_contained);
      m_stagex->push_back(candMuonMCSFwd_contained);
      m_stagex->push_back(candMuonMCSBwd_contained);
      m_stagex->push_back(candMuonMCSFwdEnergy_contained);
      m_stagex->push_back(candMuonMCSBwdEnergy_contained);
      m_stagex->push_back(candMuonRangeMomentumMuassmp_contained);
      m_stagex->push_back(candMuonRangeEnergyMuassmp_contained);
      m_stagex->push_back(candMuonEdep_contained);
      m_stagex->push_back(candMuonEdepMinusErange_mu_contained);
      m_stagex->push_back(candMuonEdepMinusErange_p_contained);
      m_stagex->push_back(candMuonLength_uncontained);
      m_stagex->push_back(candMuonTheta_uncontained);
      m_stagex->push_back(candMuonCosTheta_uncontained);
      m_stagex->push_back(candMuonPhi_uncontained);
      m_stagex->push_back(candMuonMCSFwd_uncontained);
      m_stagex->push_back(candMuonMCSBwd_uncontained);
      m_stagex->push_back(candMuonMCSFwdEnergy_uncontained);
      m_stagex->push_back(candMuonMCSBwdEnergy_uncontained);
      m_stagex->push_back(candMuonRangeMomentumMuassmp_uncontained);
      m_stagex->push_back(candMuonRangeEnergyMuassmp_uncontained);
      m_stagex->push_back(candMuonEdep_uncontained);
      m_stagex->push_back(candMuonEdepMinusErange_mu_uncontained);
      m_stagex->push_back(candMuonEdepMinusErange_p_uncontained);
      m_stagex->push_back(candProtonLength);
      m_stagex->push_back(candProtonTheta);
      m_stagex->push_back(candProtonCosTheta);
      m_stagex->push_back(candProtonPhi);
      m_stagex->push_back(candProtonMCSFwd);
      m_stagex->push_back(candProtonMCSBwd);
      m_stagex->push_back(candProtonRangeMomentumPassmp);
      m_stagex->push_back(candProtonRangeEnergyPassmp);
      m_stagex->push_back(candProtonEdep);
      m_stagex->push_back(candProtonEdepMinusErange_mu);
      m_stagex->push_back(candProtonEdepMinusErange_p);
      m_stagex->push_back(candLeadingProtonLength);
      m_stagex->push_back(candLeadingProtonTheta);
      m_stagex->push_back(candLeadingProtonCosTheta);
      m_stagex->push_back(candLeadingProtonPhi);
      m_stagex->push_back(candLeadingProtonMCSFwd);
      m_stagex->push_back(candLeadingProtonMCSBwd);
      m_stagex->push_back(candLeadingProtonRangeMomentumPassmp);
      m_stagex->push_back(candLeadingProtonRangeEnergyPassmp);
      m_stagex->push_back(candLeadingProtonEdep);
      m_stagex->push_back(candLeadingProtonEdepMinusErange_mu);
      m_stagex->push_back(candLeadingProtonEdepMinusErange_p);
      m_stagex->push_back(candNonLeadingProtonLength);
      m_stagex->push_back(candNonLeadingProtonTheta);
      m_stagex->push_back(candNonLeadingProtonCosTheta);
      m_stagex->push_back(candNonLeadingProtonPhi);
      m_stagex->push_back(candNonLeadingProtonMCSFwd);
      m_stagex->push_back(candNonLeadingProtonMCSBwd);
      m_stagex->push_back(candNonLeadingProtonRangeMomentumPassmp);
      m_stagex->push_back(candNonLeadingProtonRangeEnergyPassmp);
      m_stagex->push_back(candNonLeadingProtonEdep);
      m_stagex->push_back(candNonLeadingProtonEdepMinusErange_mu);
      m_stagex->push_back(candNonLeadingProtonEdepMinusErange_p);
      m_stagex->push_back(neutrinoReconstructedEnergyUncalib);
      m_stagex->push_back(neutrinoReconstructedEnergyCalib);

      vars->reconstructedNeutrinoEnergy = neutrinoReconstructedEnergyUncalib.at(0);
      vars->reconstructedNeutrinoEnergyCalib = neutrinoReconstructedEnergyCalib.at(0);

    }

    proton_range_energy_func->Delete();
    muon_range_energy_func->Delete();
    muon_mcs_energy_func->Delete();

  }

  void SelectionMaker::PushBack2DVectors(std::vector<std::vector<std::pair<double, double>>>* m_stagex, var_list* vars, bool isHasPID){

    // n.b. every time you add a variable here you need to
    // add the histogram name, title, and bin ranges in HistogramHandler.h
    // vectors

    numusel::AnalysisCuts _anacuts;
    numusel::Configuration _config;

    TF1* proton_range_energy_func = new TF1("proton_range_energy_func", "(x-[1])/[0]", 0, 10);
    TF1* muon_range_energy_func = new TF1("muon_range_energy_func", "(x-[1])/[0]", 0, 10);
    TF1* muon_mcs_energy_func = new TF1("muon_mcs_energy_func", "(x-[1])/[0]", 0, 10);

    proton_range_energy_func->SetParameters(_config.proton_range_m, _config.proton_range_c);
    muon_range_energy_func->SetParameters(_config.muon_range_contained_m, _config.muon_range_contained_c);
    muon_mcs_energy_func->SetParameters(_config.muon_mcs_uncontained_m, _config.muon_mcs_uncontained_c);

    std::vector<std::pair<double, double>> track_thetaphi;
    for (int i = 0; i < vars->track_theta->size(); i++){ 
      std::pair<double, double> track_thetaphipair(vars->track_theta->at(i), vars->track_phi->at(i));
      track_thetaphi.push_back(track_thetaphipair);
    }

    m_stagex->push_back(track_thetaphi);

    if (isHasPID == true){

      std::vector<std::pair<double, double>> candMuondEdxResRange;
      std::vector<std::pair<double, double>> candMuondEdxResRange_contained;
      std::vector<std::pair<double, double>> candMuondEdxResRange_uncontained;
      std::vector<std::pair<double, double>> candProtondEdxResRange;
      std::vector<std::pair<double, double>> candProtondEdxResRange_leading;
      std::vector<std::pair<double, double>> candProtondEdxResRange_nonleading;
      std::vector<std::pair<double, double>> neutrinoReconstructedEnergyCalib;

      std::pair<double, double> neutrinoReconstructedEnergyCalib_tmp(0., 0.);

      // setup pid vector for later use
      std::vector<double> pid;
      for (int i = 0; i < vars->noBragg_fwd_mip->size(); i++){
        double lmip = vars->noBragg_fwd_mip->at(i);
        double lp = std::max(vars->bragg_fwd_p->at(i), vars->bragg_bwd_p->at(i));

        if (lmip == -999)
          pid.push_back(-999);
        else pid.push_back(std::log(lmip/lp));
      }

      std::vector<std::pair<int, float>> leadingProtonFinder;
      leadingProtonFinder.resize(0);

      // loop PID vector that we've just filled
      // We want to find:
      // * contained muons
      // * uncontained muons
      // * protons
      // separately since they require different energy reconstructions
      //
      // Also going to separate protons out by leading/non leading, mostly for 
      // interests sake

      for (int i = 0; i < pid.size(); i++){

        // get muon candidate
        if (pid.at(i) > _anacuts.pid_cutvalue){

          // any muon
          for (int j = 0; j < vars->track_dedxperhit_smeared->at(i).size(); j++){
            std::pair<float, float> thisPair(vars->track_resrangeperhit->at(i).at(j), vars->track_dedxperhit_smeared->at(i).at(j));
            candMuondEdxResRange.push_back(thisPair);
          }

          // contained muons
          if (vars->track_isContained->at(i) == true){

            for (int j = 0; j < vars->track_dedxperhit_smeared->at(i).size(); j++){
              std::pair<float, float> thisPair(vars->track_resrangeperhit->at(i).at(j), vars->track_dedxperhit_smeared->at(i).at(j));
              candMuondEdxResRange_contained.push_back(thisPair);
            }

            // add reconstructed energy for contained muons
            neutrinoReconstructedEnergyCalib_tmp.second += muon_range_energy_func->Eval(vars->track_range_energy_muassumption->at(i));

          }
          // uncontained muons
          else{

            for (int j = 0; j < vars->track_dedxperhit_smeared->at(i).size(); j++){
              std::pair<float, float> thisPair(vars->track_resrangeperhit->at(i).at(j), vars->track_dedxperhit_smeared->at(i).at(j));
              candMuondEdxResRange_uncontained.push_back(thisPair);
            }

            // get MCS energy
            float mcs_energy = 0;
            if (vars->track_mcs_muassmp_fwd_loglikelihood->at(i) < vars->track_mcs_muassmp_bwd_loglikelihood->at(i)){
              mcs_energy = vars->track_mcs_muassmp_energy_fwd->at(i);
            }
            else {
              mcs_energy = vars->track_mcs_muassmp_energy_bwd->at(i);
            }

            // add reconstructed energy from MCS for uncontained muons
            neutrinoReconstructedEnergyCalib_tmp.second += muon_mcs_energy_func->Eval(mcs_energy);
          }

        }
        // protons
        else{
          if (pid.at(i) < _anacuts.pid_cutvalue){

            for (int j = 0; j < vars->track_dedxperhit_smeared->at(i).size(); j++){
              std::pair<float, float> thisPair(vars->track_resrangeperhit->at(i).at(j), vars->track_dedxperhit_smeared->at(i).at(j));

              candProtondEdxResRange.push_back(thisPair);
            }

            std::pair<int, float> thisProtonInformation;
            thisProtonInformation.first = i;
            thisProtonInformation.second = vars->track_length->at(i);
            leadingProtonFinder.push_back(thisProtonInformation);

          }

          // add reconstructed energy from range for protons
          neutrinoReconstructedEnergyCalib_tmp.second += proton_range_energy_func->Eval(vars->track_range_energy_passumption->at(i));

        }

      }

      // now going to use the vector of pairs we created earlier, 
      // sort them by length, and then we can grab the leading
      // and non-leading protons
      if (leadingProtonFinder.size() > 0){

        // sort thisProtonInformation based on length
        std::sort(leadingProtonFinder.begin(), leadingProtonFinder.end(), [](auto &left, auto &right){
            return left.second > right.second;
            });

        // get leading proton information
        for (int j = 0; j < vars->track_dedxperhit_smeared->at(leadingProtonFinder.at(0).first).size(); j++){
          std::pair<float, float> thisPair(vars->track_resrangeperhit->at(leadingProtonFinder.at(0).first).at(j), vars->track_dedxperhit_smeared->at(leadingProtonFinder.at(0).first).at(j));
          candProtondEdxResRange_leading.push_back(thisPair);
        }

        //and non-leading proton information
        for (int i = 1; i < leadingProtonFinder.size(); i++){
          for (int j = 0; j < vars->track_dedxperhit_smeared->at(leadingProtonFinder.at(i).first).size(); j++){
            std::pair<float, float> thisPair(vars->track_resrangeperhit->at(leadingProtonFinder.at(i).first).at(j), vars->track_dedxperhit_smeared->at(leadingProtonFinder.at(i).first).at(j));
            candProtondEdxResRange_nonleading.push_back(thisPair);
          }

        }
      }

      // we can only make the true neutrino energy versus deposited energy for 
      // simulation, obviously
      if (vars->isSimulation)
        neutrinoReconstructedEnergyCalib_tmp.first = vars->true_genie_starte->at(0);
      else 
        neutrinoReconstructedEnergyCalib_tmp.first = 0;
      neutrinoReconstructedEnergyCalib.push_back(neutrinoReconstructedEnergyCalib_tmp);

      m_stagex->push_back(candMuondEdxResRange);
      m_stagex->push_back(candMuondEdxResRange_contained);
      m_stagex->push_back(candMuondEdxResRange_uncontained);
      m_stagex->push_back(candProtondEdxResRange);
      m_stagex->push_back(candProtondEdxResRange_leading);
      m_stagex->push_back(candProtondEdxResRange_nonleading);
      m_stagex->push_back(neutrinoReconstructedEnergyCalib);


    }

    proton_range_energy_func->Delete();
    muon_range_energy_func->Delete();
    muon_mcs_energy_func->Delete();

  }

}
