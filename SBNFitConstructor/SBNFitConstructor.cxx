#include "SBNFitConstructor.h"

int main(int argv, char** argc){

  // read in selected events trees
  TFile *f_selevents = new TFile("/uboone/app/users/alister1/numuSelection/ubcode_v06_26_01_18/srcs/uboonecode/uboone/NuMuSelection/AnalysisFramework/selectedEvents_sim.root", "read");
  TTree* event_info_tree = (TTree*)f_selevents->Get("simulation");
  TTree* ew_tree = (TTree*)f_selevents->Get("ew");
  initialiseTrees(event_info_tree, ew_tree);

  TFile * f_out = new TFile("numu_SBNFit.root", "RECREATE");
  TTree* t_cosmic = new TTree("numu_sbnfit_cosmic", "numu sbnfit tree cosmic");
  TTree* t_mixed = new TTree("numu_sbnfit_mixed", "numu sbnfit tree mixed");
  TTree* t_oofv = new TTree("numu_sbnfit_oofv", "numu sbnfit tree oofv");
  TTree* t_nc = new TTree("numu_sbnfit_nc", "numu sbnfit tree nc");
  TTree* t_nue = new TTree("numu_sbnfit_nue", "numu sbnfit tree nue/anue");
  TTree* t_anumu = new TTree("numu_sbnfit_anumu", "numu sbnfit tree anumu");
  TTree* t_ccother = new TTree("numu_sbnfit_ccother", "numu sbnfit tree ccother");
  TTree* t_cc0pinp = new TTree("numu_sbnfit_cc0pinp", "numu sbnfit tree cc0pinp");

  setVariables(t_cosmic);
  setVariables(t_mixed);
  setVariables(t_oofv);
  setVariables(t_nc);
  setVariables(t_nue);
  setVariables(t_anumu);
  setVariables(t_ccother);
  setVariables(t_cc0pinp);

  // input files
  std::vector<std::string> file_list;
  std::string file_name;
  std::ifstream input_file(argc[1]);

  while (getline(input_file, file_name)){
    file_list.push_back(file_name);
  }

  for (gallery::Event ev(file_list); !ev.atEnd(); ev.next()){

    /**
     * Ideally, we do the matching by run, subrun, event IDs, but 
     * for some reason the event numbers are just iterated when constructing
     * the principal through TreeReader
     * 
     * I'm sorry for the way I'm doing the matching right now
     */
    int run = ev.eventAuxiliary().runID().run();
    int sub_run = ev.eventAuxiliary().subRunID().subRun();
    int event = ev.eventAuxiliary().event();


    const auto& mctruthHandle = ev.getValidHandle< std::vector<simb::MCTruth> > ("TreeReader");
    std::vector<evwgh::MCEventWeight> eventweightHandle = *ev.getValidHandle< std::vector<evwgh::MCEventWeight> > ("mcweight");

    evwgh::MCEventWeight evtwght = eventweightHandle.at(0);

    for (simb::MCTruth const& mctruth : (*mctruthHandle)){
      double nu_w = mctruth.GetNeutrino().W();
      double nu_x = mctruth.GetNeutrino().X();
      double nu_y = mctruth.GetNeutrino().Y();
      double nu_qsqr = mctruth.GetNeutrino().QSqr();

      for (int i = 0; i < ew_tree->GetEntries(); i++){

        ew_tree->GetEntry(i);

        if (nu_w == ew_nu_w && nu_x == ew_nu_x && nu_y == ew_nu_y && nu_qsqr == ew_nu_qsqr){

          event_info_tree->GetEntry(i);
 
          out_reconstructed_neutrino_energy_calib = sel_reconstructed_neutrino_energy_calib;
          weights = evtwght.fWeight; 
          
          if (eventCat->at(7) == 1)
            t_cosmic->Fill();
          
          else if (eventCat->at(6) == 1)
            t_mixed->Fill();

          else if (eventCat->at(5) == 1)
            t_oofv->Fill();      

          else if (eventCat->at(4) == 1)
            t_nc->Fill();      
          
          else if (eventCat->at(3) == 1)
            t_nue->Fill();      
          
          else if (eventCat->at(2) == 1)
            t_anumu->Fill();      

          else if (eventCat->at(1) == 1)
            t_ccother->Fill();      
          
          else if (eventCat->at(0) == 1)
            t_cc0pinp->Fill();      
         
        }

      }

    }


  }

  f_out->Write();
}

void initialiseTrees(TTree* sel, TTree* ew){

  sel->SetBranchStatus("*", 0);
  sel->SetBranchStatus("run", 1);
  sel->SetBranchStatus("subrun", 1);
  sel->SetBranchStatus("event", 1);
  sel->SetBranchStatus("reconstructed_neutrino_energy_calib", 1);
  sel->SetBranchStatus("event_cat", 1);

  sel->SetBranchAddress("run", &sel_run);
  sel->SetBranchAddress("subrun", &sel_subrun);
  sel->SetBranchAddress("event", &sel_event);
  sel->SetBranchAddress("reconstructed_neutrino_energy_calib", &sel_reconstructed_neutrino_energy_calib);
  sel->SetBranchAddress("event_cat", &eventCat);

  ew->SetBranchStatus("*", 0);
  ew->SetBranchStatus("run", 1);
  ew->SetBranchStatus("subrun", 1);
  ew->SetBranchStatus("event", 1);
  ew->SetBranchStatus("MCTruth_neutrino_W", 1);
  ew->SetBranchStatus("MCTruth_neutrino_X", 1);
  ew->SetBranchStatus("MCTruth_neutrino_Y", 1);
  ew->SetBranchStatus("MCTruth_neutrino_QSqr", 1);

  ew->SetBranchAddress("run", &ew_run);
  ew->SetBranchAddress("subrun", &ew_subrun);
  ew->SetBranchAddress("event", &ew_event);
  ew->SetBranchAddress("MCTruth_neutrino_W", &ew_nu_w);
  ew->SetBranchAddress("MCTruth_neutrino_X", &ew_nu_x);
  ew->SetBranchAddress("MCTruth_neutrino_Y", &ew_nu_y);
  ew->SetBranchAddress("MCTruth_neutrino_QSqr", &ew_nu_qsqr);

}

void setVariables(TTree* t){
  t->Branch("reconstructed_neutrino_energy_calib", &out_reconstructed_neutrino_energy_calib);
  t->Branch("weights", "std::map<std::string, std::vector<double>>", &weights);
}
