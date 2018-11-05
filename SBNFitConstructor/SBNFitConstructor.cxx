#include "SBNFitConstructor.h"

int main(int argv, char** argc){

  // read in selected events trees
  TFile *f_selevents = new TFile("/uboone/data/users/alister1/testBuildShowersAsTracks_3/selectedEvents_sim.root", "read");
  TTree* event_info_tree = (TTree*)f_selevents->Get("simulation");
  TTree* ew_tree = (TTree*)f_selevents->Get("ew");
  initialiseTrees(event_info_tree, ew_tree);

  TFile * f_out = new TFile("numu_SBNFit.root", "RECREATE");
  TTree* t_out = new TTree("numu_sbnfit", "numu sbnfit tree");
  t_out->Branch("reconstructed_neutrino_energy", &out_reconstructed_neutrino_energy);
  t_out->Branch("weights", "std::map<std::string, std::vector<double>>", &weights);

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
          
          out_reconstructed_neutrino_energy = sel_resconstructed_neutrino_energy;
          weights = evtwght.fWeight; 
          t_out->Fill();
          
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
  sel->SetBranchStatus("reconstructed_neutrino_energy", 1);

  sel->SetBranchAddress("run", &sel_run);
  sel->SetBranchAddress("subrun", &sel_subrun);
  sel->SetBranchAddress("event", &sel_event);
  sel->SetBranchAddress("reconstructed_neutrino_energy", &sel_resconstructed_neutrino_energy);

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
