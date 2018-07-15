#include "MakeDataSimulationComparisonPlots.h"

int main(){

  // pull out TTrees from provided TFiles
  TTree* t_onbeam = (TTree*)(new TFile(s_onbeam.c_str(), "read"))->Get("numuselection/analysis_tree");
  TTree* t_offbeam = (TTree*)(new TFile(s_offbeam.c_str(), "read"))->Get("numuselection/analysis_tree");
  TTree* t_simulation = (TTree*)(new TFile(s_simulation.c_str(), "read"))->Get("numuselection/analysis_tree");

  // initialise variables and trees
  var_list onbeam_vars;
  var_list offbeam_vars;
  var_list simulation_vars;

  setTreeVars(t_onbeam, &onbeam_vars, false);
  setTreeVars(t_offbeam, &offbeam_vars, false);
  setTreeVars(t_simulation, &simulation_vars, true);

  // make hitograms to fill
  TH1D* h_nTracks_sim = new TH1D("h_nTracks_sim", ";Number of Tracks;", 5, 0, 5);
  TH1D* h_nShowers_sim = new TH1D("h_nShowers_sim", ";Number of Showers;", 5, 0, 5);
  TH1D* h_trackLength_sim = new TH1D("h_trackLength_sim", ";Track Length (cm);", 50, 0, 700);
  TH1D* h_nTracks_onbeam = new TH1D("h_nTracks_onbeam", ";Number of Tracks;", 5, 0, 5);
  TH1D* h_nShowers_onbeam = new TH1D("h_nShowers_onbeam", ";Number of Showers;", 5, 0, 5);
  TH1D* h_trackLength_onbeam = new TH1D("h_trackLength_onbeam", ";Track Length (cm);", 50, 0, 700);
  TH1D* h_nTracks_offbeam = new TH1D("h_nTracks_offbeam", ";Number of Tracks;", 5, 0, 5);
  TH1D* h_nShowers_offbeam = new TH1D("h_nShowers_offbeam", ";Number of Showers;", 5, 0, 5);
  TH1D* h_trackLength_offbeam = new TH1D("h_trackLength_offbeam", ";Track Length (cm);", 50, 0, 700);

  // loop simulation
  for (int i = 0; i < t_simulation->GetEntries(); i++){
    
    t_simulation->GetEntry(i);

    h_nTracks_sim->Fill(simulation_vars.nSelectedTracks);
    h_nShowers_sim->Fill(simulation_vars.nSelectedShowers);
   // h_trackLength_sim->Fill(simulation_vars.track_length);

  }

  h_nTracks_sim->Sumw2();
  h_nShowers_sim->Sumw2();
  h_trackLength_sim->Sumw2();
  h_nTracks_sim->Scale(simscaling);
  h_nShowers_sim->Scale(simscaling);
  h_trackLength_sim->Scale(simscaling);
  
  // loop onbeam
  for (int i = 0; i < t_onbeam->GetEntries(); i++){

    t_onbeam->GetEntry(i);

    h_nTracks_onbeam->Fill(onbeam_vars.nSelectedTracks);
    h_nShowers_onbeam->Fill(onbeam_vars.nSelectedShowers);
   // h_trackLength_onbeam->Fill(onbeam_vars.track_length);

  }
 
  // loop offbeam
  for (int i = 0; i < t_offbeam->GetEntries(); i++){

    t_offbeam->GetEntry(i);

    h_nTracks_offbeam->Fill(offbeam_vars.nSelectedTracks);
    h_nShowers_offbeam->Fill(offbeam_vars.nSelectedShowers);
   // h_trackLength_offbeam->Fill(offbeam_vars.track_length);

  }

  h_nTracks_offbeam->Sumw2();
  h_nShowers_offbeam->Sumw2();
  h_trackLength_offbeam->Sumw2();
  h_nTracks_offbeam->Scale(offbeamscaling);
  h_nShowers_offbeam->Scale(offbeamscaling);
  h_trackLength_offbeam->Scale(offbeamscaling);

  THStack *hs1 = new THStack("hs1", "hs1");
  h_nTracks_offbeam->SetFillColor(kBlack);
  h_nTracks_offbeam->SetFillStyle(3345);
  hs1->Add(h_nTracks_offbeam);
  h_nTracks_sim->SetFillColor(kGreen+1);
  hs1->Add(h_nTracks_sim);
  hs1->Draw("hist");

  h_nTracks_onbeam->SetMarkerStyle(20);
  h_nTracks_onbeam->Draw("samepE1");

  std::cout << h_nTracks_onbeam->Integral() << std::endl;

  return 0;
}
