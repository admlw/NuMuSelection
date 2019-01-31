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

#include "DoEnergyRecoStudies.h"

std::pair<float,float> leastSquaresFitForMandC(TGraph* input_graph, TF1 *lin, float min_fit, float max_fit){

  float min_r = 1000000000;
  float min_m = 0;
  float min_c = 0;

  Double_t* xs = input_graph->GetX();
  Double_t* ys = input_graph->GetY();

  for (float m = 0.; m < 2.; m+=0.0001){
    for (float c = -0.1; c < 0.1; c+=0.001){

      lin->SetParameter(0, m);
      lin->SetParameter(1, c);

      float r = 0;

      for (int j = 0; j < input_graph->GetN(); j++){

        if (xs[j] < min_fit || xs[j] > max_fit) continue;

        r += std::pow(ys[j]-lin->Eval(xs[j]),2);

      }

      if (r < min_r){

        min_r = r;
        min_c = c;
        min_m = m;

      }
    }

  }

  std::pair<float, float> return_pair;
  return_pair.first = min_m;
  return_pair.second = min_c;

  return return_pair;

}

void leastSquaresFitToMedian(TH2D* input_histogram, float min_fit, float max_fit){

  std::vector<float> xs_v;
  std::vector<float> ys_v;

  for (int i = 0; i < input_histogram->GetNbinsX(); i++){

    float mpv = 0;
    int binno = 0;

    for (int j = 0; j < input_histogram->GetNbinsY(); j++){

      if (input_histogram->GetBinContent(i,j) > input_histogram->GetBinContent(i,binno)){
        mpv = input_histogram->GetYaxis()->GetBinCenter(j);
        binno = j;
      }

    }

    float bin_center = input_histogram->GetXaxis()->GetBinCenter(i);

    if ( bin_center > min_fit && bin_center < max_fit){
        xs_v.push_back(input_histogram->GetXaxis()->GetBinCenter(i));
        ys_v.push_back(mpv);
    }

  }

  float *xs = &xs_v[0];
  float *ys = &ys_v[0];

  TGraph* mpvs = new TGraph((int)xs_v.size(), xs, ys);
  mpvs->SetMarkerStyle(20);
  mpvs->SetMarkerColor(kBlack);
  mpvs->SetMarkerSize(0.6);
  mpvs->SetLineWidth(2);
  mpvs->Draw("psame");

  TF1* pol1 = new TF1("pol1", "[0]*x+[1]", 0, 10);
  std::pair<float,float> mandc = leastSquaresFitForMandC(mpvs, pol1, min_fit, max_fit);

  pol1->SetParameter(0,mandc.first);
  pol1->SetParameter(1,mandc.second);
  pol1->Draw("same");

  TString fitRange = Form("Fit Range: %0.2f - %0.2f GeV", min_fit, max_fit);
  TString fitResult = Form("Fit Result: E_{k}^{Reco} = %0.3f #times E_{k}^{True} + %0.3f", pol1->GetParameter(0), pol1->GetParameter(1));

  TPaveText *pt = new TPaveText(0.11, 0.75, 0.7, 0.85, "NDC");
  pt->AddText(fitRange);
  pt->AddText(fitResult);
  pt->SetTextColor(kWhite);
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->Draw("same");

  std::cout << "PRINTING SUMMARY OF ENERGY RECO STUDIES" << std::endl;
  std::cout << " protons by range:  " << std::endl;
  std::cout << "    M: " << pol1->GetParameter(0) << std::endl; 
  std::cout << "    C: " << pol1->GetParameter(1) << std::endl;


}

void doFit(std::vector<TH1D*> input_histograms, int pdg){

  float resolutions[input_histograms.size()];
  float resolutions_x[input_histograms.size()];
  float err_x[input_histograms.size()];
  float err_y[input_histograms.size()];

  for (int i = 0; i < input_histograms.size(); i++){

    TCanvas *c_res = new TCanvas();
    c_res->cd();
    input_histograms.at(i)->Draw();

    float peak_val = input_histograms.at(i)->GetMaximum();
    float peak_bin = input_histograms.at(i)->GetMaximumBin();
    float peak_bin_val = input_histograms.at(i)->GetBinCenter(input_histograms.at(i)->GetMaximumBin());
    float low_val = 0;
    float high_val = 0;
    for (int j = peak_bin; j > 0; j--){

      float test_val = input_histograms.at(i)->GetBinContent(j);

      if (test_val/peak_val < 0.3){
        low_val = input_histograms.at(i)->GetBinCenter(j);
        break;
      }
    }
    for (int j = peak_bin; j < 2*peak_bin; j++){

      float test_val = input_histograms.at(i)->GetBinContent(j);

      if (test_val/peak_val < 0.3){
        high_val = input_histograms.at(i)->GetBinCenter(j);
        break;
      }
    }

    TF1* fittingFunction = new TF1("fittingFunction", "gaus(0)*exp(3)", -1, 1);
    fittingFunction->SetParameters(1.0, peak_bin_val, 1.0);

    input_histograms.at(i)->Fit(fittingFunction, "I", "", low_val, high_val);


    float resolution = std::abs(fittingFunction->GetParameter(2));
    float resolution_err = std::abs(fittingFunction->GetParError(2));
    if (pdg == 2212){
      resolutions_x[i] = proton_binVals.at(i).at(0) + (proton_binVals.at(i).at(1) - proton_binVals.at(i).at(0))/2.0;
      err_x[i] = (proton_binVals.at(i).at(1)-proton_binVals.at(i).at(0))/2.0;
    }
    if (pdg == 13){
      resolutions_x[i] = muon_binVals.at(i).at(0) + (muon_binVals.at(i).at(1) - muon_binVals.at(i).at(0))/2.0;
      err_x[i] = (muon_binVals.at(i).at(1)-muon_binVals.at(i).at(0))/2.0;
    }

    resolutions[i] = (resolution*100);
    err_y[i] = resolution_err*100;

    c_res->SaveAs((std::string("plots/")+std::string(input_histograms.at(i)->GetName())+".pdf").c_str());
    c_res->SaveAs((std::string("plots/")+std::string(input_histograms.at(i)->GetName())+".png").c_str());

  }

  TCanvas *c_gr = new TCanvas();
  TGraphErrors *gr = new TGraphErrors(input_histograms.size(), resolutions_x, resolutions, err_x, err_y);
  
  if (pdg == 2212)
    gr->SetTitle(";True Proton Energy;Resolution (%)");
  else if (pdg == 13)
    gr->SetTitle(";True Muon Energy;Resolution (%)");

  gr->SetMarkerStyle(20);
  gr->GetYaxis()->SetRangeUser(0,100);
  gr->Draw("ap");

  c_gr->SaveAs((std::string("plots/")+std::string(input_histograms.at(0)->GetName())+"_resolutions.pdf").c_str());
  c_gr->SaveAs((std::string("plots/")+std::string(input_histograms.at(0)->GetName())+"_resolutions.png").c_str());

};


int main(){

  // pull out TTrees from provided TFiles
  // input trees are set in Configuration.h
  TFile* inputfile = new TFile("selectedEvents_sim.root", "READ");
  TTree* t_simulation = (TTree*)inputfile->Get("simulation");

  // initialise variables and trees
  var_list simulation_vars_tmp;

  var_list* simulation_vars = &simulation_vars_tmp;

  _treeHandler.SetTreeVars(t_simulation, &simulation_vars_tmp, true);

  //------------------------------------
  // loop simulation
  //------------------------------------

  TH2D* range_mom_mcs_mom_contained = new TH2D("range_mom_mcs_mom_contained", ";Momentum from Range (GeV); Momentum from MCS (GeV)", 50, 0, 2, 50, 0, 2);
  TH2D* true_energy_range_energy_proton = new TH2D("true_energy_range_energy_proton", "Protons;E_{k}^{True} (GeV); E_{k}^{Reco, range}", 100, 0.0, 0.6, 100, 0.0, 0.6);
  TH2D* true_energy_range_energy_muon = new TH2D("true_energy_range_energy_muon", "Muons;E_{k}^{True};E_{k}^{Reco, range}", 50, 0.0, 2.0, 50, 0.0, 2.0);
  TH2D* true_energy_range_energy_muon_contained = new TH2D("true_energy_range_energy_muon_contained", "Muons, Contained;E_{k}^{True};E_{k}^{Reco, range}", 50, 0.0, 2.0, 50, 0.0, 2.0);
  TH2D* true_energy_range_energy_muon_uncontained = new TH2D("true_energy_range_energy_muon_uncontained", "Muons, Uncontained;E_{k}^{True};E_{k}^{Reco, range}", 50, 0.0, 2.0, 50, 0.0, 2.0);
  TH2D* true_energy_mcs_energy_muon = new TH2D("true_energy_range_energy_muon", "Muons;E_{k}^{True};E_{k}^{Reco, MCS}", 50, 0.0, 2.0, 50, 0.0, 2.0);
  TH2D* true_energy_mcs_energy_muon_contained = new TH2D("true_energy_range_energy_muon_contained", "Muons, Contained;E_{k}^{True};E_{k}^{Reco, MCS}", 50, 0.0, 2.0, 50, 0.0, 2.0);
  TH2D* true_energy_mcs_energy_muon_uncontained = new TH2D("true_energy_range_energy_muon_uncontained", "Muons, Uncontained;E_{k}^{True};E_{k}^{Reco, MCS}", 50, 0.0, 2.0, 50, 0.0, 2.0);

  std::vector<TH1D*> proton_resolution_histograms;
  std::vector<TH1D*> muon_range_resolution_histograms;
  std::vector<TH1D*> muon_range_contained_resolution_histograms;
  std::vector<TH1D*> muon_range_uncontained_resolution_histograms;
  std::vector<TH1D*> muon_mcs_resolution_histograms;
  std::vector<TH1D*> muon_mcs_contained_resolution_histograms;
  std::vector<TH1D*> muon_mcs_uncontained_resolution_histograms;

  for (int i = 0; i < proton_binVals.size(); i++){

    TString binRange = Form("%0.2fto%0.2f", proton_binVals.at(i).at(0), proton_binVals.at(i).at(1));
    proton_resolution_histograms.push_back(new TH1D(std::string("proton_resolution_histograms_"+binRange).c_str(), ";(E_k^{range}-E_{k}^{true})/E_{k}^{true};N candidates", 50, -0.5, 0.5));

  }

  for (int i = 0; i < muon_binVals.size(); i++){

    TString binRange = Form("%0.2fto%0.2f", muon_binVals.at(i).at(0), muon_binVals.at(i).at(1));
    muon_range_resolution_histograms.push_back(new TH1D(std::string("muon_range_resolution_histograms_"+binRange).c_str(), ";;", 50, -0.75, 0.755));
    muon_range_contained_resolution_histograms.push_back(new TH1D(std::string("muon_range_contained_resolution_histograms_"+binRange).c_str(), ";(E_k^{range}-E_{k}^{true})/E_{k}^{true};N candidates", 50, -0.75, 0.75));
    muon_range_uncontained_resolution_histograms.push_back(new TH1D(std::string("muon_range_uncontained_resolution_histograms_"+binRange).c_str(), ";;", 50, -0.75, 0.75));
    muon_mcs_resolution_histograms.push_back(new TH1D(std::string("muon_mcs_resolution_histograms_"+binRange).c_str(), ";;", 50, -0.75, 0.75));
    muon_mcs_contained_resolution_histograms.push_back(new TH1D(std::string("muon_mcs_contained_resolution_histograms_"+binRange).c_str(), ";;", 50, -0.75, 0.75));
    muon_mcs_uncontained_resolution_histograms.push_back(new TH1D(std::string("muon_mcs_uncontained_resolution_histograms_"+binRange).c_str(), ";(E_k^{MCS}-E_{k}^{true})/E_{k}^{true};N candidates", 50, -0.75, 0.75));

  }

  std::cout << "[DES] Beginning simulation loop..." << std::endl;
  for (int i = 0; i < t_simulation->GetEntries(); i++){

    t_simulation->GetEntry(i);

    for (int j = 0; j < simulation_vars->true_match_starte->size(); j++){


      //
      // protons
      //
      if (simulation_vars->true_match_pdg->at(j) == 2212){


        float true_energy = simulation_vars->true_match_starte->at(j) - 0.938;

        float range_energy = simulation_vars->track_range_energy_passumption->at(j);
        float resolution = (range_energy - true_energy)/true_energy;

        for (int k = 0 ; k < proton_binVals.size(); k++){

          if (true_energy >= proton_binVals.at(k).at(0) && true_energy < proton_binVals.at(k).at(1))
            proton_resolution_histograms.at(k)->Fill(resolution);

        }

        true_energy_range_energy_proton->Fill(true_energy, range_energy);

      }

      //
      // muons
      //
      if (simulation_vars->true_match_pdg->at(j) == 13){
        float true_energy = simulation_vars->true_match_starte->at(j) - 0.105;

        float range_mom = simulation_vars->track_range_mom_muassumption->at(j);
        float range_energy = simulation_vars->track_range_energy_muassumption->at(j);


        // get mcs stuff
        float mcs_mom = 0;
        float mcs_energy = 0;

        if (simulation_vars->track_mcs_muassmp_fwd_loglikelihood->at(j) < simulation_vars->track_mcs_muassmp_bwd_loglikelihood->at(j)){
          mcs_mom = simulation_vars->track_mcs_muassmp_fwd->at(j);
          mcs_energy = simulation_vars->track_mcs_muassmp_energy_fwd->at(j);
        }
        else {
          mcs_mom = simulation_vars->track_mcs_muassmp_fwd->at(j);
          mcs_energy = simulation_vars->track_mcs_muassmp_energy_bwd->at(j);
        }

        true_energy_range_energy_muon->Fill(true_energy, range_energy);
        true_energy_mcs_energy_muon->Fill(true_energy, mcs_energy);

        float range_resolution = (range_energy - true_energy)/true_energy;
        float mcs_resolution = (mcs_energy - true_energy)/true_energy;

        for (int k = 0 ; k < muon_binVals.size(); k++){

          if (true_energy >= muon_binVals.at(k).at(0) && true_energy < muon_binVals.at(k).at(1)){
            muon_range_resolution_histograms.at(k)->Fill(range_resolution);
            muon_mcs_resolution_histograms.at(k)->Fill(mcs_resolution);
          }

        }

        // contained muons
        if (simulation_vars->track_isContained->at(j)){
          true_energy_range_energy_muon_contained->Fill(true_energy, range_energy);
          true_energy_mcs_energy_muon_contained->Fill(true_energy, mcs_energy);
          range_mom_mcs_mom_contained->Fill(range_mom, mcs_mom);

          for (int k = 0 ; k < muon_binVals.size(); k++){

            if (true_energy >= muon_binVals.at(k).at(0) && true_energy < muon_binVals.at(k).at(1)){
              muon_range_contained_resolution_histograms.at(k)->Fill(range_resolution);
              muon_mcs_contained_resolution_histograms.at(k)->Fill(mcs_resolution);
            }

          }



        }
        // uncontained muons
        else{
          true_energy_range_energy_muon_uncontained->Fill(true_energy, range_energy);
          true_energy_mcs_energy_muon_uncontained->Fill(true_energy, mcs_energy);

          for (int k = 0 ; k < muon_binVals.size(); k++){

            if (true_energy >= muon_binVals.at(k).at(0) && true_energy < muon_binVals.at(k).at(1)){
              muon_range_uncontained_resolution_histograms.at(k)->Fill(range_resolution);
              muon_mcs_uncontained_resolution_histograms.at(k)->Fill(mcs_resolution);
            }

          }


        }

      }

    }

  }


  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);

  doFit(proton_resolution_histograms, 2212);
  //doFit(muon_range_resolution_histograms, 13);
  doFit(muon_range_contained_resolution_histograms, 13);
  //doFit(muon_range_uncontained_resolution_histograms, 13);
  //doFit(muon_mcs_resolution_histograms, 13);
  //doFit(muon_mcs_contained_resolution_histograms, 13);
  doFit(muon_mcs_uncontained_resolution_histograms, 13);

  TCanvas *c1 = new TCanvas();
  true_energy_range_energy_proton->GetZaxis()->SetRangeUser(-0.0001, true_energy_range_energy_proton->GetMaximum());
  true_energy_range_energy_proton->SetContour(1000);
  true_energy_range_energy_proton->Draw("colz");
  leastSquaresFitToMedian(true_energy_range_energy_proton, 0.05, 0.30);

  c1->SaveAs("plots/true_energy_range_energy_proton.pdf");
  c1->SaveAs("plots/true_energy_range_energy_proton.png");

  TCanvas *c2 = new TCanvas();
  true_energy_range_energy_muon->GetZaxis()->SetRangeUser(-0.0001, true_energy_range_energy_muon->GetMaximum());
  true_energy_range_energy_muon->SetContour(1000);
  true_energy_range_energy_muon->Draw("colz");
  leastSquaresFitToMedian(true_energy_range_energy_muon, 0.05, 0.6);

  c2->SaveAs("plots/true_energy_range_energy_muon.pdf");
  c2->SaveAs("plots/true_energy_range_energy_muon.png");

  TCanvas *c3 = new TCanvas();
  true_energy_range_energy_muon_contained->GetZaxis()->SetRangeUser(-0.0001, true_energy_range_energy_muon_contained->GetMaximum());
  true_energy_range_energy_muon_contained->SetContour(1000);
  true_energy_range_energy_muon_contained->Draw("colz");
  leastSquaresFitToMedian(true_energy_range_energy_muon_contained, 0.05, 0.6);

  c3->SaveAs("plots/true_energy_range_energy_muon_contained.pdf");
  c3->SaveAs("plots/true_energy_range_energy_muon_contained.png");

  TCanvas *c5 = new TCanvas();
  true_energy_mcs_energy_muon->GetZaxis()->SetRangeUser(-0.0001, true_energy_mcs_energy_muon->GetMaximum());
  true_energy_mcs_energy_muon->SetContour(1000);
  true_energy_mcs_energy_muon->Draw("colz");
  leastSquaresFitToMedian(true_energy_mcs_energy_muon, 0.05, 0.6);

  c5->SaveAs("plots/true_energy_mcs_energy_muon.pdf");
  c5->SaveAs("plots/true_energy_mcs_energy_muon.png");

  TCanvas *c7 = new TCanvas();
  true_energy_mcs_energy_muon_uncontained->GetZaxis()->SetRangeUser(-0.0001, true_energy_mcs_energy_muon_uncontained->GetMaximum());
  true_energy_mcs_energy_muon_uncontained->SetContour(1000);
  true_energy_mcs_energy_muon_uncontained->Draw("colz");
  leastSquaresFitToMedian(true_energy_mcs_energy_muon_uncontained, 0.1, 0.9);

  c7->SaveAs("plots/true_energy_mcs_energy_muon_uncontained.pdf");
  c7->SaveAs("plots/true_energy_mcs_energy_muon_uncontained.png");

  TCanvas *c8 = new TCanvas();
  range_mom_mcs_mom_contained->GetZaxis()->SetRangeUser(-0.0001, range_mom_mcs_mom_contained->GetMaximum());
  range_mom_mcs_mom_contained->SetContour(1000);
  range_mom_mcs_mom_contained->Draw("colz");

  c8->SaveAs("plots/range_mom_mcs_mom_contained.pdf");
  c8->SaveAs("plots/range_mom_mcs_mom_contained.png");

  return 0;

}
