void makeProtonMultiplicityPlot(){

  int chosenpdg = 2212;
  string species = "Proton";

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  TPad *topPad = new TPad("topPad", "", 0.005, 0.3, 0.995, 0.995);
  TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
  topPad->SetTopMargin(0.1);
  topPad->SetBottomMargin(0.005);
  bottomPad->SetTopMargin(0.005);
  bottomPad->SetBottomMargin(0.25);
  topPad->Draw();
  bottomPad->Draw();

  TTree* t_nom = (TTree*)_file0->Get("numuselection/analysis_tree");
  TTree* t_dic = (TTree*)_file1->Get("numuselection/analysis_tree");

  TH1D* h_protonmult_nom = new TH1D("h_protonmult_nom", ";Reconstructed Proton Multiplicity;", 8, 0, 8);
  TH1D* h_protonmult_dic = new TH1D("h_protonmult_dic", ";Reconstructed Proton Multiplicity;", 8, 0, 8);

  // set branch address - nominal
  std::vector<int>* true_match_pdg_nom = nullptr;
  t_nom->SetBranchAddress("true_match_pdg", &true_match_pdg_nom);

  // set branch address - dic
  std::vector<int>* true_match_pdg_dic = nullptr;
  t_dic->SetBranchAddress("true_match_pdg", &true_match_pdg_dic);

  // nominal loop
  for (size_t i = 0; i < t_nom->GetEntries(); i++){

    t_nom->GetEntry(i);

    int nprotons = 0;

    for (int j = 0; j < true_match_pdg_nom->size(); j++){
      if (std::abs(true_match_pdg_nom->at(j)) == chosenpdg)
        nprotons++;

    }

    if (true_match_pdg_nom->at(0) != -999)
      h_protonmult_nom->Fill(nprotons);

  }

  // dic loop
  for (size_t i = 0; i < t_dic->GetEntries(); i++){

    t_dic->GetEntry(i);

    int nprotons = 0;

    for (int j = 0; j < true_match_pdg_dic->size(); j++){
      if (std::abs(true_match_pdg_dic->at(j)) == chosenpdg)
        nprotons++;

    }

    if (true_match_pdg_dic->at(0) != -999)
      h_protonmult_dic->Fill(nprotons);

  }

  topPad->cd();

  h_protonmult_nom->SetLineColor(kAzure+1);
  h_protonmult_nom->SetLineWidth(2);
  h_protonmult_nom->Sumw2();
  h_protonmult_nom->Scale(1.0);
  h_protonmult_nom->SetLabelOffset(10000);
  h_protonmult_nom->GetYaxis()->SetRangeUser(0.1, std::max(h_protonmult_nom->GetMaximum(), h_protonmult_dic->GetMaximum()) *1.25);
  h_protonmult_nom->Draw("hist");

  h_protonmult_dic->SetLineColor(kGreen+1);
  h_protonmult_dic->SetLineWidth(2);
  h_protonmult_dic->Sumw2();
  h_protonmult_dic->Scale(0.99);
  h_protonmult_dic->Draw("hist same");

  TLegend *leg = new TLegend(0.5, 0.75, 0.8, 0.85);
  leg->AddEntry(h_protonmult_nom, "Nominal");
  leg->AddEntry(h_protonmult_dic, "DIC");
  leg->Draw("same");

  bottomPad->cd();
  TH1D* protonratio = (TH1D*)h_protonmult_dic->Clone("protonratio");
  protonratio->Divide(h_protonmult_nom);
  protonratio->SetLineColor(kBlack);
  protonratio->SetTitle(std::string(";"+species+" Multiplicity;DIC/NOM").c_str());
  protonratio->GetYaxis()->SetTitleSize(0.1);
  protonratio->GetXaxis()->SetTitleSize(0.1);
  protonratio->GetYaxis()->SetLabelSize(0.1);
  protonratio->GetXaxis()->SetLabelSize(0.1);
  protonratio->GetYaxis()->SetTitleOffset(0.4);
  protonratio->GetYaxis()->SetRangeUser(0.5, 1.5);
  protonratio->Draw();


  int nomtotal = h_protonmult_nom->Integral(0,8);
  int nomnonzero = h_protonmult_nom->Integral(2,8);

  int dictotal = h_protonmult_dic->Integral(0,8);
  int dicnonzero = h_protonmult_dic->Integral(2,8);

  std::cout << "nominal ratio protons/total: " << (double)nomnonzero/nomtotal << std::endl;
  std::cout << "dic ratio protons/total: " << (double)dicnonzero/dictotal << std::endl;
  std::cout << "dic nonzer/nom nonzero: " << (double)dicnonzero/nomnonzero << std::endl;

  c1->SaveAs(std::string(species+"_multiplicity.pdf").c_str());

}
