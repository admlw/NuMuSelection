void makeShowerFractionAsFnOfThetaXZPlots(){

  int chosenPDG = 13;
  double massPDG = 0.105;
  std::string namePDG = "Muon";
  int nbins = 20;
  float binlow = -3.15;
  float binhigh = 3.15;

  std::vector<int>* pfp_pdgCode_nom = nullptr;
  std::vector<int>* true_match_pdg_nom = nullptr;
  std::vector<float>* true_match_startx_nom = nullptr;
  std::vector<float>* true_match_endx_nom = nullptr;
  std::vector<float>* true_match_startz_nom = nullptr;
  std::vector<float>* true_match_endz_nom = nullptr;

  std::vector<int>* pfp_pdgCode_dic = nullptr;
  std::vector<int>* true_match_pdg_dic = nullptr;
  std::vector<float>* true_match_startx_dic = nullptr;
  std::vector<float>* true_match_endx_dic = nullptr;
  std::vector<float>* true_match_startz_dic = nullptr;
  std::vector<float>* true_match_endz_dic = nullptr;

  TEfficiency* trackFrac_nom = new TEfficiency("trackFrac_nom", (namePDG+std::string(";True #theta_{XZ};")).c_str(), nbins, binlow, binhigh);
  TEfficiency* trackFrac_dic = new TEfficiency("trackFrac_dic", (namePDG+std::string(";True #theta_{XZ};")).c_str(), nbins, binlow, binhigh);

  TTree* t0 = (TTree*)_file0->Get("numuselection/analysis_tree");
  TTree* t1 = (TTree*)_file1->Get("numuselection/analysis_tree");

  t0->SetBranchAddress("pfp_pdgCode", &pfp_pdgCode_nom);
  t0->SetBranchAddress("true_match_pdg", &true_match_pdg_nom);
  t0->SetBranchAddress("true_match_startx", &true_match_startx_nom);
  t0->SetBranchAddress("true_match_endx", &true_match_endx_nom);
  t0->SetBranchAddress("true_match_startz", &true_match_startz_nom);
  t0->SetBranchAddress("true_match_endz", &true_match_endz_nom);


  t1->SetBranchAddress("pfp_pdgCode", &pfp_pdgCode_dic);
  t1->SetBranchAddress("true_match_pdg", &true_match_pdg_dic);
  t1->SetBranchAddress("true_match_startx", &true_match_startx_dic);
  t1->SetBranchAddress("true_match_endx", &true_match_endx_dic);
  t1->SetBranchAddress("true_match_startz", &true_match_startz_dic);
  t1->SetBranchAddress("true_match_endz", &true_match_endz_dic);

  for (int i_en = 0; i_en < t0->GetEntries(); i_en++){

    t0->GetEntry(i_en);

    for (int i = 0; i < true_match_pdg_nom->size(); i++){

        if (std::abs(true_match_pdg_nom->at(i)) == chosenPDG){
         
            if (pfp_pdgCode_nom->at(i) == 13)
                trackFrac_nom->Fill(true, std::atan2(true_match_endx_nom->at(i)-true_match_startx_nom->at(i),true_match_endz_nom->at(i)-true_match_startz_nom->at(i)));
            else trackFrac_nom->Fill(false, std::atan2(true_match_endx_nom->at(i)-true_match_startx_nom->at(i),true_match_endz_nom->at(i)-true_match_startz_nom->at(i)));
        }
    }
  }

  for (int i_en = 0; i_en < t1->GetEntries(); i_en++){

    t1->GetEntry(i_en);

    for (int i = 0; i < true_match_pdg_dic->size(); i++){

        if (std::abs(true_match_pdg_dic->at(i)) == chosenPDG){

            if (pfp_pdgCode_dic->at(i) == 13)
                trackFrac_dic->Fill(true, std::atan2(true_match_endx_dic->at(i)-true_match_startx_dic->at(i),true_match_endz_dic->at(i)-true_match_startz_dic->at(i)) );
            else trackFrac_dic->Fill(false, std::atan2(true_match_endx_dic->at(i)-true_match_startx_dic->at(i),true_match_endz_dic->at(i)-true_match_startz_dic->at(i)) );
        }
    }
  }

  TH1* nevents_nom = (TH1*)trackFrac_nom->GetCopyTotalHisto();
 
  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  TPad *topPad = new TPad("topPad", "", 0.005, 0.3, 0.995, 0.995);
  TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
  topPad->SetBottomMargin(0.02);
  bottomPad->SetTopMargin(0.0);
  bottomPad->SetBottomMargin(0.25);
  bottomPad->SetGridy();
  topPad->Draw();
  bottomPad->Draw();
  topPad->cd();

  nevents_nom->SetFillColor(kGray);
  nevents_nom->SetLineWidth(0);
  nevents_nom->Scale(1./nevents_nom->Integral());
  nevents_nom->GetYaxis()->SetRangeUser(0,1);
  nevents_nom->GetXaxis()->SetLabelOffset(1000);
  nevents_nom->SetTitle((namePDG+";;Fraction Reconstructed as Tracks").c_str());
  nevents_nom->Draw();
  trackFrac_nom->SetMarkerColor(kAzure+1);
  trackFrac_nom->SetLineColor(kAzure+1);
  trackFrac_nom->Draw("same");
  trackFrac_dic->SetMarkerColor(kGreen+1);
  trackFrac_dic->SetLineColor(kGreen+1);
  trackFrac_dic->Draw("same");

  gPad->RedrawAxis();

  bottomPad->cd();
  TH1* nevents_total_nom  = (TH1*)trackFrac_nom->GetCopyTotalHisto();
  TH1* nevents_passed_nom = (TH1*)trackFrac_nom->GetCopyPassedHisto();

  nevents_total_nom->Sumw2();
  nevents_passed_nom->Sumw2();

  nevents_passed_nom->Divide(nevents_total_nom);

  TH1* nevents_total_dic  = (TH1*)trackFrac_dic->GetCopyTotalHisto();
  TH1* nevents_passed_dic = (TH1*)trackFrac_dic->GetCopyPassedHisto();

  nevents_total_dic->Sumw2();
  nevents_passed_dic->Sumw2();

  nevents_passed_dic->Divide(nevents_total_dic);

  nevents_passed_dic->Divide(nevents_passed_nom);
  
  nevents_passed_dic->GetYaxis()->SetRangeUser(0.8, 1.2);
  nevents_passed_dic->SetTitle(";True E_{k} (GeV);DIC/Nominal");
  nevents_passed_dic->GetXaxis()->SetTitleSize(0.12);
  nevents_passed_dic->GetYaxis()->SetTitleSize(0.12);
  nevents_passed_dic->GetXaxis()->SetLabelSize(0.12);
  nevents_passed_dic->GetYaxis()->SetLabelSize(0.12);
  nevents_passed_dic->GetYaxis()->SetTitleOffset(0.38);
  nevents_passed_dic->GetYaxis()->SetNdivisions(305);

  nevents_passed_dic->Draw("sameE1");

  c1->SaveAs((namePDG+std::string("trackFracThetaXZ.pdf")).c_str());

}
