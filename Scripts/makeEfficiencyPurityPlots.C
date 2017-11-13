void printPlotToFile(TFile* inputFile, TString plotName, TString EffOrP){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

  TEfficiency* eCC = (TEfficiency*)inputFile->Get("numusel/"+plotName+EffOrP);
  TEfficiency* eCC0Pi = (TEfficiency*)inputFile->Get("numusel/"+plotName+"CC0Pi"+EffOrP);

  TGraphAsymmErrors* hCC = eCC->CreateGraph();
  TGraphAsymmErrors* hCC0Pi = eCC0Pi->CreateGraph();

  hCC->SetLineColor(kBlack);
  hCC0Pi->SetLineColor(kRed); 

  hCC->GetYaxis()->SetRangeUser(0,1);

  hCC->Draw("ap");
  hCC0Pi->Draw("psame");

  TLegend *leg = new TLegend(0.5, 0.15, 0.85, 0.3);
  leg->AddEntry(hCC, "CC");
  leg->AddEntry(hCC0Pi, "CC0Pi");
  leg->Draw("same");

  c1->SaveAs(plotName+EffOrP+".eps");

}

void makeEfficiencyPurityPlots(){

  std::vector<TString> plotNames= {"mcNuEnergy", "mcLeptonMom", "mcLeptonTheta", "mcLeptonPhi"};

  for (int i = 0; i < plotNames.size(); i++){

    printPlotToFile(_file0, plotNames.at(i), "Eff");
    printPlotToFile(_file0, plotNames.at(i), "Pur");
 
  }

}
