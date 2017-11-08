void printPlotToFile(TFile* inputFile, TString plotName){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

  TEfficiency* eCC = (TEfficiency*)inputFile->Get("numusel/"+plotName+"Eff");
  TEfficiency* eCC0Pi = (TEfficiency*)inputFile->Get("numusel/"+plotName+"CC0PiEff");

  TGraphAsymmErrors* hCC = eCC->CreateGraph();
  TGraphAsymmErrors* hCC0Pi = eCC0Pi->CreateGraph();

  hCC->SetLineColor(kBlack);
  hCC0Pi->SetLineColor(kRed); 

  hCC->GetYaxis()->SetRangeUser(0,1);

  hCC->Draw("ap");
  hCC0Pi->Draw("psame");

  c1->SaveAs(plotName+"Eff.eps");

}

void makeEfficiencyPlots(){

  std::vector<TString> plotNames= {"mcNuEnergy", "mcLeptonMom", "mcLeptonTheta", "mcLeptonPhi"};

  for (int i = 0; i < plotNames.size(); i++){

    printPlotToFile(_file0, plotNames.at(i));
  
  }

}
