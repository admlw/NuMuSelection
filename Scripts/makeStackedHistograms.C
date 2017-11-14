void printPlotToFile(TFile* inputFile, TString plotName){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

  THStack *hs = new THStack("hs", "");

  TH1D* hOther = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"Other");
  TH1D* hNC = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"NC");
  TH1D* hAntiNuE = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"AntiNuE");
  TH1D* hNuE = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"NuE");
  TH1D* hAntiNuMu = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"AntiNuMu");
  TH1D* hOOFV = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"OutOfFV");
  TH1D* hMixed = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"Mixed");
  TH1D* hCosmic = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"Cosmic");
  TH1D* hCC = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"TrueCC");
  TH1D* hCC0Pi = (TH1D*)inputFile->Get("numusel/selectedTrack"+plotName+"TrueCC0Pi");

  
  hCC0Pi->SetFillColor(TColor::GetColor( 215, 48, 39 ));
  hCC->SetFillColor(TColor::GetColor(165, 0, 38));
  hAntiNuMu->SetFillColor(TColor::GetColor(0, 104, 55));
  hNuE->SetFillColor(TColor::GetColor(166, 217, 106));//(102, 189, 99));
  hAntiNuE->SetFillColor(TColor::GetColor(166, 217, 106));
  hNC->SetFillColor(TColor::GetColor(113, 1, 98));
  hOOFV->SetFillColor(TColor::GetColor(78, 179, 211));
  hMixed->SetFillColor(TColor::GetColor(8, 104, 172));
  hCosmic->SetFillColor(TColor::GetColor(8, 64, 129));
  hOther->SetFillColor(TColor::GetColor(197, 197, 197));

  hCC->Add(hCC0Pi, -1);

  hCC0Pi->SetLineWidth(0);
  hCC->SetLineWidth(0);
  hAntiNuMu->SetLineWidth(0);
  hNuE->SetLineWidth(0);
  hAntiNuE->SetLineWidth(0);
  hNC->SetLineWidth(0);
  hOOFV->SetLineWidth(0);
  hMixed->SetLineWidth(0);
  hCosmic->SetLineWidth(0);
  hOther->SetLineWidth(0);


  hs->Add(hOther);
  hs->Add(hCosmic);
  hs->Add(hMixed);
  hs->Add(hOOFV);
  hs->Add(hNC);
  hs->Add(hAntiNuE);
  hs->Add(hNuE);
  hs->Add(hAntiNuMu);
  hs->Add(hCC);
  hs->Add(hCC0Pi);

  TH1D* hTotal = new TH1D("hTotal", "", 
      hCC0Pi->GetNbinsX(), 
      hCC0Pi->GetBinLowEdge(1), 
      hCC0Pi->GetBinLowEdge(hCC0Pi->GetNbinsX()+1));


  hTotal->SetFillStyle(3005);
  hTotal->Add(hOther);
  hTotal->Add(hCosmic);
  hTotal->Add(hMixed);
  hTotal->Add(hOOFV);
  hTotal->Add(hNC);
  hTotal->Add(hAntiNuE);
  hTotal->Add(hNuE);
  hTotal->Add(hAntiNuMu);
  hTotal->Add(hCC);
  hTotal->Add(hCC0Pi);
  hTotal->SetFillColor(kBlack);

  hTotal->Draw();

  hs->SetMaximum(hTotal->GetMaximum()*1.1);
  hs->Draw();

  hs->GetXaxis()->SetTitle("Candidate Muon Track Length (cm)");
  hs->GetYaxis()->SetTitleOffset(1.4); 
  hs->GetYaxis()->SetTitle("Number of Selected Events");
  
  hTotal->Draw("E2same");

  TLegend *leg = new TLegend(0.5, 0.3, 0.85, 0.85);
  leg->AddEntry(hCC0Pi, "#nu_{#mu} CC0Pi");
  leg->AddEntry(hCC, "Other #nu_{#mu} CC");
  leg->AddEntry(hAntiNuMu, "#bar{#nu_{#mu}}");
  leg->AddEntry(hNuE, "#nu_{e}, #bar{#nu_{e}}");
//  leg->AddEntry(hAntiNuE, "#bar{#nu_{e}}");
  leg->AddEntry(hNC, "NC");
  leg->AddEntry(hOOFV, "OOFV");
  leg->AddEntry(hMixed, "Mixed");
  leg->AddEntry(hCosmic, "Cosmic");
  leg->AddEntry(hOther, "Other");
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->Draw("same");

  c1->SaveAs(plotName+".eps");

}

void makeStackedHistograms(){

  std::vector<TString> plotNames= {"Length"};

  for (int i = 0; i < plotNames.size(); i++){

    printPlotToFile(_file0, plotNames.at(i));

  }

}
