#include "HistogramHandler.h"

namespace numusel{

  std::pair<double,int> HistogramHandler::CalcChi2(TH1D* E, TH1D* O){

    if (E->GetNbinsX() != O->GetNbinsX()) throw std::exception();

    double chi2 = 0;
    int nbins = 0;
    for (int i = 0; i < E->GetNbinsX()+1; i++){
      if (O->GetBinContent(i) == 0 || E->GetBinContent(i) == 0) continue;
      chi2 += std::pow(O->GetBinContent(i) - E->GetBinContent(i),2)/E->GetBinContent(i);
      nbins++;
    }

    std::pair<double, int> chi2ndof(chi2, nbins);

    return chi2ndof;

  }

  void HistogramHandler::Fill2DHistMC(hists_2d* h2d, std::vector<std::pair<double, double>> variable){

    for (int i = 0; i < variable.size(); i++){

      h2d->h_mc->Fill(variable.at(i).first, variable.at(i).second);

    }

  };

  void HistogramHandler::Fill2DHistOnbeam(hists_2d* h2d, std::vector<std::pair<double, double>> variable){

    for (int i = 0; i < variable.size(); i++){

      h2d->h_onbeam->Fill(variable.at(i).first, variable.at(i).second);

    }

  };

  void HistogramHandler::Fill2DHistOffbeam(hists_2d* h2d, std::vector<std::pair<double, double>> variable){

    for (int i = 0; i < variable.size(); i++){

      h2d->h_offbeam->Fill(variable.at(i).first, variable.at(i).second);

    }

  };


  void HistogramHandler::FillHistMC(hists_1d* h1d, std::vector<double> variable, std::bitset<8> eventCat){

    for (int i = 0; i < variable.size(); i++){
      // cosmic
      if (eventCat[7] == 1)
        h1d->h_mccosmic->Fill(variable.at(i));

      // mixed
      else if (eventCat[6] == 1)
        h1d->h_mcmixed->Fill(variable.at(i));

      // oofv
      else if (eventCat[5] == 1)
        h1d->h_mcoofv->Fill(variable.at(i));

      // nc
      else if (eventCat[4] == 1)
        h1d->h_mcnc->Fill(variable.at(i));

      // nuenuebar
      else if (eventCat[3] == 1)
        h1d->h_mcnuenuebar->Fill(variable.at(i));

      // numubar
      else if (eventCat[2] == 1)
        h1d->h_mcnumubar->Fill(variable.at(i));

      // numu-other
      else if (eventCat[1] == 1)
        h1d->h_mcnumuccother->Fill(variable.at(i));

      // numucc0pi
      else if (eventCat[0] == 1)
        h1d->h_mcnumucc0pinp->Fill(variable.at(i));

    }

  };

  void HistogramHandler::FillHistOnBeam(hists_1d* h1d, std::vector<double> variable){

    for (int i = 0; i < variable.size(); i++){
      h1d->h_onbeam->Fill(variable.at(i));
    }

  };

  void HistogramHandler::FillHistOffBeam(hists_1d* h1d, std::vector<double> variable){

    for (int i = 0; i < variable.size(); i++){
      h1d->h_offbeam->Fill(variable.at(i));
    }

  };

  void HistogramHandler::FillTrackHistMC(trackhists_1d* h1d, std::vector<double> variable, std::vector<double> pid, std::vector<double> isCutPassed){

    HistogramHandler _histoHandler;
    Configuration _config;

    if (_config.MakeTrackPlots == false) return;

    for (int i = 0; i < variable.size(); i++){

      if ((_config.UseTrackCut == true && isCutPassed.at(i) == 1) || (_config.UseTrackCut == false)){

        if (pid.at(i) == 13)
          h1d->h_muon->Fill(variable.at(i));
  
        else if (pid.at(i) == 2212)
          h1d->h_proton->Fill(variable.at(i));
  
        else if (pid.at(i) == 211)
          h1d->h_pion->Fill(variable.at(i));
  
        else if (pid.at(i) == 321)
          h1d->h_kaon->Fill(variable.at(i));
  
        else if (pid.at(i) == 11)
          h1d->h_electron->Fill(variable.at(i));
  
        else
          h1d->h_other->Fill(variable.at(i));
      }
    }

  };

  void HistogramHandler::FillTrackHistOnBeam(trackhists_1d* h1d, std::vector<double> variable, std::vector<double> isCutPassed){

    HistogramHandler _histoHandler;
    Configuration _config;

    if (_config.MakeTrackPlots == false) return;

    for (int i = 0; i < variable.size(); i++){
      
      if ((_config.UseTrackCut == true && isCutPassed.at(i) == 1) || (_config.UseTrackCut == false)){
         h1d->h_onbeam->Fill(variable.at(i));
      }
    }
  };

  void HistogramHandler::FillTrackHistOffBeam(trackhists_1d* h1d, std::vector<double> variable, std::vector<double> isCutPassed){

    HistogramHandler _histoHandler;
    Configuration _config;

    if (_config.MakeTrackPlots == false) return;

    for (int i = 0; i < variable.size(); i++){
      if ((_config.UseTrackCut == true && isCutPassed.at(i) == 1) || (_config.UseTrackCut == false)){
        h1d->h_offbeam->Fill(variable.at(i));
      }
    }

  };


  void HistogramHandler::InitialiseHistoVec(std::vector< std::vector<hists_1d*> > *plots_to_make, int n_plots){

    numusel::HistogramHandler _histohandler;

    for (int i_st = 0; i_st < _config.n_stages; i_st++){
      for (int i_pl = 0; i_pl < n_plots; i_pl++){

        plots_to_make->at(i_st).at(i_pl) = new hists_1d(
            std::string("h_"+_histohandler.histoNames.at(i_pl)+"_stage"+std::to_string(i_st)),
            _histohandler.histoLabels.at(i_pl),
            _histohandler.histoBins.at(i_pl).at(0),
            _histohandler.histoBins.at(i_pl).at(1),
            _histohandler.histoBins.at(i_pl).at(2)
            );
      }
    }
  };

  void HistogramHandler::InitialiseTrackHistoVec(std::vector< std::vector<trackhists_1d*> > *plots_to_make, int n_plots){

    numusel::HistogramHandler _histohandler;

    for (int i_st = 0; i_st < _config.n_stages; i_st++){
      for (int i_pl = 0; i_pl < n_plots; i_pl++){

        if (_histohandler.histoMakeTrackPlot.at(i_pl) == false) continue;

        plots_to_make->at(i_st).at(i_pl) = new trackhists_1d(
            std::string("h_"+_histohandler.histoNames.at(i_pl)+"tracks_stage"+std::to_string(i_st)),
            _histohandler.histoLabels.at(i_pl),
            _histohandler.histoBins.at(i_pl).at(0),
            _histohandler.histoBins.at(i_pl).at(1),
            _histohandler.histoBins.at(i_pl).at(2)
            );
      }
    }
  };

  void HistogramHandler::InitialiseHistoVec(std::vector< std::vector<eff_1d*> > *plots_to_make, int n_plots){

    numusel::HistogramHandler _histohandler;

    for (int i_st = 0; i_st < _config.n_stages; i_st++){
      for (int i_pl = 0; i_pl < n_plots; i_pl++){

        plots_to_make->at(i_st).at(i_pl) = new eff_1d(
            std::string("h_"+_histohandler.effNames.at(i_pl)+"_stage"+std::to_string(i_st)),
            _histohandler.effLabels.at(i_pl),
            _histohandler.effBins.at(i_pl).at(0),
            _histohandler.effBins.at(i_pl).at(1),
            _histohandler.effBins.at(i_pl).at(2)
            );
      }
    }
  };

  void HistogramHandler::InitialiseHistoVec(std::vector< std::vector<hists_2d*> > *plots_to_make, int n_plots){

    numusel::HistogramHandler _histohandler;

    for (int i_st = 0; i_st < _config.n_stages; i_st++){
      for (int i_pl = 0; i_pl < n_plots; i_pl++){

        plots_to_make->at(i_st).at(i_pl) = new hists_2d(
            std::string("h_"+_histohandler.histoNames_2D.at(i_pl)+"_stage"+std::to_string(i_st)),
            _histohandler.histoLabels_2D.at(i_pl),
            _histohandler.histoBins_2D.at(i_pl).at(0),
            _histohandler.histoBins_2D.at(i_pl).at(1),
            _histohandler.histoBins_2D.at(i_pl).at(2),
            _histohandler.histoBins_2D.at(i_pl).at(3),
            _histohandler.histoBins_2D.at(i_pl).at(4),
            _histohandler.histoBins_2D.at(i_pl).at(5)
            );
      }
    }
  };

  void HistogramHandler::StyleHistograms(hists_1d* hists){

    hists->h_mccosmic->SetFillColor(TColor::GetColor(8,64,129));
    hists->h_mcmixed->SetFillColor(TColor::GetColor(8,104,172));
    hists->h_mcoofv->SetFillColor(TColor::GetColor(78,179,211));
    hists->h_mcnc->SetFillColor(TColor::GetColor(113,1,98));
    hists->h_mcnuenuebar->SetFillColor(TColor::GetColor(166,217,106));
    hists->h_mcnumubar->SetFillColor(TColor::GetColor(0,104,55));
    hists->h_mcnumuccother->SetFillColor(TColor::GetColor(165,0,38));
    hists->h_mcnumucc0pinp->SetFillColor(TColor::GetColor(215,48,39));

    hists->h_mccosmic->SetMarkerColor(TColor::GetColor(8,64,129));
    hists->h_mcmixed->SetMarkerColor(TColor::GetColor(8,104,172));
    hists->h_mcoofv->SetMarkerColor(TColor::GetColor(78,179,211));
    hists->h_mcnc->SetMarkerColor(TColor::GetColor(113,1,98));
    hists->h_mcnuenuebar->SetMarkerColor(TColor::GetColor(166,217,106));
    hists->h_mcnumubar->SetMarkerColor(TColor::GetColor(0,104,55));
    hists->h_mcnumuccother->SetMarkerColor(TColor::GetColor(165,0,38));
    hists->h_mcnumucc0pinp->SetMarkerColor(TColor::GetColor(215,48,39));

    hists->h_mccosmic->SetLineWidth(0);
    hists->h_mcmixed->SetLineWidth(0);
    hists->h_mcoofv->SetLineWidth(0);
    hists->h_mcnc->SetLineWidth(0);
    hists->h_mcnuenuebar->SetLineWidth(0);
    hists->h_mcnumubar->SetLineWidth(0);
    hists->h_mcnumuccother->SetLineWidth(0);
    hists->h_mcnumucc0pinp->SetLineWidth(0);


    hists->h_offbeam->SetFillColor(kBlack);
    hists->h_offbeam->SetFillStyle(3345);

    hists->h_onbeam->SetLineColor(kBlack);
    hists->h_onbeam->SetMarkerStyle(20);
    hists->h_onbeam->SetMarkerSize(0.6);
  };

  void HistogramHandler::StyleTrackHistograms(trackhists_1d* hists){

    hists->h_muon->SetFillColor(TColor::GetColor(8, 64, 129));
    hists->h_proton->SetFillColor(TColor::GetColor(215, 48, 39));
    hists->h_pion->SetFillColor(TColor::GetColor(166, 217, 106));
    hists->h_kaon->SetFillColor(TColor::GetColor(133, 1, 98));
    hists->h_electron->SetFillColor(TColor::GetColor(251, 233, 69));
    hists->h_other->SetFillColor(TColor::GetColor(197, 197, 197));
    
    hists->h_muon->SetMarkerColor(TColor::GetColor(8, 64, 129));
    hists->h_proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
    hists->h_pion->SetMarkerColor(TColor::GetColor(166, 217, 106));
    hists->h_kaon->SetMarkerColor(TColor::GetColor(133, 1, 98));
    hists->h_electron->SetMarkerColor(TColor::GetColor(251, 233, 69));
    hists->h_other->SetMarkerColor(TColor::GetColor(197, 197, 197));

    hists->h_muon->SetLineWidth(0);
    hists->h_proton->SetLineWidth(0);
    hists->h_pion->SetLineWidth(0);
    hists->h_kaon->SetLineWidth(0);
    hists->h_electron->SetLineWidth(0);
    hists->h_other->SetLineWidth(0);

    hists->h_offbeam->SetFillColor(kBlack);
    hists->h_offbeam->SetFillStyle(3345);

    hists->h_onbeam->SetLineColor(kBlack);
    hists->h_onbeam->SetMarkerStyle(20);
    hists->h_onbeam->SetMarkerSize(0.6);
  };


  void HistogramHandler::ScaleHistograms(std::vector< std::vector<hists_1d*> > hists){

    for (int i_st = 0; i_st < hists.size(); i_st++){
      for (int i_pl = 0; i_pl < hists.at(i_st).size(); i_pl++){

        hists_1d* thisHistSet = hists.at(i_st).at(i_pl);

        thisHistSet->h_mccosmic->Sumw2();
        thisHistSet->h_mcmixed->Sumw2();
        thisHistSet->h_mcoofv->Sumw2();
        thisHistSet->h_mcnc->Sumw2();
        thisHistSet->h_mcnuenuebar->Sumw2();
        thisHistSet->h_mcnumubar->Sumw2();
        thisHistSet->h_mcnumuccother->Sumw2();
        thisHistSet->h_mcnumucc0pinp->Sumw2();
        thisHistSet->h_offbeam->Sumw2();

        thisHistSet->h_mccosmic->Scale(_config.simscaling);
        thisHistSet->h_mcmixed->Scale(_config.simscaling);
        thisHistSet->h_mcoofv->Scale(_config.simscaling);
        thisHistSet->h_mcnc->Scale(_config.simscaling);
        thisHistSet->h_mcnuenuebar->Scale(_config.simscaling);
        thisHistSet->h_mcnumubar->Scale(_config.simscaling);
        thisHistSet->h_mcnumuccother->Scale(_config.simscaling);
        thisHistSet->h_mcnumucc0pinp->Scale(_config.simscaling);
        thisHistSet->h_offbeam->Scale(_config.offbeamscaling);
      }
    }

  };

  void HistogramHandler::ScaleTrackHistograms(std::vector< std::vector<trackhists_1d*> > hists){

    HistogramHandler _histoHandler;
     
    for (int i_st = 0; i_st < hists.size(); i_st++){
      for (int i_pl = 0; i_pl < hists.at(i_st).size(); i_pl++){
      
        if (_histoHandler.histoMakeTrackPlot.at(i_pl) == false) continue;
       
        trackhists_1d* thisHistSet = hists.at(i_st).at(i_pl);

        thisHistSet->h_muon->Sumw2();
        thisHistSet->h_proton->Sumw2();
        thisHistSet->h_pion->Sumw2();
        thisHistSet->h_kaon->Sumw2();
        thisHistSet->h_electron->Sumw2();
        thisHistSet->h_other->Sumw2();
        thisHistSet->h_offbeam->Sumw2();

        thisHistSet->h_muon->Scale(_config.simscaling);
        thisHistSet->h_proton->Scale(_config.simscaling);
        thisHistSet->h_pion->Scale(_config.simscaling);
        thisHistSet->h_kaon->Scale(_config.simscaling);
        thisHistSet->h_electron->Scale(_config.simscaling);
        thisHistSet->h_other->Scale(_config.simscaling);
        thisHistSet->h_offbeam->Scale(_config.offbeamscaling);
      }
    }

  };

  void HistogramHandler::MakeStackedTrackHistogramsAndSave(std::vector< std::vector<trackhists_1d*> > hists){

    numusel::HistogramHandler _histohandler;
    numusel::Configuration _config;

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    TPad *topPad = new TPad("topPad", "", 0.005, 0.3, 0.995, 0.995);
    TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
    topPad->SetTopMargin(0.28);
    topPad->SetBottomMargin(0.005);
    bottomPad->SetTopMargin(0.005);
    bottomPad->SetBottomMargin(0.25);
    topPad->Draw();
    bottomPad->Draw();

    for (int i_st = 0; i_st < hists.size(); i_st++){
      for (int i_pl = 0; i_pl < hists.at(i_st).size(); i_pl++){

        if (_histohandler.histoMakeTrackPlot.at(i_pl) == false) continue;

        trackhists_1d* thisHistSet = hists.at(i_st).at(i_pl);

        if (thisHistSet->h_onbeam->Integral() == 0) continue;

        _histohandler.StyleTrackHistograms(thisHistSet);

        topPad->cd();

        THStack *hs = new THStack("hs", std::string(_histohandler.histoLabels.at(i_pl)).c_str());
        hs->Add(thisHistSet->h_offbeam);
        hs->Add(thisHistSet->h_muon);
        hs->Add(thisHistSet->h_proton);
        hs->Add(thisHistSet->h_pion);
        hs->Add(thisHistSet->h_kaon);
        hs->Add(thisHistSet->h_electron);
        hs->Add(thisHistSet->h_other);

        TH1D *h_tot = (TH1D*)thisHistSet->h_offbeam->Clone("h_tot");
        h_tot->Add(thisHistSet->h_muon);
        h_tot->Add(thisHistSet->h_proton);
        h_tot->Add(thisHistSet->h_pion);
        h_tot->Add(thisHistSet->h_kaon);
        h_tot->Add(thisHistSet->h_electron);
        h_tot->Add(thisHistSet->h_other);

        h_tot->SetFillStyle(3354);
        h_tot->SetFillColor(kBlack);

        h_tot->GetYaxis()->SetRangeUser(0.001, h_tot->GetMaximum()*1.1);
        h_tot->GetYaxis()->SetTitle("Selected Events");
        h_tot->GetYaxis()->SetTitleSize(0.045);
        h_tot->GetYaxis()->SetTitleOffset(1.1);

        h_tot->GetYaxis()->SetRangeUser(0.001, std::max(h_tot->GetMaximum(),thisHistSet->h_onbeam->GetMaximum())*1.25);

        h_tot->Draw();
        hs->Draw("histsame");
        h_tot->Draw("sameE2");

        thisHistSet->h_onbeam->Draw("samepE1");

        // Write POT to TPaveText
        std::ostringstream streamObj;
        streamObj << _config.onbeam_tor860_wcut;
        std::string data_pot = std::string("POT: ")+streamObj.str();
        TPaveText *pt = new TPaveText(0.5, 0.66, 0.89, 0.71,"NDC");
        pt->AddText(data_pot.c_str());
        pt->SetFillColor(kWhite);
        pt->SetTextAlign(32);
        pt->Draw("same");

        std::pair<double, int> chi2ndof = _histohandler.CalcChi2(h_tot, thisHistSet->h_onbeam);
        TString chi2_string = Form("%0.2f", chi2ndof.first);

        std::string chi2ndof_string = std::string("#chi^{2}: ")+std::string(chi2_string.Data())+std::string("/")+std::to_string(chi2ndof.second);
        TPaveText *pt2 = new TPaveText(0.12, 0.66, 0.5, 0.71,"NDC");
        pt2->AddText(chi2ndof_string.c_str());
        pt2->SetFillColor(kWhite);
        pt2->SetTextAlign(12);
        pt2->Draw("same");

        TLegend *leg_1 = new TLegend(0.1, 0.75, 0.5, 0.95);
        leg_1->AddEntry(thisHistSet->h_muon,   Form("Muon, %g entries", thisHistSet->h_muon->Integral()));
        leg_1->AddEntry(thisHistSet->h_proton, Form("Proton, %g entries", thisHistSet->h_proton->Integral()));
        leg_1->AddEntry(thisHistSet->h_pion,   Form("Pion, %g entries", thisHistSet->h_pion->Integral()));
        leg_1->AddEntry(thisHistSet->h_kaon,   Form("Kaon, %g entries", thisHistSet->h_kaon->Integral()));

        TLegend *leg_2 = new TLegend(0.5, 0.75, 0.9, 0.95);
        leg_2->AddEntry(thisHistSet->h_electron, Form("Electron, %g entries", thisHistSet->h_electron->Integral()));
        leg_2->AddEntry(thisHistSet->h_other,    Form("Other, %g entries", thisHistSet->h_other->Integral()));
        leg_2->AddEntry(thisHistSet->h_offbeam,  Form("Off-beam Data, %g entries", thisHistSet->h_offbeam->Integral()));
        leg_2->AddEntry(thisHistSet->h_onbeam,   Form("On-beam Data, %g entries", thisHistSet->h_onbeam->Integral()));

        leg_1->SetLineWidth(0);
        leg_1->SetFillStyle(0);
        leg_1->Draw("same");
        leg_2->SetLineWidth(0);
        leg_2->SetFillStyle(0);
        leg_2->Draw("same");

        h_tot->GetYaxis()->Draw("same");

        gPad->RedrawAxis();

        bottomPad->cd();
        gStyle->SetOptStat(0);
        TH1D* h_onbeam_clone = (TH1D*)thisHistSet->h_onbeam->Clone("h_onbeam_clone");
        h_onbeam_clone->Sumw2();
        h_onbeam_clone->Add(h_tot, -1);
        h_onbeam_clone->Divide(h_tot);

        TF1* f0 = new TF1("f0", "[0]", 
            _histohandler.histoBins.at(i_pl).at(1), 
            _histohandler.histoBins.at(i_pl).at(2));

        f0->SetTitle(_histohandler.histoLabels.at(i_pl).c_str());
        f0->SetParameter(0,0.0);
        f0->SetNpx(100);
        f0->SetLineStyle(2);
        f0->SetLineColor(kGray+2);
        f0->GetYaxis()->SetRangeUser(-1, 1);
        f0->GetXaxis()->SetRangeUser(_histohandler.histoBins.at(i_pl).at(1), _histohandler.histoBins.at(i_pl).at(2));
        f0->GetXaxis()->SetLabelSize(0.08);
        f0->GetYaxis()->SetLabelSize(0.08);
        f0->GetXaxis()->SetTitleSize(0.1);
        f0->GetXaxis()->SetTitleOffset(1.0);
        f0->GetYaxis()->SetNdivisions(303);
        f0->GetYaxis()->SetTitle("#frac{(ONBEAM-STACK)}{STACK}");
        f0->GetYaxis()->SetTitleSize(0.1);
        f0->GetYaxis()->SetTitleOffset(0.4);
        f0->Draw();

        TF1* f50low = new TF1("f50low", "[0]", 
            _histohandler.histoBins.at(i_pl).at(1), 
            _histohandler.histoBins.at(i_pl).at(2) + (_histohandler.histoBins.at(i_pl).at(2)-_histohandler.histoBins.at(i_pl).at(1))/(2*_histohandler.histoBins.at(i_pl).at(0)));

        f50low->SetParameter(0,0.5);
        f50low->SetTitle(_histohandler.histoLabels.at(i_pl).c_str());

        f50low->SetLineStyle(2);
        f50low->SetLineColor(kGray);
        f50low->GetYaxis()->SetRangeUser(-1, 1);
        f50low->GetXaxis()->SetLabelSize(0.08);
        f50low->GetYaxis()->SetLabelSize(0.08);
        f50low->SetNpx(99);
        f50low->Draw("same");

        TF1* f50high = new TF1("f50high", "[0]", 
            _histohandler.histoBins.at(i_pl).at(1), 
            _histohandler.histoBins.at(i_pl).at(2) + (_histohandler.histoBins.at(i_pl).at(2)-_histohandler.histoBins.at(i_pl).at(1))/(2*_histohandler.histoBins.at(i_pl).at(0)));

        f50high->SetParameter(0,-0.5);
        f50high->SetTitle(_histohandler.histoLabels.at(i_pl).c_str());

        f50high->SetLineStyle(2);
        f50high->SetLineColor(kGray);
        f50high->GetYaxis()->SetRangeUser(-1, 1);
        f50high->GetXaxis()->SetLabelSize(0.08);
        f50high->GetYaxis()->SetLabelSize(0.08);
        f50high->SetNpx(99);
        f50high->Draw("same");

        h_onbeam_clone->Draw("same");


        c1->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames.at(i_pl)
              +std::string("Track_stage")
              +std::to_string(i_st)
              +std::string(".pdf")).c_str());

        c1->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames.at(i_pl)
              +std::string("Track_stage")
              +std::to_string(i_st)
              +std::string(".png")).c_str());

       }

    }

  };



  void HistogramHandler::MakeStackedHistogramsAndSave(std::vector< std::vector<hists_1d*> > hists){

    numusel::HistogramHandler _histohandler;
    numusel::Configuration    _config;

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
    TPad *topPad = new TPad("topPad", "", 0.005, 0.3, 0.995, 0.995);
    TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
    topPad->SetTopMargin(0.28);
    topPad->SetBottomMargin(0.005);
    bottomPad->SetTopMargin(0.005);
    bottomPad->SetBottomMargin(0.25);
    topPad->Draw();
    bottomPad->Draw();

    for (int i_st = 0; i_st < hists.size(); i_st++){
      for (int i_pl = 0; i_pl < hists.at(i_st).size(); i_pl++){

        hists_1d* thisHistSet = hists.at(i_st).at(i_pl);

        if (thisHistSet->h_onbeam->Integral() == 0) continue;

        _histohandler.StyleHistograms(thisHistSet);

        topPad->cd();

        THStack *hs = new THStack("hs", std::string(_histohandler.histoLabels.at(i_pl)).c_str());
        hs->Add(thisHistSet->h_offbeam);
        hs->Add(thisHistSet->h_mccosmic);
        hs->Add(thisHistSet->h_mcmixed);
        hs->Add(thisHistSet->h_mcoofv);
        hs->Add(thisHistSet->h_mcnc);
        hs->Add(thisHistSet->h_mcnuenuebar);
        hs->Add(thisHistSet->h_mcnumubar);
        hs->Add(thisHistSet->h_mcnumuccother);
        hs->Add(thisHistSet->h_mcnumucc0pinp);

        TH1D *h_tot = (TH1D*)thisHistSet->h_offbeam->Clone("h_tot");
        h_tot->Add(thisHistSet->h_mccosmic);
        h_tot->Add(thisHistSet->h_mcmixed);
        h_tot->Add(thisHistSet->h_mcoofv);
        h_tot->Add(thisHistSet->h_mcnc);
        h_tot->Add(thisHistSet->h_mcnuenuebar);
        h_tot->Add(thisHistSet->h_mcnumubar);
        h_tot->Add(thisHistSet->h_mcnumuccother);
        h_tot->Add(thisHistSet->h_mcnumucc0pinp);

        h_tot->SetFillStyle(3354);
        h_tot->SetFillColor(kBlack);

        h_tot->GetYaxis()->SetRangeUser(0.001, h_tot->GetMaximum()*1.1);
        h_tot->GetYaxis()->SetTitle("Selected Events");
        h_tot->GetYaxis()->SetTitleSize(0.045);
        h_tot->GetYaxis()->SetTitleOffset(1.1);

        h_tot->GetYaxis()->SetRangeUser(0.001, std::max(h_tot->GetMaximum(),thisHistSet->h_onbeam->GetMaximum())*1.25);

        h_tot->Draw();
        hs->Draw("histsame");
        h_tot->Draw("sameE2");

        thisHistSet->h_onbeam->Draw("samepE1");

        // Write POT to TPaveText
        std::ostringstream streamObj;
        streamObj << _config.onbeam_tor860_wcut;
        std::string data_pot = std::string("POT: ")+streamObj.str();
        TPaveText *pt = new TPaveText(0.5, 0.66, 0.89, 0.71,"NDC");
        pt->AddText(data_pot.c_str());
        pt->SetFillColor(kWhite);
        pt->SetTextAlign(32);
        pt->Draw("same");

        std::pair<double, int> chi2ndof = _histohandler.CalcChi2(h_tot, thisHistSet->h_onbeam);
        TString chi2_string = Form("%0.2f", chi2ndof.first);

        std::string chi2ndof_string = std::string("#chi^{2}: ")+std::string(chi2_string.Data())+std::string("/")+std::to_string(chi2ndof.second);
        TPaveText *pt2 = new TPaveText(0.12, 0.66, 0.5, 0.71,"NDC");
        pt2->AddText(chi2ndof_string.c_str());
        pt2->SetFillColor(kWhite);
        pt2->SetTextAlign(12);
        pt2->Draw("same");

        TLegend *leg_1 = new TLegend(0.1, 0.75, 0.5, 0.98);
        leg_1->AddEntry(thisHistSet->h_mccosmic, Form("Cosmic, %g entries", thisHistSet->h_mccosmic->Integral()));
        leg_1->AddEntry(thisHistSet->h_mcmixed, Form("Mixed, %g entries", thisHistSet->h_mcmixed->Integral()));
        leg_1->AddEntry(thisHistSet->h_mcoofv, Form("OOFV, %g entries", thisHistSet->h_mcoofv->Integral()));
        leg_1->AddEntry(thisHistSet->h_mcnc, Form("NC, %g entries", thisHistSet->h_mcnc->Integral()));
        leg_1->AddEntry(thisHistSet->h_mcnuenuebar, Form("#nu_{e}/#bar{#nu_{e}}, %g entries", thisHistSet->h_mcnuenuebar->Integral()));

        TLegend *leg_2 = new TLegend(0.5, 0.75, 0.9, 0.98);
        leg_2->AddEntry(thisHistSet->h_mcnumubar, Form("#bar{#nu_{#mu}}, %g entries", thisHistSet->h_mcnumubar->Integral()));
        leg_2->AddEntry(thisHistSet->h_mcnumuccother, Form("#nu_{#mu}CC-Other, %g entries", thisHistSet->h_mcnumuccother->Integral()));
        leg_2->AddEntry(thisHistSet->h_mcnumucc0pinp, Form("#nu_{#mu}CC0#piNP, %g entries", thisHistSet->h_mcnumucc0pinp->Integral()));
        leg_2->AddEntry(thisHistSet->h_offbeam, Form("Off-beam Data, %g entries", thisHistSet->h_offbeam->Integral()));
        leg_2->AddEntry(thisHistSet->h_onbeam, Form("On-beam Data, %g entries", thisHistSet->h_onbeam->Integral()));

        leg_1->SetLineWidth(0);
        leg_1->SetFillStyle(0);
        leg_1->Draw("same");
        leg_2->SetLineWidth(0);
        leg_2->SetFillStyle(0);
        leg_2->Draw("same");

        h_tot->GetYaxis()->Draw("same");

        gPad->RedrawAxis();

        bottomPad->cd();
        gStyle->SetOptStat(0);
        TH1D* h_onbeam_clone = (TH1D*)thisHistSet->h_onbeam->Clone("h_onbeam_clone");
        h_onbeam_clone->Sumw2();
        h_onbeam_clone->Add(h_tot, -1);
        h_onbeam_clone->Divide(h_tot);

        TF1* f0 = new TF1("f0", "[0]", 
            _histohandler.histoBins.at(i_pl).at(1), 
            _histohandler.histoBins.at(i_pl).at(2));

        f0->SetTitle(_histohandler.histoLabels.at(i_pl).c_str());
        f0->SetParameter(0,0.0);
        f0->SetNpx(100);
        f0->SetLineStyle(2);
        f0->SetLineColor(kGray+2);
        f0->GetYaxis()->SetRangeUser(-1, 1);
        f0->GetXaxis()->SetRangeUser(_histohandler.histoBins.at(i_pl).at(1), _histohandler.histoBins.at(i_pl).at(2));
        f0->GetXaxis()->SetLabelSize(0.08);
        f0->GetYaxis()->SetLabelSize(0.08);
        f0->GetXaxis()->SetTitleSize(0.1);
        f0->GetXaxis()->SetTitleOffset(1.0);
        f0->GetYaxis()->SetNdivisions(303);
        f0->GetYaxis()->SetTitle("#frac{(ONBEAM-STACK)}{STACK}");
        f0->GetYaxis()->SetTitleSize(0.1);
        f0->GetYaxis()->SetTitleOffset(0.4);
        f0->Draw();

        TF1* f50low = new TF1("f50low", "[0]", 
            _histohandler.histoBins.at(i_pl).at(1), 
            _histohandler.histoBins.at(i_pl).at(2) + (_histohandler.histoBins.at(i_pl).at(2)-_histohandler.histoBins.at(i_pl).at(1))/(2*_histohandler.histoBins.at(i_pl).at(0)));

        f50low->SetParameter(0,0.5);
        f50low->SetTitle(_histohandler.histoLabels.at(i_pl).c_str());

        f50low->SetLineStyle(2);
        f50low->SetLineColor(kGray);
        f50low->GetYaxis()->SetRangeUser(-1, 1);
        f50low->GetXaxis()->SetLabelSize(0.08);
        f50low->GetYaxis()->SetLabelSize(0.08);
        f50low->SetNpx(99);
        f50low->Draw("same");

        TF1* f50high = new TF1("f50high", "[0]", 
            _histohandler.histoBins.at(i_pl).at(1), 
            _histohandler.histoBins.at(i_pl).at(2) + (_histohandler.histoBins.at(i_pl).at(2)-_histohandler.histoBins.at(i_pl).at(1))/(2*_histohandler.histoBins.at(i_pl).at(0)));

        f50high->SetParameter(0,-0.5);
        f50high->SetTitle(_histohandler.histoLabels.at(i_pl).c_str());

        f50high->SetLineStyle(2);
        f50high->SetLineColor(kGray);
        f50high->GetYaxis()->SetRangeUser(-1, 1);
        f50high->GetXaxis()->SetLabelSize(0.08);
        f50high->GetYaxis()->SetLabelSize(0.08);
        f50high->SetNpx(99);
        f50high->Draw("same");

        h_onbeam_clone->Draw("same");


        c1->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string(".pdf")).c_str());

        c1->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string(".png")).c_str());

        TCanvas *c2 = new TCanvas();
        c2->cd();

        TH1D* h_purity = (TH1D*)thisHistSet->h_mcnumucc0pinp->Clone("h_purity");
        h_tot->Sumw2();
        h_purity->Divide(h_tot);
        h_purity->GetYaxis()->SetRangeUser(0,1);

        h_purity->SetMarkerStyle(20);
        h_purity->SetMarkerSize(0.3);
        h_purity->SetLineWidth(1);
        h_purity->SetLineColor(TColor::GetColor(12,128,243));
        h_purity->SetMarkerColor(TColor::GetColor(12,128,243));
        h_purity->GetYaxis()->SetTitle("Purity");
        h_purity->Draw("E1p][same");

        gPad->RedrawAxis();
        c2->SetGridy();

        c2->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames.at(i_pl)
              +std::string("_purity_stage")
              +std::to_string(i_st)
              +std::string(".pdf")).c_str());

        c2->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames.at(i_pl)
              +std::string("_purity_stage")
              +std::to_string(i_st)
              +std::string(".png")).c_str());

      }

    }

  };

  void HistogramHandler::Make2DHistogramsAndSave(std::vector< std::vector<hists_2d*> > hists){

    numusel::HistogramHandler _histohandler;

    for (int i_st = 0; i_st < hists.size(); i_st++){
      for (int i_pl = 0; i_pl < hists.at(i_st).size(); i_pl++){

        hists_2d* thisHistSet = hists.at(i_st).at(i_pl);
        if (thisHistSet->h_onbeam->Integral() == 0) continue;

        TCanvas *c_mc = new TCanvas();
        hists.at(i_st).at(i_pl)->h_mc->Draw("colz");

        // if it's the plot of E_nu versus E_dep then fit a straight line
        // to it
        if (_histohandler.histoNames_2D.at(i_pl) == "h_trueenu_recoenu"){

          hists.at(i_st).at(i_pl)->h_mc->Fit("pol1", "", "", 0.4, 1.0);
          float par0 = hists.at(i_st).at(i_pl)->h_mc->GetFunction("pol1")->GetParameter(0);
          float par1 = hists.at(i_st).at(i_pl)->h_mc->GetFunction("pol1")->GetParameter(1);

          TF1* func = new TF1("func", "[0]*x+[1]", 0, 3);
          func->SetParameters(par1, par0);
          func->Draw("same");

        }

        c_mc->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames_2D.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string("_simulation.png")).c_str());

        TCanvas *c_onbeam = new TCanvas();
        hists.at(i_st).at(i_pl)->h_onbeam->Draw("colz");

        c_onbeam->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames_2D.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string("_onbeam.png")).c_str());

        TCanvas *c_offbeam = new TCanvas();
        hists.at(i_st).at(i_pl)->h_offbeam->Draw("colz");

        c_offbeam->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.histoNames_2D.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string("_offbeam.png")).c_str());
      }
    }
  };

  void HistogramHandler::MakeEfficiencyHistogramsAndSave(std::vector<std::vector<eff_1d*>> effplots){

    numusel::HistogramHandler _histohandler;

    TFile* outfile = new TFile("efficiencyPlots.root", "update");
    outfile->cd();

    TCanvas *c1 = new TCanvas();

    for (int i_st = 0; i_st < effplots.size(); i_st++){
      for (int i_pl = 0; i_pl < effplots.at(i_st).size(); i_pl++){

        eff_1d* eff1d = effplots.at(i_st).at(i_pl);

        TH1D* effplot = (TH1D*)eff1d->h_num->Clone("effplot");
        effplot->Sumw2();
        effplot->Divide(effplots.at(0).at(i_pl)->h_denom);
        effplot->SetLineColor(TColor::GetColor(166,217,106));
        effplot->SetFillColor(TColor::GetColor(166,217,106));
        effplot->GetYaxis()->SetRangeUser(0,1.0);
        effplot->GetYaxis()->SetTitle("Efficiency");

        c1->SetGridy();

        effplot->Draw("E1");

        effplot->SetName(std::string(
              _histohandler.effNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)).c_str());
        effplot->Write();

        c1->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.effNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string(".pdf")).c_str());

        c1->SaveAs(std::string(
              std::string("plots/")
              +_histohandler.effNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string(".png")).c_str());

      }
    }
  };
}
