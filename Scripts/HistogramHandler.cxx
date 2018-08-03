#include "HistogramHandler.h"

namespace numusel{

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

  void HistogramHandler::ScaleHistograms(hists_1d* hists){

    hists->h_mccosmic->Sumw2();
    hists->h_mcmixed->Sumw2();
    hists->h_mcoofv->Sumw2();
    hists->h_mcnc->Sumw2();
    hists->h_mcnuenuebar->Sumw2();
    hists->h_mcnumubar->Sumw2();
    hists->h_mcnumuccother->Sumw2();
    hists->h_mcnumucc0pinp->Sumw2();
    hists->h_offbeam->Sumw2();

    hists->h_mccosmic->Scale(_config.simscaling);
    hists->h_mcmixed->Scale(_config.simscaling);
    hists->h_mcoofv->Scale(_config.simscaling);
    hists->h_mcnc->Scale(_config.simscaling);
    hists->h_mcnuenuebar->Scale(_config.simscaling);
    hists->h_mcnumubar->Scale(_config.simscaling);
    hists->h_mcnumuccother->Scale(_config.simscaling);
    hists->h_mcnumucc0pinp->Scale(_config.simscaling);
    hists->h_offbeam->Scale(_config.offbeamscaling);

  };

  void HistogramHandler::MakeStackedHistogramsAndSave(std::vector< std::vector<hists_1d*> > hists){

    numusel::HistogramHandler _histhandler;
    
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

        _histhandler.StyleHistograms(thisHistSet);
        _histhandler.ScaleHistograms(thisHistSet);

        topPad->cd();

        THStack *hs = new THStack("hs", std::string(_histhandler.histoLabels.at(i_pl)).c_str());
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

        h_tot->Draw();
        hs->Draw("histsame");
        h_tot->Draw("sameE2");

        thisHistSet->h_onbeam->Draw("samepE1");

        TLegend *leg_1 = new TLegend(0.1, 0.75, 0.5, 0.98);
        leg_1->AddEntry(thisHistSet->h_mccosmic, "Cosmic");
        leg_1->AddEntry(thisHistSet->h_mcmixed, "Mixed");
        leg_1->AddEntry(thisHistSet->h_mcoofv, "OOFV");
        leg_1->AddEntry(thisHistSet->h_mcnc, "NC");
        leg_1->AddEntry(thisHistSet->h_mcnuenuebar, "#nu_{e}/#bar{#nu_{e}}");

        TLegend *leg_2 = new TLegend(0.5, 0.75, 0.9, 0.98);
        leg_2->AddEntry(thisHistSet->h_mcnumubar, "#bar{#nu_{#mu}}");
        leg_2->AddEntry(thisHistSet->h_mcnumuccother, "#nu_{#mu}CC-Other");
        leg_2->AddEntry(thisHistSet->h_mcnumucc0pinp, "#nu_{#mu}CC0#piNP");
        leg_2->AddEntry(thisHistSet->h_offbeam, "Off-beam Data");
        leg_2->AddEntry(thisHistSet->h_onbeam, "On-beam Data");

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
            _histhandler.histoBins.at(i_pl).at(1), 
            _histhandler.histoBins.at(i_pl).at(2));

        f0->SetTitle(_histhandler.histoLabels.at(i_pl).c_str());
        f0->SetParameter(0,0.0);
        f0->SetNpx(100);
        f0->SetLineStyle(2);
        f0->SetLineColor(kGray+2);
        f0->GetYaxis()->SetRangeUser(-1, 1);
        f0->GetXaxis()->SetRangeUser(_histhandler.histoBins.at(i_pl).at(1), _histhandler.histoBins.at(i_pl).at(2));
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
            _histhandler.histoBins.at(i_pl).at(1), 
            _histhandler.histoBins.at(i_pl).at(2) + (_histhandler.histoBins.at(i_pl).at(2)-_histhandler.histoBins.at(i_pl).at(1))/(2*_histhandler.histoBins.at(i_pl).at(0)));

        f50low->SetParameter(0,0.5);
        f50low->SetTitle(_histhandler.histoLabels.at(i_pl).c_str());

        f50low->SetLineStyle(2);
        f50low->SetLineColor(kGray);
        f50low->GetYaxis()->SetRangeUser(-1, 1);
        f50low->GetXaxis()->SetLabelSize(0.08);
        f50low->GetYaxis()->SetLabelSize(0.08);
        f50low->SetNpx(99);
        f50low->Draw("same");

        TF1* f50high = new TF1("f50high", "[0]", 
            _histhandler.histoBins.at(i_pl).at(1), 
            _histhandler.histoBins.at(i_pl).at(2) + (_histhandler.histoBins.at(i_pl).at(2)-_histhandler.histoBins.at(i_pl).at(1))/(2*_histhandler.histoBins.at(i_pl).at(0)));

        f50high->SetParameter(0,-0.5);
        f50high->SetTitle(_histhandler.histoLabels.at(i_pl).c_str());

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
              +_histhandler.histoNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string(".pdf")).c_str());

        c1->SaveAs(std::string(
              std::string("plots/")
              +_histhandler.histoNames.at(i_pl)
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
              +_histhandler.histoNames.at(i_pl)
              +std::string("_purity_stage")
              +std::to_string(i_st)
              +std::string(".pdf")).c_str());

        c2->SaveAs(std::string(
              std::string("plots/")
              +_histhandler.histoNames.at(i_pl)
              +std::string("_purity_stage")
              +std::to_string(i_st)
              +std::string(".png")).c_str());

      }

    }

  };

  void HistogramHandler::MakeEfficiencyHistogramsAndSave(std::vector<std::vector<eff_1d*>> effplots){

    numusel::HistogramHandler _histhandler;

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

        c1->SaveAs(std::string(
              std::string("plots/")
              +_histhandler.effpurNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string(".pdf")).c_str());

        c1->SaveAs(std::string(
              std::string("plots/")
              +_histhandler.effpurNames.at(i_pl)
              +std::string("_stage")
              +std::to_string(i_st)
              +std::string(".png")).c_str());

      }
    }
  };
}
