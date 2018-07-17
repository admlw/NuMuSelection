// ROOT
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"

// local
#include "EventCategoriser.h"
#include "DataTypes.h"
#include "HistogramsToProduce.h"
#include "AnalysisCuts.h"
#include "SelectionMaker.h"

// cpp
#include <iostream>
#include <bitset>
#include <string>

// vars
std::string s_onbeam     = "/uboone/data/users/alister1/numuSelection/18-07-16-NuMuSelection/numusel_onbeam.root";
std::string s_offbeam    = "/uboone/data/users/alister1/numuSelection/18-07-16-NuMuSelection/numusel_offbeam.root";
std::string s_simulation = "/uboone/data/users/alister1/numuSelection/18-07-16-NuMuSelection/numusel_bnbcos.root";

double offbeamscaling = 1.05;
double simscaling = 2.3; 

// there are four stages to the selection
// * no cuts (pure UBXSec)
// * n tracks cut
// * n showers cut
// * pid cut
int n_stages = 4; 

std::vector<std::vector<hists_1d*>> plots_to_make;

void setTreeVars(TTree* tree, var_list* varstoset, bool isSimulation){

  tree->SetBranchAddress("isUBXSecSelected" , &(varstoset->isUBXSecSelected));
  tree->SetBranchAddress("nSelectedTracks"  , &(varstoset->nSelectedTracks));
  tree->SetBranchAddress("nSelectedShowers" , &(varstoset->nSelectedShowers));
  tree->SetBranchAddress("bragg_fwd_p"      , &(varstoset->bragg_fwd_p));
  tree->SetBranchAddress("bragg_bwd_p"      , &(varstoset->bragg_bwd_p));
  tree->SetBranchAddress("noBragg_fwd_mip"  , &(varstoset->noBragg_fwd_mip));
  tree->SetBranchAddress("track_length"     , &(varstoset->track_length));
  tree->SetBranchAddress("vertex_x"         , &(varstoset->vertex_x));
  tree->SetBranchAddress("vertex_y"         , &(varstoset->vertex_y));
  tree->SetBranchAddress("vertex_z"         , &(varstoset->vertex_z));

  if (isSimulation){
    tree->SetBranchAddress("isBeamNeutrino"   , &(varstoset->isBeamNeutrino));
    tree->SetBranchAddress("isCosmic"         , &(varstoset->isCosmic));
    tree->SetBranchAddress("isMixed"          , &(varstoset->isMixed));
    tree->SetBranchAddress("isInFV"           , &(varstoset->isInFV));
    tree->SetBranchAddress("true_nu_ccnc"     , &(varstoset->true_nu_ccnc));
    tree->SetBranchAddress("true_match_pdg"   , &(varstoset->true_match_pdg));
    tree->SetBranchAddress("true_mcp_pdg"     , &(varstoset->true_mcp_pdg));
    tree->SetBranchAddress("true_mcp_process" , &(varstoset->true_mcp_process));
  }

};

void FillHistMC(hists_1d* h1d, std::vector<double> variable, std::bitset<8> eventCat){

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

void FillHistOnBeam(hists_1d* h1d, std::vector<double> variable){

  for (int i = 0; i < variable.size(); i++){
    h1d->h_onbeam->Fill(variable.at(i));
  }

};

void FillHistOffBeam(hists_1d* h1d, std::vector<double> variable){

  for (int i = 0; i < variable.size(); i++){
    h1d->h_offbeam->Fill(variable.at(i));
  }

};

void StyleHistograms(hists_1d* hists){

  hists->h_mccosmic->SetFillColor(TColor::GetColor(8,64,129));
  hists->h_mcmixed->SetFillColor(TColor::GetColor(8,104,172));
  hists->h_mcoofv->SetFillColor(TColor::GetColor(78,179,211));
  hists->h_mcnc->SetFillColor(TColor::GetColor(113,1,98));
  hists->h_mcnuenuebar->SetFillColor(TColor::GetColor(166,217,211));
  hists->h_mcnumubar->SetFillColor(TColor::GetColor(0,104,55));
  hists->h_mcnumuccother->SetFillColor(TColor::GetColor(165,0,38));
  hists->h_mcnumucc0pinp->SetFillColor(TColor::GetColor(215,48,39));

  hists->h_mccosmic->SetMarkerColor(TColor::GetColor(8,64,129));
  hists->h_mcmixed->SetMarkerColor(TColor::GetColor(8,104,172));
  hists->h_mcoofv->SetMarkerColor(TColor::GetColor(78,179,211));
  hists->h_mcnc->SetMarkerColor(TColor::GetColor(113,1,98));
  hists->h_mcnuenuebar->SetMarkerColor(TColor::GetColor(166,217,211));
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
}

void ScaleHistograms(hists_1d* hists){

  hists->h_mccosmic->Sumw2();
  hists->h_mcmixed->Sumw2();
  hists->h_mcoofv->Sumw2();
  hists->h_mcnc->Sumw2();
  hists->h_mcnuenuebar->Sumw2();
  hists->h_mcnumubar->Sumw2();
  hists->h_mcnumuccother->Sumw2();
  hists->h_mcnumucc0pinp->Sumw2();
  hists->h_offbeam->Sumw2();

  hists->h_mccosmic->Scale(simscaling);
  hists->h_mcmixed->Scale(simscaling);
  hists->h_mcoofv->Scale(simscaling);
  hists->h_mcnc->Scale(simscaling);
  hists->h_mcnuenuebar->Scale(simscaling);
  hists->h_mcnumubar->Scale(simscaling);
  hists->h_mcnumuccother->Scale(simscaling);
  hists->h_mcnumucc0pinp->Scale(simscaling);
  hists->h_offbeam->Scale(offbeamscaling);

}
void MakeStackedHistogramsAndSave(std::vector< std::vector<hists_1d*> > hists){

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetTopMargin(0.20);

  for (int i_st = 0; i_st < hists.size(); i_st++){
    for (int i_pl = 0; i_pl < hists.at(i_st).size(); i_pl++){

      hists_1d* thisHistSet = hists.at(i_st).at(i_pl);

      StyleHistograms(thisHistSet);
      ScaleHistograms(thisHistSet);


      THStack *hs = new THStack("hs", std::string(histoLabels.at(i_pl)).c_str());
      hs->Add(thisHistSet->h_offbeam);
      hs->Add(thisHistSet->h_mccosmic);
      hs->Add(thisHistSet->h_mcmixed);
      hs->Add(thisHistSet->h_mcoofv);
      hs->Add(thisHistSet->h_mcnc);
      hs->Add(thisHistSet->h_mcnuenuebar);
      hs->Add(thisHistSet->h_mcnumubar);
      hs->Add(thisHistSet->h_mcnumuccother);
      hs->Add(thisHistSet->h_mcnumucc0pinp);

      hs->Draw("hist");

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
      h_tot->Draw("sameE2");

      thisHistSet->h_onbeam->Draw("samepE1");

      TLegend *leg_1 = new TLegend(0.1, 0.82, 0.5, 0.98);
      leg_1->AddEntry(thisHistSet->h_mccosmic, "Cosmic");
      leg_1->AddEntry(thisHistSet->h_mcmixed, "Mixed");
      leg_1->AddEntry(thisHistSet->h_mcoofv, "OOFV");
      leg_1->AddEntry(thisHistSet->h_mcnc, "NC");
      leg_1->AddEntry(thisHistSet->h_mcnuenuebar, "#nu_{e}/#bar{#nu_{e}}");

      TLegend *leg_2 = new TLegend(0.5, 0.82, 0.9, 0.98);
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


      c1->SaveAs(std::string(
            std::string("plots/")
            +histoNames.at(i_pl)
            +std::string("_stage")
            +std::to_string(i_st)
            +std::string(".png")).c_str());
      
    }

  }

}


