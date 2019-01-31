void make_data_trackfrac_phi_plot(){

    std::string varToPlot("track_phi");
    int nbins = 9;
    float bins[] = {-3.15, -3.1499, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.1499, 3.15};
 
    // files
    TFile * f_data = new TFile("/uboone/data/users/alister1/numuSelection/files/190125/onbeam_selectionInformation.root");
    TFile * f_cv   = new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_cv.root");

    std::vector<TFile*> detvarsfiles={
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_sce.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_larg4bugfix.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dlup.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dldown.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dtup.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dtdown.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_noiseampup.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_noiseampdown.root"),
        new TFile("/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dic.root")
    };

    TTree* t_data = (TTree*)f_data->Get("numuselection/analysis_tree");
    TTree* t_cv   = (TTree*)f_cv->Get("numuselection/analysis_tree");
   
    std::vector<TTree*> detvarstrees;
    for (int i = 0; i < detvarsfiles.size(); i++){
        detvarstrees.push_back((TTree*)detvarsfiles.at(i)->Get("numuselection/analysis_tree"));
    }

    
    // deal with onbeam
    TH1D* num_onbeam   = new TH1D("num_onbeam",   ";Proton candidate #phi;Track-like fraction",nbins, bins);
    TH1D* denom_onbeam = new TH1D("denom_onbeam", ";Proton candidate #phi;Track-like fraction",nbins, bins);

    t_data->Draw(std::string(varToPlot+">> num_onbeam").c_str(), "log(noBragg_fwd_mip/max(bragg_fwd_p, bragg_bwd_p)) < -1 && pfp_pdgCode == 13");
    t_data->Draw(std::string(varToPlot+">> denom_onbeam").c_str(), "log(noBragg_fwd_mip/max(bragg_fwd_p, bragg_bwd_p)) < -1");
    
    num_onbeam->Sumw2();
    denom_onbeam->Sumw2();
    num_onbeam->Divide(denom_onbeam);

    // deal with cv
    TH1D* num_cv = new TH1D("num_cv", ";Proton candidate #phi;Track-like fraction",nbins,bins );
    TH1D* denom_cv = new TH1D("denom_cv", ";Proton candidate #phi;Track-like fraction",nbins,bins );

    t_cv->Draw(std::string(varToPlot+">> num_cv").c_str(), "log(noBragg_fwd_mip/max(bragg_fwd_p, bragg_bwd_p)) < -1 && pfp_pdgCode == 13");
    t_cv->Draw(std::string(varToPlot+">> denom_cv").c_str(), "log(noBragg_fwd_mip/max(bragg_fwd_p, bragg_bwd_p)) < -1");
    num_cv->Sumw2();
    denom_cv->Sumw2();

    num_cv->Divide(denom_cv);
    
    // and deal with the detector variations

    std::vector<TH1D*> ratioHistos;
    for (int i = 0; i < detvarstrees.size(); i++){

       TH1D* h_nom = new TH1D("h_nom", ";;", nbins, bins);
       TH1D* h_denom = new TH1D("h_denom", ";;", nbins, bins);
       
       detvarstrees.at(i)->Draw(std::string(varToPlot+">> h_nom").c_str(), "log(noBragg_fwd_mip/max(bragg_fwd_p, bragg_bwd_p)) < -1 && pfp_pdgCode == 13");
       detvarstrees.at(i)->Draw(std::string(varToPlot+">> h_denom").c_str(), "log(noBragg_fwd_mip/max(bragg_fwd_p, bragg_bwd_p)) < -1");
       
       h_nom->Sumw2();
       h_denom->Sumw2();

       h_nom->Divide(h_denom);

       ratioHistos.push_back((TH1D*)h_nom->Clone(std::string(Form("h%i", i)).c_str()));

       h_nom->Delete();
       h_denom->Delete();
    }


    float errors[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    for (int i = 0; i < detvarstrees.size(); i++){

       for (int j = 0; j < num_cv->GetNbinsX(); j++){
           errors[j] += std::pow(num_cv->GetBinContent(j)-ratioHistos.at(i)->GetBinContent(j), 2);
       }

    }

    for (int i = 0 ; i < num_cv->GetNbinsX(); i++){
        num_cv->SetBinError(i, std::sqrt(std::pow(num_cv->GetBinError(i),2)+errors[i]));
    }

    TCanvas *c1 = new TCanvas();

    num_cv->SetLineColor(kGreen+1);
    num_cv->SetFillColor(kGreen-6);
    num_cv->SetFillStyle(3345);
    num_cv->SetMarkerStyle(20);
    num_cv->SetMarkerColor(kGreen+1);
    num_cv->SetMarkerSize(0.6);
    num_cv->GetYaxis()->SetRangeUser(0,1);
    num_cv->SetBinContent(1, num_cv->GetBinContent(2));
    num_cv->SetBinError(1, num_cv->GetBinError(2)*1.3);
    num_cv->SetBinContent(nbins, num_cv->GetBinContent(nbins-1));
    num_cv->SetBinError(nbins, num_cv->GetBinError(8)*1.3);
    num_cv->Draw("E3");
    num_cv->DrawClone("psameX0E1");
    num_onbeam->SetMarkerStyle(20);
    num_onbeam->SetMarkerSize(0.6);
    num_onbeam->Draw("sameE1");

    TLegend *leg = new TLegend(0.6, 0.18, 0.85, 0.28);
    leg->AddEntry(num_onbeam, "On-beam data");
    leg->AddEntry(num_cv, "CV + stat. + syst.");
    leg->Draw("same");

    c1->SaveAs("trackfrac_phi.pdf");

}
