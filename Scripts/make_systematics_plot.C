void make_systematics_plot(){

    bool flux_or_xsec = 1; // flux = 0, xsec = 1;

    int n_uni;
    if (flux_or_xsec == 0)
        n_uni = 1000;
    else if (flux_or_xsec == 1)
        n_uni = 50;

    TFile *_file0 = new TFile("/uboone/data/users/alister1/SBNFit_files/filtered_events_mc_numu_sbnfit.root");
    TTree *t_cc0pinp = (TTree*)_file0->Get("numu_sbnfit_cc0pinp");
    TTree *t_ccother = (TTree*)_file0->Get("numu_sbnfit_ccother");
    TTree *t_oofv = (TTree*)_file0->Get("numu_sbnfit_oofv");
    TTree *t_cosmic = (TTree*)_file0->Get("numu_sbnfit_cosmic");
    TTree *t_nc = (TTree*)_file0->Get("numu_sbnfit_nc");
    TTree *t_mixed = (TTree*)_file0->Get("numu_sbnfit_mixed");
    TTree *t_nue = (TTree*)_file0->Get("numu_sbnfit_nue");
    TTree *t_anumu = (TTree*)_file0->Get("numu_sbnfit_anumu");

    std::vector<TTree*> tree_vector = {
        t_cc0pinp,
        t_ccother,
        t_oofv,
        t_cosmic,
        t_nc,
        t_mixed,
        t_nue,
        t_anumu
    };

    std::vector< std::map< std::string, std::vector<double> >* > weights_vector = {
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr
    };

    std::vector<double> recoenergy_vector = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    };

    t_cc0pinp->SetBranchAddress("weights", &(weights_vector.at(0)));
    t_ccother->SetBranchAddress("weights", &(weights_vector.at(1)));
    t_oofv->SetBranchAddress("weights", &(weights_vector.at(2)));
    t_cosmic->SetBranchAddress("weights", &(weights_vector.at(3)));
    t_nc->SetBranchAddress("weights", &(weights_vector.at(4)));
    t_mixed->SetBranchAddress("weights", &(weights_vector.at(5)));
    t_nue->SetBranchAddress("weights", &(weights_vector.at(6)));
    t_anumu->SetBranchAddress("weights", &(weights_vector.at(7)));
    t_cc0pinp->SetBranchAddress("reconstructed_neutrino_energy_calib", &recoenergy_vector.at(0)); 
    t_ccother->SetBranchAddress("reconstructed_neutrino_energy_calib", &recoenergy_vector.at(1)); 
    t_oofv->SetBranchAddress("reconstructed_neutrino_energy_calib",  &recoenergy_vector.at(2)); 
    t_cosmic->SetBranchAddress("reconstructed_neutrino_energy_calib",  &recoenergy_vector.at(3)); 
    t_nc->SetBranchAddress("reconstructed_neutrino_energy_calib",  &recoenergy_vector.at(4)); 
    t_mixed->SetBranchAddress("reconstructed_neutrino_energy_calib", &recoenergy_vector.at(5)); 
    t_nue->SetBranchAddress("reconstructed_neutrino_energy_calib",  &recoenergy_vector.at(6)); 
    t_anumu->SetBranchAddress("reconstructed_neutrino_energy_calib", &recoenergy_vector.at(7)); 

    //
    // cv universe!
    //

    // define histograms
    TH1D* h_cc0pinp_cv = new TH1D("h_cc0pinp_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);
    TH1D* h_ccother_cv = new TH1D("h_ccother_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);
    TH1D* h_oofv_cv = new TH1D("h_oofv_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);
    TH1D* h_cosmic_cv = new TH1D("h_cosmic_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);
    TH1D* h_nc_cv = new TH1D("h_nc_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);
    TH1D* h_mixed_cv = new TH1D("h_mixed_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);
    TH1D* h_nue_cv = new TH1D("h_nue_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);
    TH1D* h_anumu_cv = new TH1D("h_anumu_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);

    // make nominal histogram
    t_cc0pinp->Draw("reconstructed_neutrino_energy_calib >> h_cc0pinp_cv");
    t_ccother->Draw("reconstructed_neutrino_energy_calib >> h_ccother_cv");
    t_oofv->Draw("reconstructed_neutrino_energy_calib >> h_oofv_cv");
    t_cosmic->Draw("reconstructed_neutrino_energy_calib >> h_cosmic_cv");
    t_nc->Draw("reconstructed_neutrino_energy_calib >> h_nc_cv");
    t_mixed->Draw("reconstructed_neutrino_energy_calib >> h_mixed_cv");
    t_nue->Draw("reconstructed_neutrino_energy_calib >> h_nue_cv");
    t_anumu->Draw("reconstructed_neutrino_energy_calib >> h_anumu_cv");

    TH1D *tot_cv = new TH1D("tot_cv", ";Reconstructed neutrino energy (GeV);", 25, 0, 3);

    tot_cv->Add(h_cc0pinp_cv);
    tot_cv->Add(h_ccother_cv);
    tot_cv->Add(h_oofv_cv);
    tot_cv->Add(h_cosmic_cv);
    tot_cv->Add(h_nc_cv);
    tot_cv->Add(h_mixed_cv);
    tot_cv->Add(h_nue_cv);
    tot_cv->Add(h_anumu_cv);

    //
    // alt universes
    //

    // set up vectors of histograms
    std::vector<TH1D*> h_cc0pinp_uni;
    std::vector<TH1D*> h_ccother_uni;
    std::vector<TH1D*> h_oofv_uni;
    std::vector<TH1D*> h_cosmic_uni;
    std::vector<TH1D*> h_nc_uni;
    std::vector<TH1D*> h_mixed_uni;
    std::vector<TH1D*> h_nue_uni;
    std::vector<TH1D*> h_anumu_uni;
    std::vector<TH1D*> h_total;

    for (int i = 0; i < n_uni; i++){

        TString uni_no = ("_%i", i);

        h_cc0pinp_uni.push_back(new TH1D((std::string("h_cc0pinp")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));
        h_ccother_uni.push_back(new TH1D((std::string("h_ccother")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));
        h_oofv_uni.push_back(new TH1D((std::string("h_oofv")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));
        h_cosmic_uni.push_back(new TH1D((std::string("h_cosmic")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));
        h_nc_uni.push_back(new TH1D((std::string("h_nc")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));
        h_mixed_uni.push_back(new TH1D((std::string("h_mixed")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));
        h_nue_uni.push_back(new TH1D((std::string("h_nue")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));
        h_anumu_uni.push_back(new TH1D((std::string("h_anumu")+std::string(uni_no)).c_str(), ";Reonstructed neutrino energy (GeV);", 25, 0, 3));  
        h_total.push_back(new TH1D((std::string("h_total")+std::string(uni_no)).c_str(), ";Reconstructed neutrino energy (Gev);", 25, 0, 3));

    }

    std::vector< std::vector<TH1D*> > histogram_vector = {
        h_cc0pinp_uni,
        h_ccother_uni,
        h_oofv_uni,
        h_cosmic_uni,
        h_nc_uni,
        h_mixed_uni,
        h_nue_uni,
        h_anumu_uni
    };

    // loop trees

    for (int i_tree = 0; i_tree < tree_vector.size(); i_tree++){

        // now loop over events in the each
        int nEntries = tree_vector.at(i_tree)->GetEntries();
        std::cout << "beginning loop of " << nEntries  << " events" << std::endl;
        for (int i = 0; i < nEntries; i++){

            if(i%1000 == 0)std::cout << i << "/" << nEntries << std::endl;

            tree_vector.at(i_tree)->GetEntry(i);

            std::vector<double> collapsed_weights(n_uni, 1.0);
            for ( auto const weight : (*(weights_vector.at(i_tree)))){
                std::vector< double > weight_uni = weight.second;

                if ((int)weight_uni.size() != n_uni) continue;

                // loop universes
                for (int i_u = 0; i_u < weight_uni.size(); i_u++){
                    collapsed_weights.at(i_u) *= weight_uni.at(i_u);
                }

            }

            for (int i_u = 0; i_u < collapsed_weights.size(); i_u++){
                histogram_vector.at(i_tree).at(i_u)->Fill(recoenergy_vector.at(i_tree), collapsed_weights.at(i_u));
            }
        }
        std::cout << "end loop of " << nEntries << " events" << std::endl;
    }
    tot_cv->GetYaxis()->SetRangeUser(0, tot_cv->GetMaximum()*2.0);
    tot_cv->Draw();

    // loop universes
    for (int i = 0; i < n_uni; i++){
    
        //loop histogram
        for (int i_hist = 0; i_hist < histogram_vector.size()-1; i_hist++){
            histogram_vector.at(histogram_vector.size()-1).at(i)->Add(histogram_vector.at(i_hist).at(i));
        }

    }

    for (int i = 0; i < n_uni; i++){
        histogram_vector.at(histogram_vector.size()-1).at(i)->SetLineColor(kGreen+1);
        histogram_vector.at(histogram_vector.size()-1).at(i)->Draw("same");
    }
    tot_cv->DrawClone("same");


}
