void make_systematics_plot(){

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
    c1->cd();
    TH1D* h_cv = new TH1D("h_cv", ";Reconstructed Neutrino Energy (GeV);", 20, 0, 3);
    TH1D* universes[5000];
    for (int i = 0; i < 5000; i++){

        TString uni_name = Form("uni_%i", i);
        universes[i] = new TH1D(uni_name, "", 20, 0, 3);
        universes[i]->SetLineColor(kAzure-4);

    }

    double reco_energy;
    std::map<std::string, std::vector<double>>* weights = nullptr;

    std::vector<std::string> tree_names = {
        "numu_sbnfit_cc0pinp",
        "numu_sbnfit_ccother",
        "numu_sbnfit_cosmic",
        "numu_sbnfit_mixed",
        "numu_sbnfit_oofv",
        "numu_sbnfit_nc",
        "numu_sbnfit_nue",
        "numu_sbnfit_anumu"
    };

    for (int i_tree = 0; i_tree < tree_names.size(); i_tree++){

        std::cout << tree_names.at(i_tree) << std::endl;

        TTree* t = (TTree*)_file0->Get(tree_names.at(i_tree).c_str());
        t->SetBranchAddress("reconstructed_neutrino_energy_calib", &reco_energy);
        t->SetBranchAddress("weights", &weights);

        std::cout << "N_entries: " << t->GetEntries() << std::endl;

        for (int i = 0; i < t->GetEntries(); i++){

            t->GetEntry(i);

            h_cv->Fill(reco_energy*1.25);

            std::vector<double> total_weight;
            total_weight.resize(5000, 1.);
            for (auto const& weight : *weights){
                if (weight.second.size() != 5000) continue;

                // loop universes
                for (int j = 0; j < 5000; j++){ 
                    total_weight.at(j) = total_weight.at(j)*weight.second.at(j);
                }
            }

            for (int i = 0; i < 5000; i++){

                universes[i]->Fill(reco_energy*1.25, total_weight.at(i));

            }
        }

    }

    h_cv->SetLineColor(kAzure-1);
    h_cv->SetLineWidth(2);
    h_cv->Scale(0.152283232);
    h_cv->SetMaximum(h_cv->GetMaximum()*1.6);
    h_cv->Draw();
    for (int i = 0; i < 5000; i++){
        universes[i]->Scale(0.152283232);
        universes[i]->Draw("samehisto");
    }
    h_cv->DrawClone("samehisto");

    TLegend *leg = new TLegend(0.3, 0.8, 0.83, 0.89);
    leg->AddEntry(h_cv, "Central Value Universe");
    leg->AddEntry(universes[0], "Alternate Universes");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw("same");

}
