void fix_robertos_shit(){

    TFile* _file0 = new TFile("/uboone/data/users/alister1/SBNFit_files/filtered_events_mc_nueintrinsic_sbnfit_skimmed.root");

    TTree* t = (TTree*)_file0->Get("nu_e_cc0pinp");

    double nu_E;
    double reco_energy;

    t->SetBranchAddress("nu_E", &nu_E);
    t->SetBranchAddress("reco_energy", &reco_energy);

    TH2D* h = new TH2D("h", ";true_nu_e;reco_nu_e", 50, 0, 3, 50, 0, 3);

    for (int i = 0; i < t->GetEntries()/2; i++){

        double this_nu_E;
        double this_reco_energy;

        t->GetEntry(i);
        
        this_nu_E = nu_E;

        t->GetEntry(i*2);

        this_reco_energy = reco_energy;

        h->Fill(this_reco_energy, this_nu_E);

    }

    h->Draw("colz");

}
