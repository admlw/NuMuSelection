void fix_robertos_shit(){


    TTree* t = (TTree*)_file0->Get("nu_e_cc0pinp");

    double nu_E;
    double reco_energy;
    std::map< std::string, std::vector<double> >* weights = nullptr;

    double nu_E_out;
    double reco_energy_out;
    std::map< std::string, std::vector<double> >* weights_out = nullptr;

    std::string out_file_name(std::string(_file0->GetName()));
    out_file_name.insert(out_file_name.length()-5, "_fixed");
    std::cout << "Output file is: " << out_file_name << std::endl;
    TFile* file_out = new TFile(out_file_name.c_str(), "RECREATE");

    TList* list = _file0->GetListOfKeys();
    if (!list) {
        std::cout << "no keys found in file" << std::endl;
        exit(1);
    }

    TObject* obj;

    for (int i = 0; i < list->GetSize(); i++){

        TKey* key = (TKey*)list->At(i);
        obj = key->ReadObj();

        if ( !obj->InheritsFrom("TTree")){
            std::cout << "Object does not inherit from TTree" << std::endl;
            continue;
        }

        TTree *t = (TTree*)_file0->Get(obj->GetName());
        t->SetBranchAddress("nu_E", &nu_E);
        t->SetBranchAddress("reco_energy", &reco_energy);
        t->SetBranchAddress("weights", &weights);

        TTree* tree_out = new TTree(t->GetName(), "");
        tree_out->Branch("nu_E", &nu_E_out);
        tree_out->Branch("reco_energy", &reco_energy_out);
        tree_out->Branch("weights", &weights_out);

        for (int i = 0; i < t->GetEntries()/2; i++){

            t->GetEntry(i);

            nu_E_out = nu_E;
            weights_out = weights;

            t->GetEntry(i*2);

            reco_energy_out = reco_energy;


            tree_out->Fill();
        }
    }

    file_out->Write();

}
