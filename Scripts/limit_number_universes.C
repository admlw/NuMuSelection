void limit_number_universes(){

    int pre_filter_total_universes = 5000;
    int post_filter_total_universes = 1000;

    double reconstructed_neutrino_energy_calib;
    std::string reco_energy_name("reco_energy");
    //std::string reco_energy_name("reconstructed_neutrino_energy_calib");
    std::map< std::string, std::vector<double> >* weights = nullptr; 

    double reconstructed_neutrino_energy_calib_out;
    std::map< std::string, std::vector<double> >* weights_out = nullptr;
    
    // output stuff
    std::string out_file_name(std::string(_file0->GetName()));
    out_file_name.insert(out_file_name.length()-5, "_skimmed");
    std::cout << "Output file is: " << out_file_name << std::endl;
    TFile* file_out = new TFile(out_file_name.c_str(), "RECREATE");

    TList* list = _file0->GetListOfKeys();
    if (!list) {
        std::cout << "No keys found in file" << std::endl;
        exit(1);
    }

    TObject* obj;
    // first loop all of the trees in the file
    for (int i = 0; i < list->GetSize(); i++){

        TKey* key = (TKey*)list->At(i);
        obj = key->ReadObj();

        if ( !obj->InheritsFrom("TTree")){
            std::cout << "Object not a tree!" << std::endl;
            continue;
        }

        // get tree
        TTree* tree = (TTree*)_file0->Get(obj->GetName());

        tree->SetBranchAddress(reco_energy_name.c_str(), &reconstructed_neutrino_energy_calib);
        tree->SetBranchAddress("weights", &weights);

        TTree* tree_out = new TTree(tree->GetName(), "");
        tree_out->Branch(reco_energy_name.c_str(), &reconstructed_neutrino_energy_calib_out);
        tree_out->Branch("weights", &weights_out);
      
        // loop entries of tree
        for (int j = 0; j < tree->GetEntries(); j++){

            tree->GetEntry(j);

            reconstructed_neutrino_energy_calib_out = reconstructed_neutrino_energy_calib;

            // loop map
            (*weights_out).clear();
            for ( auto const this_weight : *weights){
                if (this_weight.second.size() == pre_filter_total_universes){
                    std::vector<double> new_universes;
                    for (int k = 0; k < post_filter_total_universes; k++){
                        new_universes.push_back(this_weight.second.at(k));
                    }
                    (*weights_out)[std::string(this_weight.first)] = new_universes;
                }
                else{
                    (*weights_out)[std::string(this_weight.first)] = this_weight.second;
                }

            }

            tree_out->Fill();

        }
       

    }

   file_out->Write();

}
