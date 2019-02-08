void produce_nue_backgrounds(){

    double reco_energy;
    std::map< std::string, std::vector<double> >* weights = nullptr;
    
    double reco_energy_out;
    std::map< std::string, std::vector<double> >* weights_out = nullptr;

    TFile* f_out = new TFile("nue_bnb_backgrounds.root", "RECREATE");

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

        TH1D* h = new TH1D("h", ";;", 5, 0, 3);
        
        TTree *t = (TTree*)_file0->Get(obj->GetName());
        t->SetBranchAddress("weights", &weights);

        TTree* t_out = new TTree(t->GetName(), "");
        t_out->Branch("reco_energy", &reco_energy_out);
        t_out->Branch("weights", &weights_out);

        TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
        c1->cd();

        if (std::string(obj->GetName()) != "cosmic"){
            t->Draw("reco_energy >> h");
        }
        
        h->Draw();

        TF1 *f1 = new TF1("f1","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1)", 0, 10);
        f1->SetParameters(1.0, 1.0, 1.0);
        h->Fit("f1");
      
        int nEvents = h->Integral()*26;
        float bins[12] = {0.2, 0.3, 0.375, 0.475, 0.550, 0.675, 0.800, 0.950, 1.10, 1.300, 1.500, 3.000};
        TH1D* h_new = new TH1D("h_new", ";;", 11, bins);

        t->GetEntry(0);
        for (int i = 0; i < nEvents; i++){

            reco_energy_out = f1->GetRandom(); 

            (*weights_out).clear();
            std::cout << (*weights_out).size() << std::endl;
            for ( auto const this_weight : *weights){
                std::cout << i << std::endl;
                std::vector<double> new_universes;
                for (int k = 0; k < this_weight.second.size(); k++){
                    new_universes.push_back(1.0);
                }
                (*weights_out)[std::string(this_weight.first)] = new_universes;
            }

            t_out->Fill();

            h_new->Fill(reco_energy_out);
        }

        
        h_new->SetLineColor(kBlack);
        c1->SaveAs((std::string(obj->GetName())+std::string(".pdf")).c_str());

        TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
        h_new->Draw();
        c2->SaveAs((std::string(obj->GetName())+std::string("_new.pdf")).c_str());
        
        h->Delete();
        h_new->Delete();

    }

    f_out->Write();    

}
