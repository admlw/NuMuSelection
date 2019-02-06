#include "SortTree.h"

void SortTree(){

    TFile* _file0 = new TFile("/uboone/data/users/alister1/numuSelection/files/190125/ew_dirt.root", "READ");

    TTree* _tree0 = (TTree*)_file0->Get("tree");

    TFile* f = new TFile("out.root", "recreate");
    TTree* t = (TTree*)_tree0->CloneTree(0);

    _tree0->BuildIndex("run","event");
    TTreeIndex* index = (TTreeIndex*)_tree0->GetTreeIndex();


    for (int i = 0; i < index->GetN(); i++){
        if (i%10000 == 0) std::cout << i << "/" << index->GetN() << std::endl;
        Long64_t local = _tree0->LoadTree(index->GetIndex()[i]); 
        _tree0->GetEntry(local); 
        t->Fill();
    }
/*
    TTree* t1 = (TTree*)t->CloneTree(0);
    t->BuildIndex("subrun");
    TTreeIndex* index1 = (TTreeIndex*)t->GetTreeIndex();

    for (int i = 0; i < index1->GetN(); i++){
        if (i%10000 == 0) std::cout << i << "/" << index1->GetN() << std::endl;
        Long64_t local = t->LoadTree(index1->GetIndex()[i]); 
        t->GetEntry(local); 
        t1->Fill();
    }

    TTree* t2 = (TTree*)t1->CloneTree(0);
    t1->BuildIndex("run");
    TTreeIndex* index2 = (TTreeIndex*)t1->GetTreeIndex();

    for (int i = 0; i < index2->GetN(); i++){
        if (i%10000 == 0) std::cout << i << "/" << index2->GetN() << std::endl;
        Long64_t local = t1->LoadTree(index2->GetIndex()[i]); 
        t1->GetEntry(local); 
        t2->Fill();
    }
*/
    std::cout << "about to write..." << std::endl;

    f->cd();
    t->Write();
}

int main(){

    SortTree();
    return 0;

}
