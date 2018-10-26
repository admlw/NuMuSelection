void test(){

  std::cout << _file0->GetName() << std::endl;

  TTree* tree = (TTree*)_file0->Get("numuselection/analysis_tree");

  for (int i = 0; i < tree->GetEntries(); i++){
    
    tree->GetEntry(i);

  }

}
