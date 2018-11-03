#include "SplitTree.h"

void SplitTree(){

  numusel::TreeHandler _treehandler;

  TFile* _file0 = new TFile("selectedEvents_sim.root", "READ");

  TTree* selection_tree = (TTree*)_file0->Get("simulation");
  TTree* ew_tree = (TTree*)_file0->Get("ew");

  var_list simulation_vars_tmp;
  var_list* simulation_vars = &simulation_vars_tmp;
  ew_list ew_vars_tmp;
  ew_list* ew_vars = &ew_vars_tmp;

  _treehandler.SetTreeVars(selection_tree, &simulation_vars_tmp, true);
  _treehandler.SetEWTreeVars(ew_tree, &ew_vars_tmp);

  int n_entries_in_file = 50;
  int n_files = std::ceil((double)selection_tree->GetEntries()/n_entries_in_file);

  std::cout << "[TreeSplitter] Number of entries: " << selection_tree->GetEntries() << std::endl;
  std::cout << "[TreeSplitter] Splitting trees into " << n_files << " files" << std::endl;

  for (int i_file = 0; i_file < n_files; i_file++){

    TString file_name = Form("split_tree/output_file_%i.root", i_file);
    TFile *file_out = new TFile(file_name, "RECREATE");
    file_out->cd();

    TTree* selection_tree_out = (TTree*)selection_tree->CloneTree(0);
    TTree* ew_tree_out = (TTree*)ew_tree->CloneTree(0);

    for (int i_entry = i_file*n_entries_in_file; i_entry < (i_file*n_entries_in_file) + n_entries_in_file; i_entry++){

      if (i_entry >= selection_tree->GetEntries()) continue;

      selection_tree->GetEntry(i_entry);
      selection_tree_out->Fill();

      ew_tree->GetEntry(i_entry);
      ew_tree_out->Fill();

    }

    selection_tree_out->Write();
    ew_tree_out->Write();

    file_out->Close();
  }
}

int main(){
  
  SplitTree();
  return 0;

}
