/**
 * \author Adam Lister
 *
 * \email a.lister1@lancaster.ac.uk
 *
 * \description this framework is designed to work with the output of the 
 *              NuMuSelection module contained in the Module directory
 *              in this repository
 *              it pulls out many different data/mc comparisons at
 *              each stage of the selection
 *              the stages are currently set to be
 *                0. pure UBXSec CC-inclusive
 *                1. Topology cut: N Pfparticles
 *                2. Topology cut: N tracks
 *                3. Topology cut: N showers
 *                4. ParticleID cut
 */

#include "DoEnergyRecoStudies.h"

int main(){

  // pull out TTrees from provided TFiles
  // input trees are set in Configuration.h
  TFile* inputfile = new TFile("selectedEvents.root", "READ");
  TTree* t_simulation = (TTree*)inputfile->Get("simulation");

  // initialise variables and trees
  var_list simulation_vars_tmp;

  var_list* simulation_vars = &simulation_vars_tmp;

  _treeHandler.SetTreeVars(t_simulation, &simulation_vars_tmp, true);

  //------------------------------------
  // loop simulation
  //------------------------------------
 
  std::cout << "[DES] Beginning simulation loop..." << std::endl;
  for (int i = 0; i < t_simulation->GetEntries(); i++){

    t_simulation->GetEntry(i);

    for (int j = 0; j < simulation_vars->true_match_starte->size(); j++){

      if (simulation_vars->true_match_pdg->at(j) == 2212){
      


      }
    
    }

  }

  return 0;

}
