#ifndef CONFIGURATION_H
#define CONFIGURATION_H

// cpp
#include <iostream>

namespace numusel{

  class Configuration{

    public:

      std::string s_onbeam     = "/uboone/data/users/alister1/numuSelection/files/onbeam_selectionInformation.root";
      std::string s_offbeam    = "/uboone/data/users/alister1/numuSelection/files/offbeam_selectionInformation.root";
      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/BNBCOS_dic_devdataset_selectionInformation.root"; // dic
      std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/BNBCOS_nominal_devdataset_selectionInformation.root";
      std::string s_ew         = "/uboone/data/users/alister1/numuSelection/files/BNBCOS_nominal_devdataset_eventweight.root";

      double bnbcosPOT = 1.96673e+20; //nominal
      //double bnbcosPOT = 1.95899e+20; // dic
      double onbeam_tor860_wcut = 2.995e+19;
      double onbeam_E1DCNT_wcut = 7074111;
      double offbeam_EXT = 10199374;

      double offbeamscaling = onbeam_E1DCNT_wcut/offbeam_EXT;
      double simscaling = onbeam_tor860_wcut/bnbcosPOT;

      bool DoPIDForShowers = true;
      bool DoPIDForTracks = true;

      bool DoEventWeightMatching = false;

      // makes plots for track distributions separated by true pdg
      bool MakeTrackPlots = true;
      // this places a cut on the track variables, which is defined in 
      // SelectionMaker::PushBackTrackCutVar
      bool UseTrackCut = false;

      // means we only run over 1000 events for sim, onbeam and offbeam
      bool QuickDev = false;
      int QuickDevEvents = 100;

      // there are four stages to the selection
      // * no cuts (pure UBXSec)
      // * n tracks cut
      // * n showers cut
      // * n pfparticles cut
      // * pid cut
      int n_stages = 5; 

      // Energy reconstruction fit parameters
      // These are obtained by running 
      // ./DoSelection
      // ./DoEnergyRecoStudies
      // and then filling in these variables.
      // I know it's not clean.

      float proton_range_m = 1.00946;
      float proton_range_c = -0.00399997;
      float muon_range_contained_m = 0.999954;
      float muon_range_contained_c = 2.5331e-8;
      float muon_mcs_uncontained_m = 0.97469;
      float muon_mcs_uncontained_c = -0.016;

  };

}

#endif
