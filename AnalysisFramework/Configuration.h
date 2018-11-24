#ifndef CONFIGURATION_H
#define CONFIGURATION_H

// cpp
#include <iostream>

namespace numusel{

  class Configuration{

    public:
/*
      // vars
      std::string s_onbeam     = "/uboone/data/users/alister1/numuSelection/18-08-20-NuMuSelection/numusel_onbeam_info.root";
      std::string s_offbeam    = "/uboone/data/users/alister1/numuSelection/18-08-20-NuMuSelection/numusel_offbeam_info.root";
      std::string s_simulation = "/uboone/data/users/alister1/numuSelection/18-08-20-NuMuSelection/numusel_bnbcos_info.root";

      // bnbcos POT 1.4708e20
      // onbeam tor860_wcut 1.555e20
      // onbeam E1DCNT_wcut 34627337
      // offbeam EXT 33667849

      double offbeamscaling = 1.028;
      double simscaling = 1.06; 
*/

      /** testBuildShowersAsTracks_2
      std::string s_onbeam     = "/uboone/data/users/alister1/testBuildShowersAsTracks_2/selectionInformation_onbeam.root";
      std::string s_offbeam    = "/uboone/data/users/alister1/testBuildShowersAsTracks_2/selectionInformation_offbeam.root";
      std::string s_simulation = "/uboone/data/users/alister1/testBuildShowersAsTracks_2/selectionInformation_bnbcos.root";

      double bnbcosPOT = 1.9715e+20; //nominal
      //double bnbcosPOT = 9.34638e+19; // dic
      double onbeam_tor860_wcut = 3.245e+19;
      double onbeam_E1DCNT_wcut = 7244484;
      double offbeam_EXT = 10946246;
*/

      // testBuildShowersAsTracks_3
      std::string s_onbeam     = "/uboone/data/users/alister1/testBuildShowersAsTracks_3/selectionInformation_onbeam.root";
      std::string s_offbeam    = "/uboone/data/users/alister1/testBuildShowersAsTracks_3/selectionInformation_offbeam.root";
      std::string s_simulation = "/uboone/data/users/alister1/testBuildShowersAsTracks_4/selectionInformation_bnbcos.root";
      std::string s_ew = "/uboone/data/users/alister1/testBuildShowersAsTracks_3/ew_bnbcos.root";

      double bnbcosPOT = 1.97315e+20; //nominal
      double onbeam_tor860_wcut = 3.245e+19;
      double onbeam_E1DCNT_wcut = 7244484;
      double offbeam_EXT = 10946246;

      double offbeamscaling = onbeam_E1DCNT_wcut/offbeam_EXT;
      double simscaling = onbeam_tor860_wcut/bnbcosPOT;

      bool DoPIDForShowers = true;
      bool DoPIDForTracks = true;

      bool DoEventWeightMatching = false;

      // makes plots for track distributions separated by true pdg
      bool MakeTrackPlots = true;
      // this places a cut on the track variables, which is defined in 
      // SelectionMaker::PushBackTrackCutVar
      bool UseTrackCut = true;

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
