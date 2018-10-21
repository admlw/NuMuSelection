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

      // test building showers as tracks
      std::string s_onbeam     = "/uboone/data/users/alister1/testBuildShowersAsTracks/numusel_showersastracks_onbeam.root";
      std::string s_offbeam    = "/uboone/data/users/alister1/testBuildShowersAsTracks/numusel_showersastracks_offbeam.root";
      std::string s_simulation = "/uboone/data/users/alister1/testBuildShowersAsTracks/numusel_showersastracks_bnbcosmic.root";
      //std::string s_simulation = "/uboone/data/users/alister1/testBuildShowersAsTracks/numusel_showersastracks_bnbcosmic_dic3.root";

      double bnbcosPOT = 1.96e+20; //nominal
      //double bnbcosPOT = 9.34638e+19; // dic
      double onbeam_tor860_wcut = 3.225e+19;
      double onbeam_E1DCNT_wcut = 7199010;
      double offbeam_EXT = 10891089;

      double offbeamscaling = onbeam_E1DCNT_wcut/offbeam_EXT;
      double simscaling = onbeam_tor860_wcut/bnbcosPOT;

      bool DoPIDForShowers = true;
      bool DoPIDForTracks = true;

      // makes plots for track distributions separated by true pdg
      bool MakeTrackPlots = true;
      // this places a cut on the track variables, which is defined in 
      // SelectionMaker::PushBackTrackCutVar
      bool UseTrackCut = false;

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