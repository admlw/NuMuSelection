#ifndef CONFIGURATION_H
#define CONFIGURATION_H

// cpp
#include <iostream>

namespace numusel{

  class Configuration{

    public:

      // vars
      std::string s_onbeam     = "/uboone/data/users/alister1/numuSelection/18-08-19-NuMuSelection/numusel_onbeam_info.root";
      std::string s_offbeam    = "/uboone/data/users/alister1/numuSelection/18-08-19-NuMuSelection/numusel_offbeam_info.root";
      std::string s_simulation = "/uboone/data/users/alister1/numuSelection/18-08-19-NuMuSelection/numusel_bnbcos_info.root";

      // bnbcos POT 1.87429e20
      // onbeam tor860_wcut 1.555e20
      // onbeam E1DCNT_wcut 34627337
      // offbeam EXT 33583135

      double offbeamscaling = 1.03;
      double simscaling = 0.83; 

      // there are four stages to the selection
      // * no cuts (pure UBXSec)
      // * n tracks cut
      // * n showers cut
      // * n pfparticles cut
      // * pid cut
      int n_stages = 5; 

  };

}

#endif
