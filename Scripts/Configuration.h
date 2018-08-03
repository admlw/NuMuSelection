#ifndef CONFIGURATION_H
#define CONFIGURATION_H

// cpp
#include <iostream>

namespace numusel{

  class Configuration{

    public:

      // vars
      std::string s_onbeam     = "/uboone/data/users/alister1/numuSelection/18-07-16-NuMuSelection/numusel_onbeam.root";
      std::string s_offbeam    = "/uboone/data/users/alister1/numuSelection/18-07-16-NuMuSelection/numusel_offbeam.root";
      std::string s_simulation = "/uboone/data/users/alister1/numuSelection/18-07-16-NuMuSelection/numusel_bnbcos2.root";

      double offbeamscaling = 1.05;
      double simscaling = 2.3; 

      // there are four stages to the selection
      // * no cuts (pure UBXSec)
      // * n tracks cut
      // * n showers cut
      // * pid cut
      int n_stages = 4; 

  };

}

#endif
