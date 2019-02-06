#ifndef CONFIGURATION_H
#define CONFIGURATION_H

// cpp
#include <iostream>

namespace numusel{

  class Configuration{

    public:

      std::string s_onbeam     = "/uboone/data/users/alister1/numuSelection/files/190125/onbeam_selectionInformation.root";
      double onbeam_tor860_wcut = 4.89e+19;
      double onbeam_E1DCNT_wcut = 10905211;

      std::string s_offbeam    = "/uboone/data/users/alister1/numuSelection/files/190125/offbeam_selectionInformation.root";
      double offbeam_EXT = 77219137;
   
      std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_bnbcos.root";
      double bnbcosPOT = 1.85731e+21; 
      std::string var = "bnbcos";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_cv.root";
      //double bnbcosPOT = 1.93572e+20;
      //std::string var = "cv";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_sce.root";
      //double bnbcosPOT = 3.92382e+20;
      //std::string var = "sce";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_larg4bugfix.root";
      //double bnbcosPOT = 1.97068e+20;
      //std::string var = "larg4bugfix";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dlup.root";
      //double bnbcosPOT = 1.94757e+20;
      //std::string var = "dlup";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dldown.root";
      //double bnbcosPOT = 1.96824e+20;
      //std::string var = "dldown";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dtup.root";
      //double bnbcosPOT = 1.9537e+20;
      //std::string var = "dtup";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dtdown.root";
      //double bnbcosPOT = 1.97528e+20;
      //std::string var = "dtdown";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_noiseampup.root";
      //double bnbcosPOT = 1.94358e+20;
      //std::string var = "noiseampup";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_noiseampdown.root";
      //double bnbcosPOT = 1.95603e+20;
      //std::string var = "noiseampdown";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_penoiseup.root";
      //double bnbcosPOT = 1.95691e+20;
      //std::string var = "penoiseup";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_penoisedown.root";
      //double bnbcosPOT = 1.95691e+20;
      //std::string var = "penoisedown";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_dic.root";
      //double bnbcosPOT = 1.95567e+20;
      //std::string var = "dic"

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_squeezeresp.root";
      //double bnbcosPOT = 1.93515e+20;
      //std::string var = "squeezeresp";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_stretchresp.root";
      //double bnbcosPOT = 1.93901e+20;
      //std::string var = "stretchresp";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_deadsaturatechannels.root";
      //double bnbcosPOT = 1.97028e+20;
      //std::string var = "deadsaturatedchannels";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_altdeadchannels.root";
      //double bnbcosPOT = 1.93426e+20;
      //std::string var = "altdeadchannels";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_enhancedexttpcvis.root";
      //double bnbcosPOT = 1.96627e+20;
      //std::string var = "enhancedexttpcvis";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_lifetime10ms.root";
      //double bnbcosPOT = 1.96407e+20;
      //std::string var = "lifetime10ms";

      //std::string s_simulation = "/uboone/data/users/alister1/numuSelection/files/190125/selectionInformation_birks.root";
      //double bnbcosPOT = 1.97028e+20;
      //std::string var = "birks";

      std::string s_dirt       = "/uboone/data/users/alister1/numuSelection/files/190125/dirt_selectionInformation.root";
      double dirtPOT = 1.68796e+21;
      
      std::string s_ew         = "/uboone/data/users/alister1/numuSelection/files/190125/ew_bnbcos.root";

      std::string s_ew_dirt    = "/uboone/data/users/alister1/numuSelection/files/190125/ew_bnbcos.root";


      double offbeamscaling = onbeam_E1DCNT_wcut/offbeam_EXT;
      double simscaling = onbeam_tor860_wcut/bnbcosPOT;
      double dirtscaling = onbeam_tor860_wcut/dirtPOT;

      bool DoPIDForShowers = true;
      bool DoPIDForTracks = true;

      bool DoEventWeightMatching = true;

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

      float proton_range_m = 1.030;
      float proton_range_c = -0.007;
      float muon_range_contained_m = 1.00;
      float muon_range_contained_c = 0.00;
      float muon_mcs_uncontained_m = 0.983;
      float muon_mcs_uncontained_c = 0.002;

  };

}

#endif
