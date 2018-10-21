#ifndef ANALYSISCUTS_H
#define ANALYSISCUTS_H

// local
#include "DataTypes.h"
#include "Configuration.h"

// cpp
#include <iostream>

namespace numusel{

  class AnalysisCuts{

    public:

      int n_pfParticles_min_cut_val = 2;
      int n_pfParticles_max_cut_val = 100;
      int n_tracks_min_cut_val = 0;
      int n_tracks_max_cut_val = 100;
      int n_showers_min_cut_val = 0;
      int n_showers_max_cut_val = 100;
      float pid_cutvalue = -1.0;

      /**
       * Does event meet number of PFPs specification?
       * this is mostly redundant because of previous
       * cuts but is useful to tests
       */
      bool isPassNPfparticleSelection(var_list* vars);

      /**
       * Does event meet number of track specification?
       */
      bool isPassNTrackSelection(var_list* vars);

      /**
       * Does event meet number of shower specification?
       */
      bool isPassNShowerSelection(var_list* vars);

      /**
       * Does event pass PID cut?
       *
       * pair.first: bool which returns whether event passes
       * pair.second: vector of pdg codes returned from PID
       *
       * this assumes the default pid cut values
       */
      std::pair< bool, std::vector<int> > isPassParticleIDSelection(var_list *vars);

      /**
       *
       * same as above but assumes we pass new cut values
       */
      std::pair< bool, std::vector<int> > isPassParticleIDSelection(var_list *vars, float cutval);

  };

}

#endif
