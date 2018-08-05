#ifndef ENERGYRECO_H
#define ENERGYRECO_H

// cpp
#include <iostream>

// local
#include "DataTypes.h"
#include "Configuration.h"

namespace numusel{

  class EnergyReco{

    public:

      /**
       * Sets accessors for tree variables to be contained in the var_list
       */
      void EnergyByRange(var_list* vars);

  };

}
#endif
