#ifndef TREEHANDLER_H
#define TREEHANDLER_H

// cpp
#include <iostream>
#include <bitset>
#include <string>
#include <vector>

// ROOT
#include "TTree.h"

// local
#include "DataTypes.h"
#include "Configuration.h"

namespace numusel{

  class TreeHandler{

    public:

      /**
       * Sets accessors for tree variables to be contained in the var_list
       */
      void SetTreeVars(TTree* tree, var_list* varstoset, bool b_isSimulation);

  };

}
#endif
