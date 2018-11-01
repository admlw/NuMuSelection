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

      /**
       * Sets accessors for event weight tree
       */
      void SetEWTreeVars(TTree* tree, ew_list* varstoset);

      int FindEntryFromEvent(TTree* ewin, ew_list* ewvars, int run, int subrun, int event);

  };

}
#endif
