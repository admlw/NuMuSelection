#ifndef TREEHANDLER_H
#define TREEHANDLER_H

// cpp
#include <iostream>
#include <bitset>
#include <string>
#include <vector>

// ROOT
#include "TTree.h"
#include "TDirectory.h"

// local
#include "DataTypes.h"
#include "Configuration.h"

namespace numusel{

  class TreeHandler{

    public:

      /**
       * Set all branch status to 0 except for run, subrun, event, 
       * which are used as searcing criteria. This results in a 
       * significant speed up of the code
       */
      void PrepareTreeForSearching(TTree* tree);

      /**
       * Set all branch status to 1
       */
      void PrepareTreeForWriting(TTree* tree);

      /**
       * Sets accessors for tree variables to be contained in the var_list
       */
      void SetTreeVars(TTree* tree, var_list* varstoset, bool b_isSimulation);

      /**
       * Sets accessors for event weight tree
       */
      void SetEWTreeVars(TTree* tree, ew_list* varstoset);

      /**
       * Used to find the entry in the eventweight tree corresponding to the 
       * run, subrun, and event of the analysis tree
       */
      int FindEntryFromEvent(TTree* ewin, ew_list* ewvars, int run, int subrun, int event, int startentry);

  };

}
#endif
