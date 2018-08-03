#ifndef TREEHANDLER_H
#define TREEHANDLER_H

#include <bitset>

// ROOT
#include "TTree.h"

// local
#include "DataTypes.h"
#include "Configuration.h"

namespace numusel{

  class TreeHandler{

    public:

      /**
       * Sets accessors for tree to be appropriate vars in the var struct
       */
      void SetTreeVars(TTree* tree, var_list* varstoset, bool b_isSimulation);
  
  };

}

#endif
