#ifndef SELECTIONMAKER_H
#define SELECTIONMAKER_H

// local
#include "AnalysisCuts.h"
#include "DataTypes.h"

//cpp
#include <vector>
#include <iostream>

namespace numusel{

  class SelectionMaker{

    public:
      std::vector<std::vector<double>> thisMatrix;

      std::vector<double> m_stage0;
      std::vector<double> m_stage1;
      std::vector<double> m_stage2;
      std::vector<double> m_stage3;

      /**
       * gets variables to plot for each stage of the selection
       */
      std::vector<std::vector<double>> GetPlottingVariables(var_list* vars);

  };

}

#endif
