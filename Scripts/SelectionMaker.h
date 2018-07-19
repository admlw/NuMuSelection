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
      std::vector<std::vector<std::vector<double>>> thisMatrix;

      std::vector<std::vector<double>> m_stage0;
      std::vector<std::vector<double>> m_stage1;
      std::vector<std::vector<double>> m_stage2;
      std::vector<std::vector<double>> m_stage3;

      std::vector<std::vector<double>>* s_stage0 = &m_stage0;
      std::vector<std::vector<double>>* s_stage1 = &m_stage1;
      std::vector<std::vector<double>>* s_stage2 = &m_stage2;
      std::vector<std::vector<double>>* s_stage3 = &m_stage3;

      /**
       * gets variables to plot for each stage of the selection with input cut value
       */
      std::vector<std::vector<std::vector<double>>> GetPlottingVariables(var_list* vars, bool isEffPur, float cutvalue);
 
      /**
       * gets variables to plot for each stage of the selection
       */
      std::vector<std::vector<std::vector<double>>> GetPlottingVariables(var_list* vars, bool isEffPur);
      
      void PushBackEPVectors(std::vector<std::vector<double>>* vec, var_list* vars);

      void PushBackVVectors(std::vector<std::vector<double>>* vec, var_list* vars, bool isHasPID);

  };

}

#endif
