#ifndef SELECTIONMAKER_H
#define SELECTIONMAKER_H

// local
#include "AnalysisCuts.h"
#include "DataTypes.h"
#include "HistogramHandler.h"

// cpp
#include <vector>
#include <iostream>
#include <algorithm>

// ROOT
#include "TTree.h"

namespace numusel{

  class SelectionMaker{

    public:
      std::vector<std::vector<std::vector<double>>> thisMatrix;

      std::vector<std::vector<double>> m_stage0;
      std::vector<std::vector<double>> m_stage1;
      std::vector<std::vector<double>> m_stage2;
      std::vector<std::vector<double>> m_stage3;
      std::vector<std::vector<double>> m_stage4;

      std::vector<std::vector<double>>* s_stage0 = &m_stage0;
      std::vector<std::vector<double>>* s_stage1 = &m_stage1;
      std::vector<std::vector<double>>* s_stage2 = &m_stage2;
      std::vector<std::vector<double>>* s_stage3 = &m_stage3;
      std::vector<std::vector<double>>* s_stage4 = &m_stage4;
 
      std::vector<std::vector<std::vector<std::pair<double,double>>>> thisMatrix_2d;

      std::vector<std::vector<std::pair<double,double>>> m_stage0_2d;
      std::vector<std::vector<std::pair<double,double>>> m_stage1_2d;
      std::vector<std::vector<std::pair<double,double>>> m_stage2_2d;
      std::vector<std::vector<std::pair<double,double>>> m_stage3_2d;
      std::vector<std::vector<std::pair<double,double>>> m_stage4_2d;

      std::vector<std::vector<std::pair<double,double>>>* s_stage0_2d = &m_stage0_2d;
      std::vector<std::vector<std::pair<double,double>>>* s_stage1_2d = &m_stage1_2d;
      std::vector<std::vector<std::pair<double,double>>>* s_stage2_2d = &m_stage2_2d;
      std::vector<std::vector<std::pair<double,double>>>* s_stage3_2d = &m_stage3_2d;
      std::vector<std::vector<std::pair<double,double>>>* s_stage4_2d = &m_stage4_2d;
     
      /**
       * gets variables to plot for each stage of the selection with input cut value
       */
      std::vector<std::vector<std::vector<double>>> GetPlottingVariables(var_list* vars, kVarType var_type, float cutvalue, TTree* infile=nullptr, TTree* outfile=nullptr, int entry=-1);
 
      /**
       * gets variables to plot for each stage of the selection
       */
      std::vector<std::vector<std::vector<double>>> GetPlottingVariables(var_list* vars, kVarType var_type, TTree* infile = nullptr, TTree* outfile = nullptr, int entry=-1);
 
      /**
       * gets variables to plot in 2d
       */
      std::vector<std::vector<std::vector<std::pair<double, double>>>> Get2DPlottingVariables(var_list* vars, kVarType var_type, float cutval, TTree* infile, TTree* outfile, int entry=-1);
 
      /**
       * gets variables to plot in 2d
       */
      std::vector<std::vector<std::vector<std::pair<double, double>>>> Get2DPlottingVariables(var_list* vars, kVarType var_type, TTree* infile, TTree* outfile, int entry=-1);
 
      /**
       * This function is used to place cuts on the individual tracks for testing 
       * purposes, but is only applied to some plots (i.e. track level), so be careful using it
       *
       * You can turn this feature on/off in the Configuration class
       */
      void PushBackTrackCutVar(std::vector<std::vector<double>>* vec, var_list* vars);

      /**
       * Fill Efficiency Vectors
       */
      void PushBackEPVectors(std::vector<std::vector<double>>* vec, var_list* vars);

      /**
       * Fill PID vectors, i.e. for each track what is the true PDG code of the matched MCParticle
       */
      void PushBackPIDVectors(std::vector<std::vector<double>>* vec, var_list* vars);

      /**
       * Fill 1D vectors
       */
      void PushBack1DVectors(std::vector<std::vector<double>>* vec, var_list* vars, bool isHasPID);

      /**
       * Fill 2D vectors
       */
      void PushBack2DVectors(std::vector<std::vector<std::pair<double, double>>>* vec, var_list* vars, bool isHasPID);

  };

}

#endif
