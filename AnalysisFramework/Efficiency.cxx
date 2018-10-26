#ifndef EFFICIENCY_CXX
#define EFFICIENCY_CXX

#include "Efficiency.h"

namespace numusel{

  void Efficiency::FillEfficiencyNumerator(eff_1d* eff_plots, std::vector<double> val, std::bitset<8> evcat, var_list* vars){

    // if is CC0PiNP in FV and selected
    if (evcat[0] == 1 && evcat[5] == 0 && vars->isUBXSecSelected == true){

      for (int i = 0; i < val.size(); i++){
        
        eff_plots->h_num->Fill(val.at(i));
   

      }
    }

  };

  void Efficiency::FillEfficiencyDenominator(eff_1d* eff_plots, std::vector<double> val, std::bitset<8> evcat, var_list* vars){

    // if is CC0PiNP in FV
    if (evcat[0] == 1 && evcat[5] == 0){

      for (int i = 0; i < val.size(); i++){

        eff_plots->h_denom->Fill(val.at(i));
      
      }
    }

  };


}

#endif
