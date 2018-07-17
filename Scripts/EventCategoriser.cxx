#ifndef EVENTCATEGORISER_CXX
#define EVENTCATEGORISER_CXX

#include "EventCategoriser.h"

namespace numusel{

  std::bitset<8> EventCategoriser::CategoriseEvent(var_list* vars){

    std::bitset<8> thisBitSet;

    // cosmic
    if (vars->isCosmic)
      thisBitSet.set(7);

    // mixed
    if (vars->isMixed)
      thisBitSet.set(6);

    // oofv
    if (!vars->isInFV)
      thisBitSet.set(5);

    // nc
    if (vars->true_nu_ccnc == 1)
      thisBitSet.set(4);

    // antinue or nue
    if (std::abs(vars->true_mcp_pdg->at(0)) == 12)
      thisBitSet.set(3);

    // antinumu
    if (vars->true_mcp_pdg->at(0) == -14)
      thisBitSet.set(2);

    bool isCC0Pi = EventCategoriser::IsCC0PiEvent(vars->true_mcp_pdg, vars->true_mcp_process);

    // cc other
    if (!isCC0Pi)
      thisBitSet.set(1);

    // cc0pi
    if (isCC0Pi)
      thisBitSet.set(0);

    return thisBitSet;

  }

  bool EventCategoriser::IsCC0PiEvent(std::vector<double>* pdgs, std::vector<std::string>* processes){

    int nPrimaryPions = 0;

    for (int i = 0; i < (int)pdgs->size(); i++){

      if ((std::abs(pdgs->at(i)) == 211 || std::abs(pdgs->at(i)) == 111) && processes->at(i) == "primary"){
        nPrimaryPions++;

      }
    }

    if (nPrimaryPions == 0)
      return true;
    else return false;

  }

}

#endif
