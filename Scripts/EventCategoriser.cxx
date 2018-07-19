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

    bool isCC0Pi = EventCategoriser::IsCC0PiNPEvent(vars->true_mcp_pdg, vars->true_mcp_process, vars->true_nu_ccnc, vars->true_mcp_starte, vars->isBeamNeutrino);

    // cc other
    if (!isCC0Pi)
      thisBitSet.set(1);

    // cc0pi
    if (isCC0Pi)
      thisBitSet.set(0);

    return thisBitSet;

  }

  bool EventCategoriser::IsCC0PiNPEvent(std::vector<double>* pdgs, std::vector<std::string>* processes, int ccnc, std::vector<double>* starte, bool isBeamNeutrino){

    int nPrimaryPions = 0;
    int nPrimaryProtonsAboveThreshold = 0;

    for (int i = 0; i < (int)pdgs->size(); i++){

      if ((std::abs(pdgs->at(i)) == 211 || std::abs(pdgs->at(i)) == 111) && processes->at(i) == "primary"){
        nPrimaryPions++;

      }
      
      if ((std::abs(pdgs->at(i)) == 2212) && starte->at(i) > 0.04 && processes->at(i) == "primary")
        nPrimaryProtonsAboveThreshold++;

    }

    if (nPrimaryPions == 0 && ccnc == 0 && nPrimaryProtonsAboveThreshold > 0 )
      return true;
    else return false;

  }

}

#endif
