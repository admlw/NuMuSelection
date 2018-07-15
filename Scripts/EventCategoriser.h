#ifndef EVENTCATEGORISER_H
#define EVENTCATEGORISER_H

// cpp
#include <bitset>
#include <vector>
#include <string>

// local
#include "DataTypes.h"

namespace numusel{

  class EventCategoriser{

    public:

      std::bitset<8> CategoriseEvent(var_list vars);

      bool IsCC0PiEvent(std::vector<double>* pdgs, std::vector<double>* processes);

  };

}

#endif
