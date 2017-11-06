#ifndef GEOMETRYUTILITY_H
#define GEOMETRYUTILITY_H

#include "larcore/Geometry/Geometry.h"

namespace xsecutils{

  class geometryUtility{

    private:
      geo::Geometry geo;
      float detWidth;
      float detHeight;
      float detLength;

    public:
      std::vector<float> getTpcDimensions();

  };

}

#endif
