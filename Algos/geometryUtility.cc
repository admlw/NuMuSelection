#include "geometryUtility.h"

namespace xsecutils{

  class geometryUtils{

    std::vector<float> getTpcDimensions(){

      detWidth = geo.DetHalfWidth() * 2.0;
      detHeight = geo.DetHalfHeight() * 2.0;
      detLength = geo.DetLength() * 2.0;

      std::vector<float> tpcDimensions = {
        detWidth,
        detHeight,
        detLength
      };

      return tpcDimensions;

    }


  }

}
