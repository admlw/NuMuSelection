#ifndef PARTICLEIDUTILITY_H
#define PARTICLEIDUTILITY_H

// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// larsoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

namespace pidutil{

  class particleIdUtility{

    public:

      double getAveragedQdX(recob::Track const& track, std::vector< art::Ptr< recob::Hit > > hitCollection);

  };

}

#endif
