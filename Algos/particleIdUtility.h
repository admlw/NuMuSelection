#ifndef PARTICLEIDUTILITY_H
#define PARTICLEIDUTILITY_H

// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// larsoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT includes
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

namespace pidutil{

  class particleIdUtility{

    public:

      using Point_t = ROOT::Math::PositionVector3D< ROOT::Math::Cartesian3D< Coord_t >>;

      static bool pairCompare(std::pair<int, int> firstEl, std::pair<int, int> secondEl) {
        bool isLarger = firstEl.first < secondEl.first;
        return isLarger;
      }

      std::pair<double, double> getAveragedQdX(recob::Track const& track, std::vector< art::Ptr< recob::Hit > > hitCollection, std::vector< art::Ptr< raw::RawDigit > > rawDCollection, bool isDQdXFromRawD);

      std::vector< std::vector< std::pair<double, double> > > getDeadRegions(std::vector< std::pair<int, int> > hitCollection);

      double getUntrackedLength(recob::Track const& track, std::vector< std::vector< std::pair<double,double> > > deadRegions);

      bool isMuonLike(double averagedQdX, double len);

  };

}

#endif
