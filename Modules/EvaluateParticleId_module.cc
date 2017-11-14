////////////////////////////////////////////////////////////////////////
// Class:       EvaluateParticleId
// Plugin Type: analyzer (art v2_05_00)
// File:        EvaluateParticleId_module.cc
//
// Generated at Mon Nov 13 20:18:00 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// base includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// local includes
#include "uboone/NuMuSelection/Algos/particleIdUtility.h"

class EvaluateParticleId;


class EvaluateParticleId : public art::EDAnalyzer {
public:
  explicit EvaluateParticleId(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EvaluateParticleId(EvaluateParticleId const &) = delete;
  EvaluateParticleId(EvaluateParticleId &&) = delete;
  EvaluateParticleId & operator = (EvaluateParticleId const &) = delete;
  EvaluateParticleId & operator = (EvaluateParticleId &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  double averagedQdX;
  pidutil::particleIdUtility pidutils;

};


EvaluateParticleId::EvaluateParticleId(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void EvaluateParticleId::analyze(art::Event const & e)
{

  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel("pandoraNu", trackHandle);

  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, "pandoraNu");

  for (auto const& track : (*trackHandle)){

    // get associated hits
    std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(track.ID());

    averagedQdX = pidutils.getAveragedQdX(track, hits); 

  }

}

void EvaluateParticleId::beginJob()
{
  // Implementation of optional member function here.
}

void EvaluateParticleId::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(EvaluateParticleId)
