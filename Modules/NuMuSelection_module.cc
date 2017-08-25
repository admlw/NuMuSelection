////////////////////////////////////////////////////////////////////////
// Class:       NuMuSelection
// Plugin Type: analyzer (art v2_05_00)
// File:        NuMuSelection_module.cc
//
// Generated at Thu Aug 24 21:55:04 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class NuMuSelection;


class NuMuSelection : public art::EDAnalyzer {
public:
  explicit NuMuSelection(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuMuSelection(NuMuSelection const &) = delete;
  NuMuSelection(NuMuSelection &&) = delete;
  NuMuSelection & operator = (NuMuSelection const &) = delete;
  NuMuSelection & operator = (NuMuSelection &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

};


NuMuSelection::NuMuSelection(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void NuMuSelection::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

void NuMuSelection::beginJob()
{
  // Implementation of optional member function here.
}

void NuMuSelection::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

void NuMuSelection::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NuMuSelection)
