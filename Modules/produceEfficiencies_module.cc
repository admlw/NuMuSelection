////////////////////////////////////////////////////////////////////////
// Class:       produceEfficiencies
// Plugin Type: analyzer (art v2_05_00)
// File:        produceEfficiencies_module.cc
//
// Generated at Fri Nov  3 16:49:57 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
//
// To be run on output of the CC Inclusive selection by Marco del Tutto
//
////////////////////////////////////////////////////////////////////////

// Base Includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art Includes
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft Includes
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT Includes
#include "TTree.h"
#include "TEfficiency.h"

// UBXSec Includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"

// FMWK Includes

class produceEfficiencies;


class produceEfficiencies : public art::EDAnalyzer {
  public:
    explicit produceEfficiencies(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    produceEfficiencies(produceEfficiencies const &) = delete;
    produceEfficiencies(produceEfficiencies &&) = delete;
    produceEfficiencies & operator = (produceEfficiencies const &) = delete;
    produceEfficiencies & operator = (produceEfficiencies &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:

    art::ServiceHandle< art::TFileService > tfs;
    TTree* selectionEfficiency;

    // vars
    bool isEventPassed = false;
    double mcNuEnergy;

    // Efficiency histograms
    TEfficiency* mcNuEnergyEff;
};


produceEfficiencies::produceEfficiencies(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{}

void produceEfficiencies::analyze(art::Event const & e)
{

  // get MC neutrino
  art::Handle< std::vector<simb::MCNeutrino> > mcNuHandle;
  e.getByLabel("generator", mcNuHandle);
  if (!mcNuHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCNeutrino> > mcNuVec;
  art::fill_ptr_vector(mcNuVec, mcNuHandle);

  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel("UBXSec", selectionHandle);
  if (!selectionHandle.isValid()){

    mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found."
      << std::endl;
    throw std::exception(); 

  }

  // get true neutrino energy
  if (mcNuVec.size() == 0) return;
  art::Ptr<simb::MCNeutrino> mcNu = mcNuVec.at(0);

  mcNuEnergy = mcNu->Nu().E();

  for (auto const& selectionStatus : (*selectionHandle)){

    if (selectionStatus.GetSelectionStatus()){

      // event passes
      isEventPassed = true;

      mcNuEnergyEff->Fill(isEventPassed, mcNuEnergy);

    }
    else {
    
      // event fails

    }

  }

}

void produceEfficiencies::beginJob()
{
  // Implementation of optional member function here.
  selectionEfficiency = tfs->make<TTree>("selectionEfficiency", "selectionEfficiency");
  mcNuEnergyEff = tfs->make<TEfficiency>("mcNuEnergyEff", ";#nu_{E}^{true}; #epsilon", 15, 0, 3);

}

void produceEfficiencies::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(produceEfficiencies)
