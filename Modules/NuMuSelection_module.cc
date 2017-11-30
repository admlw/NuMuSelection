////////////////////////////////////////////////////////////////////////
// Class:       NuMuSelection
// Plugin Type: producer (art v2_05_00)
// File:        NuMuSelection_module.cc
//
// Generated at Thu Aug 24 21:55:04 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// Base Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art Includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft Includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"

// ROOT Includes
#include "TTree.h"

// UBXSec Includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

// local Includes
#include "uboone/NuMuSelection/Algos/particleIdUtility.h"

class NuMuSelection;


class NuMuSelection : public art::EDProducer {
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
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:

    pidutil::particleIdUtility pidutils;

    ubana::SelectionResult selResult;

    // fcl pars
    bool isData;
    bool isMc;

    int fRun;
    int fSubRun;
    int fEvent;

};


NuMuSelection::NuMuSelection(fhicl::ParameterSet const & p)
  //:
  //  EDProducer(p)  // ,
  // More initializers here.
{

  isData = p.get< bool > ("IsData");

  produces< std::vector<ubana::SelectionResult> >();
//  produces< art::Assns<ubana::SelectionResult, ubana::SelectionResult> >();

}

void NuMuSelection::beginJob()
{

}

void NuMuSelection::produce(art::Event & e)
{

  fRun = e.run();
  fSubRun = e.subRun();
  fEvent = e.event();

  isData = e.isRealData();
  isMc = !isData;

  selResult.SetSelectionStatus(false);
  selResult.SetFailureReason("NoUBXSec");
  std::unique_ptr< std::vector<ubana::SelectionResult> > selectionCollection( new std::vector<ubana::SelectionResult> );

  // selection info
  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel("UBXSec", selectionHandle);
  if (!selectionHandle.isValid()){

    std::cout << " >> SelectionResult product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found."
      << std::endl;
    throw std::exception(); 

  }

  // get Selected NuMuCC Inclusive TPCObject from selection handle
  art::FindManyP<ubana::TPCObject> selectedTpcObjects(selectionHandle, e, "UBXSec");
  art::Ptr<ubana::TPCObject> selectedTpcObject; 

  // track handle
  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel("pandoraNu", trackHandle);

  // hit information
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, "pandoraNu");

  // get RawDigits to pass to dQdX calculation
  art::Handle< std::vector< raw::RawDigit > > rawDHandle;
  e.getByLabel("wcNoiseFilter", rawDHandle);
  std::vector< art::Ptr< raw::RawDigit > > rawDVec;
  art::fill_ptr_vector(rawDVec, rawDHandle);

  std::cout << "\n>| PRINTING INFORMATION FOR EVENT " << fRun << "." << fSubRun << "." << fEvent << std::endl;

  if (isData)
    std::cout << ">|>| Running over PHYSICS DATA" << std::endl;
  else
    std::cout << ">|>| Running over SIMULATED DATA" << std::endl;

  if (selectedTpcObjects.at(0).size() == 1){ 

    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const std::vector<recob::Track>& selectedTracks = selectedTpcObject->GetTracks();
    const size_t nSelectedShowers = selectedTpcObject->GetNShowers();

    int muonLikeCounter = 0;

    for (size_t i = 0; i < selectedTracks.size(); i++){

      const recob::Track& track = selectedTracks.at(i);

      std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(track.ID());

      std::pair<double, double> averageDqdx = pidutils.getAveragedQdX(track, hits, rawDVec, false);

      bool isMuLike = pidutils.isMuonLike(averageDqdx.first, averageDqdx.second);

      if (isMuLike){ 
        muonLikeCounter++;
      }

    }

    std::cout << "\n>|>| Found " << muonLikeCounter << " muon candidates! " << std::endl;

    if (muonLikeCounter == 0){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("NoMu");
    }
    else if (muonLikeCounter > 1){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("PionCandidate");
    }
    else if (muonLikeCounter == 1){
      selResult.SetSelectionStatus(true);
      selResult.SetFailureReason("Passed");
    }


    if (nSelectedShowers != 0){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("RecoShower");
    }



  }
  else if (selectedTpcObjects.at(0).size() == 0){

    std::cout << ">|>| No TPC object selected, on to the next event!" << std::endl;

  }
  else {

    std::cout << ">|>| Uh-oh, there's more than one selected TPC object, that shouldn't be right..." << std::endl;

  }

  std::cout << "\n>|>| PRINTING SELECTION INFORMATION" << std::endl;
  std::cout << ">|>| SELECTION RESULT: " << selResult.GetSelectionStatus() << std::endl;
  if (selResult.GetSelectionStatus() == 0)
    std::cout << ">|>| REASON: " << selResult.GetFailureReason() << "\n" << std::endl;

  selectionCollection->push_back(selResult);

  e.put(std::move(selectionCollection));
}

void NuMuSelection::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NuMuSelection)
