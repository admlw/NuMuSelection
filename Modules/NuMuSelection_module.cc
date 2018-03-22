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
#include "lardata/Utilities/AssociationUtil.h"

// ROOT Includes
#include "TTree.h"

// UBXSec Includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

// EventWeight includes
#include "uboone/EventWeight/EventWeightTreeUtility.h"

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
  produces< std::vector<ubana::TPCObject> >();
  produces< art::Assns<ubana::TPCObject, ubana::SelectionResult> >();

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
 
  /**
   * Module produces a new ubana::SelectionResult, a new ubana::TPCObkect
   * and an association between the two
   */
  std::unique_ptr< std::vector<ubana::SelectionResult> > selectionCollection( new std::vector<ubana::SelectionResult> );
  std::unique_ptr< std::vector<ubana::TPCObject> > tpcObjectCollection( new std::vector<ubana::TPCObject> );
  std::unique_ptr< art::Assns<ubana::TPCObject, ubana::SelectionResult> > selTpcObjAssn( new art::Assns<ubana::TPCObject, ubana::SelectionResult>);

  // selection info
  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel("UBXSec", selectionHandle);
  if (!selectionHandle.isValid()){

    std::cout << "[NUMUSEL] SelectionResult product not found. " << std::endl;
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

  std::cout << "[NUMUSEL] PRINTING INFORMATION FOR EVENT " 
    << fRun << "." << fSubRun << "." << fEvent << std::endl;

  /**
   * If there is more than one selected object, then 
   * move to the next event.
   */
  if (selectedTpcObjects.at(0).size() == 1){ 

    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const std::vector<recob::Track>& selectedTracks = selectedTpcObject->GetTracks();
    const size_t nSelectedShowers = selectedTpcObject->GetNShowers();

    int muonLikeCounter = 0;

    /**
     * Loop tracks and check whether each trak is compatibile with
     * being a muon. If more than one track is muon-like then 
     * the event does not pass.
     *
     * This is currently a stop-gap fix until working PID can be put
     * here
     */
    for (size_t i = 0; i < selectedTracks.size(); i++){

      const recob::Track& track = selectedTracks.at(i);

      std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(track.ID());

      std::pair<double, double> averageDqdx = 
        pidutils.getAveragedQdX(track, hits, rawDVec, false);

      bool isMuLike = 
        pidutils.isMuonLike(averageDqdx.first, averageDqdx.second);

      if (isMuLike){ 
        muonLikeCounter++;
      }

    }

    std::cout << "[NUMUSEL] Found " 
      << muonLikeCounter << " muon candidates! " << std::endl;

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

    tpcObjectCollection->push_back(*(selectedTpcObject.get()));
  }
  else if (selectedTpcObjects.at(0).size() == 0){
    std::cout 
      << "[NUMUSEL] No TPC object selected, on to the next event!" << std::endl;
  }
  else {
    std::cout 
      << "[NUMUSEL] Uh-oh, there's more than one selected TPC object, that shouldn't be right..." << std::endl;
  }

  std::cout << "[NUMUSEL] PRINTING SELECTION INFORMATION" << std::endl;
  std::cout << "[NUMUSEL] SELECTION RESULT: " 
    << selResult.GetSelectionStatus() << std::endl;
  if (selResult.GetSelectionStatus() == 0){
    std::cout << "[NUMUSEL] REASON: " 
      << selResult.GetFailureReason() << "\n" << std::endl;
  }

  selectionCollection->push_back(selResult);
  util::CreateAssn(*this, 
      e, 
      *tpcObjectCollection, 
      *selectionCollection, 
      *selTpcObjAssn, 
      0, 
      selectionCollection->size());

  e.put(std::move(selectionCollection));
  e.put(std::move(tpcObjectCollection));
  e.put(std::move(selTpcObjAssn));
}

void NuMuSelection::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NuMuSelection)
