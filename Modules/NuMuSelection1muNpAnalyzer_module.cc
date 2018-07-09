////////////////////////////////////////////////////////////////////////
// Class:       NuMuSelection1muNpAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        NuMuSelection1muNpAnalyzer_module.cc
//
// Generated at Mon Jul  9 10:35:20 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// This is a rewriting of the NuMuSelection producer module as an analyzer 
// module. 
//
// With the introduction of the new PID algorithms it became clear that the 
// ability to change the PID cuts on the fly would be extremely useful without
// having the full weight of LArSoft.
//
// Hindsight is 20/20.

// Base includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ART includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"        
#include "lardataobj/RecoBase/Hit.h"           
#include "lardata/Utilities/AssociationUtil.h" 
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h" 

// ROOT includes
#include "TTree.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"        
#include "uboone/UBXSec/DataTypes/TPCObject.h"              
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"        

// EventWeight includes                                     
#include "uboone/EventWeight/EventWeightTreeUtility.h"      

// local Includes                                           
#include "uboone/NuMuSelection/Algos/particleIdUtility.h"   


class NuMuSelection1muNpAnalyzer;


class NuMuSelection1muNpAnalyzer : public art::EDAnalyzer {
  public:
    explicit NuMuSelection1muNpAnalyzer(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    NuMuSelection1muNpAnalyzer(NuMuSelection1muNpAnalyzer const &) = delete;
    NuMuSelection1muNpAnalyzer(NuMuSelection1muNpAnalyzer &&) = delete;
    NuMuSelection1muNpAnalyzer & operator = (NuMuSelection1muNpAnalyzer const &) = delete;
    NuMuSelection1muNpAnalyzer & operator = (NuMuSelection1muNpAnalyzer &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endSubRun(art::SubRun const & sr) override;

    // User Functions

    void initialiseAnalysisTree(TTree*, bool fIsData);

  private:

    // initialise services
    art::ServiceHandle< art::TFileService > tfs;

    // initialise classes
    ubana::SelectionResult selResult;
    ubana::FiducialVolume fiducialVolume;
    uboone::EWTreeUtil ewutil;

    // fhicl parameters
    std::string fSelectionLabel;
    std::string fHitLabel;
    std::string fTrackLabel;
    std::string fPIDLabel;
    std::string fSelectionTPCObjAssn;
    std::string fHitTrackAssn;
    std::string fHitTruthAssn;
    bool fIsData;

    // tree and tree variables
    TTree *ana_tree;
    int run;
    int subrun;
    int event;
    bool isData;
    bool isSimulation;
    bool isUBXSecSelected;
    bool isBeamNeutrino;
    bool isMixed;
    bool isCosmic;
    int nSelectedTracks;
    int nSelectedShowers;
};


NuMuSelection1muNpAnalyzer::NuMuSelection1muNpAnalyzer(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
{

  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");
  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolumeSettings");

  fSelectionLabel      = p_labels.get<std::string>("SelectionLabel", "UBXSec");
  fHitLabel            = p_labels.get<std::string>("HitLabel", "pandoraCosmicHitRemoval::UBXSec");
  fTrackLabel          = p_labels.get<std::string>("TrackLabel", "pandoraNu::UBXSec");
  fPIDLabel            = p_labels.get<std::string>("fPIDLabel");

  fSelectionTPCObjAssn = p_labels.get<std::string>("SelectionTPCObjAssn", "UBXSec");
  fHitTrackAssn        = p_labels.get<std::string>("hitTrackAssn", "pandoraNu::UBXSec");
  fHitTruthAssn        = p_labels.get<std::string>("hitTruthAssn", "pandoraCosmicHitRemoval::UBXSec");

  fIsData = p.get<bool>("IsData", "false");

}

void NuMuSelection1muNpAnalyzer::analyze(art::Event const & e)
{

  run    = e.run();
  subrun = e.subRun();
  event  = e.event();
  isData = e.isRealData();
  isSimulation = !isData;

  // get selection information for this event
  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel(fSelectionLabel, selectionHandle);

  if (!selectionHandle.isValid()){

    std::cout << "[NuMuSelection] ubana::SelectionResult product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "[NuMuSelection] ubana::SelectionResult product not found." << std::endl;
    throw std::exception();

  }

  art::FindManyP<ubana::TPCObject> selectedTPCObjects(selectionHandle, e, fSelectionTPCObjAssn);
  art::Ptr<ubana::TPCObject> selectedTPCObject;

  // other handles to information we're going to need
  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);

  art::Handle< std::vector< recob::Hit > > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  // the assns we need
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, fHitTrackAssn);
  art::FindManyP< anab::ParticleID > pidFromTrack(trackHandle, e, fPIDLabel);

                                                                                  
  // finally, get  MCTruth, MCFlux, and GTruth information for reweighting                      
  
  art::Handle< std::vector< simb::MCTruth > > mcTruthHandle;                          
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;                                  

  art::Handle< std::vector< simb::MCFlux > > mcFluxHandle;                            
  std::vector< art::Ptr<simb::MCFlux> > mcFluxVec;                                    

  art::Handle< std::vector< simb::GTruth > > gTruthHandle;                            
  std::vector< art::Ptr<simb::GTruth> > gTruthVec;                                    


  if (isSimulation){                                                                     
    e.getByLabel("generator", mcTruthHandle);                                         
    if (!mcTruthHandle.isValid()) return;                                             
    art::fill_ptr_vector(mcTruthVec, mcTruthHandle);                                  
    if (mcTruthVec.size() == 0){                                                      
      std::cout << "\n[NuMuSelection] No MCTruth Information" << std::endl;                 
      return;                                                                         
    }                                                                                 

    e.getByLabel("generator", gTruthHandle);                                          
    if (!gTruthHandle.isValid()) return;                                              
    art::fill_ptr_vector(gTruthVec, gTruthHandle);                                    
    if (gTruthVec.size() == 0){                                                       
      std::cout << "\n[NuMuSelection] No GTruth Information" << std::endl;                  
      return;                                                                         
    }                                                                                 

    e.getByLabel("generator", mcFluxHandle);                                          
    if (!mcFluxHandle.isValid()) return;                                              
    art::fill_ptr_vector(mcFluxVec, mcFluxHandle);                                    
    if (mcFluxVec.size() == 0){                                                       
      std::cout << "\n[NuMuSelection] No MCFlux Information" << std::endl;                  
      return;                                                                         
    }                                                                                 
  }                                                                                   

  art::Ptr<simb::MCFlux>  mcFlux;
  art::Ptr<simb::GTruth>  gTruth;
  art::Ptr<simb::MCTruth> mcTruth;

  if (isSimulation){
    mcFlux  = mcFluxVec.at(0);
    gTruth  = gTruthVec.at(0);
    mcTruth = mcTruthVec.at(0);
  }

  std::cout << "[NuMuSelection] --- run.subrun.event: " 
    << run << "." << subrun << "." << event << std::endl;

  // initialise variables which need scope
  recob::Vertex selectedVertex;
  nSelectedShowers = 0;
  nSelectedTracks = 0;

  int nSelectedTPCObjects = selectedTPCObjects.at(0).size();

  // we expect a single TPC object
  if (nSelectedTPCObjects == 1){
    std::cout << "[NuMuSelection] Event is UBXSec Selected." << std::endl;
    isUBXSecSelected = true;

    // get selected TPC object and associated information
    selectedTPCObject = selectedTPCObjects.at(0).at(0);
    const std::vector<recob::Track>& selectedTracks = selectedTPCObject->GetTracks();
    selectedVertex   = selectedTPCObject->GetVertex();
    nSelectedShowers = (int)selectedTPCObject->GetNShowers();
    nSelectedTracks  = (int)selectedTPCObject->GetNTracks();
    const ubana::TPCObjectOrigin& selectedOrigin = selectedTPCObject->GetOrigin();

    if (isSimulation){

      isBeamNeutrino = false;
      isMixed        = false;
      isCosmic       = false;

      if (selectedOrigin == ubana::kBeamNeutrino) isBeamNeutrino = true;
      else isBeamNeutrino = false;

      if (selectedOrigin == ubana::kMixed) isMixed = true;
      else isMixed = false;

      if (selectedOrigin == ubana::kCosmicRay) isCosmic = true;
      else isCosmic = false;

    }

    // loop tracks in TPC object and get information
    for (size_t i = 0; i < selectedTracks.size(); i++){

      const recob::Track& track = selectedTracks.at(i);

      std::cout << track.Length() << std::endl;

    }


  }
  // this is probably fine
  else if (nSelectedTPCObjects == 0){

    std::cout << "[NuMuSelection] Event has " << nSelectedTPCObjects << std::endl;
    isUBXSecSelected = false;
    return;

  }
  // this could indicate a problem
  else if (nSelectedTPCObjects > 1){

    std::cout << "[NuMuSelection] Event has " << nSelectedTPCObjects << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "[NuMuSelection] Event has " << nSelectedTPCObjects << std::endl;
    throw std::exception();
  
  }

  ana_tree->Fill();

}

void NuMuSelection1muNpAnalyzer::beginJob()
{

  ana_tree = tfs->make<TTree>("analysis_tree", "analysis tree");
  initialiseAnalysisTree(ana_tree, fIsData);

}

void NuMuSelection1muNpAnalyzer::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void NuMuSelection1muNpAnalyzer::initialiseAnalysisTree(TTree *tree, bool fIsData){

  tree->Branch("run"    , &run);
  tree->Branch("subrun" , &subrun);
  tree->Branch("event"  , &event);
  tree->Branch("isData" , &isData);
  tree->Branch("isSimulation", &isSimulation);
  tree->Branch("isUBXSecSelected", &isUBXSecSelected);

}

DEFINE_ART_MODULE(NuMuSelection1muNpAnalyzer)
