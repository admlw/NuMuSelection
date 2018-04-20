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
#include "canvas/Persistency/Common/FindMany.h"

// LArSoft Includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"

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

    /**
     * initialises branches in output tree
     */
    void initialiseTree(TTree *t);

    /**
     * clear vectors
     */
    void clearVectors();

    /**
     * resize vectors
     */
    void resizeVectors();

    /**
     * check whether event is truly CC0Pi
     */
    bool is0PiEvent(art::Ptr<simb::MCTruth> mcTruth);

    /**
     * is the event truly have 2 tracks above threshold
     */
    bool isGT2TrackAboveThresholdEvent(art::Ptr<simb::MCTruth> mcTruth, double fProtonEThreshold);

    /**
     * is the event truly a 0-shower event
     */
    bool isZeroShowerEvent(art::Ptr<simb::MCTruth> mcTruth, double fElectronEThreshold);

    /**
     * check whether event has N protons above threshold
     */
    bool isNPEvent(art::Ptr<simb::MCTruth> mcTruth, double fProtonEThreshold, int fMinNProtons);

    /**
     * find mcp mached to a given track
     */
    simb::MCParticle const* GetAssociatedMCParticle(std::vector< art::Ptr<recob::Hit> > hits, art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit); 

  private:

    art::ServiceHandle< art::TFileService > tfs;
    TTree * tree;

    // intialise classes
    pidutil::particleIdUtility pidutils;
    ubana::SelectionResult selResult;
    ubana::FiducialVolume fiducialVolume;
    uboone::EWTreeUtil ewutil;
    art::ServiceHandle< geo::Geometry > geo;

    // fcl pars
    bool fIsData;
    bool isMc;
    int fMaxNShowers;
    int fMinNTracks;
    int fMinNProtons;
    double fProtonEThreshold;
    double fElectronEThreshold;
    std::string fSelectionLabel;
    std::string fPIDLabel;
    std::string fTrackLabel;
    std::string fHitLabel;
    std::string fSelectionTPCObjAssn;
    std::string fHitTrackAssn;
    std::string fHitTruthAssns;

    // signal bools
    bool isSelected;
    bool isCC;
    bool isNuMu;
    bool isNuMuBar;
    bool isNuE;
    bool isNuEBar;
    bool isTrueVtxInFV;
    bool isSignal;
    bool isBeamNeutrino;
    bool isMixed;
    bool isCosmic;
    bool isSelectedSignal;
    bool isGTNProtonsAboveThreshold;

    // bools for efficiency calculations
    bool isUBXSecSelected;
    bool isUBXSecSignal; //<! NuMu CC interaction inside the TPC

    /** does the selected TPC object have more than two tracks? */
    bool isTwoTrackSelected;
    bool isTwoTrackTruth;

    /** does the selected TPC object have zeo showers? */
    bool isZeroShowerSelected;
    bool isZeroShowerTruth;

    /** does the selected TPC object have any pions? */
    bool is0PiSelected;
    bool is0PiTruth;

    // vars
    int run;
    int subrun;
    int event;

    // truth
    double mcNuVx;
    double mcNuVy;
    double mcNuVz;
    int mcNuCCNC;
    int mcNuMode;
    int mcNuInteractionType;
    double mcq3Transfer;
    double mcq0Transfer;
    double mcNuEnergy;
    double mcNuMom;
    double mcNuPx;
    double mcNuPy;
    double mcNuPz;
    double mcNuTheta;
    double mcLeptonEnergy;
    double mcLeptonMom;
    double mcLeptonPx;
    double mcLeptonPy;
    double mcLeptonPz;
    double mcLeptonTheta; //<! Theta from lepton initial momentum
    double mcLeptonCosTheta; //<! Cosine theta from lepton intial momentum
    double mcLeptonPhi; //<! Phi from lepton initial momentum

    /**
     * MCParticle information of MCP which is matched to muon candidate
     */
    int mcpMuonCandTrackId;
    int mcpMuonCandStatusCode;
    int mcpMuonCandPdgCode;
    int mcpMuonCandMother;
    std::string mcpMuonCandProcess;
    std::string mcpMuonCandEndProcess;
    int mcpMuonCandNDaughters;
    std::vector<int> mcpMuonCandDaughter;
    unsigned int mcpMuonCandNumTrajPoint;
    double mcpMuonCandVx;
    double mcpMuonCandVy;
    double mcpMuonCandVz;
    double mcpMuonCandEndX;
    double mcpMuonCandEndY;
    double mcpMuonCandEndZ;
    double mcpMuonCandT;
    double mcpMuonCandPx;
    double mcpMuonCandPy;
    double mcpMuonCandPz;
    double mcpMuonCandE;
    double mcpMuonCandP;
    double mcpMuonCandPT;
    double mcpMuonCandMass;
    double mcpMuonCandEndPx;
    double mcpMuonCandEndPy;
    double mcpMuonCandEndPz;
    double mcpMuonCandEndE;
    double mcpMuonCandTheta; //<! Theta from lepton initial momentum
    double mcpMuonCandCosTheta; //<! Cosine theta from lepton intial momentum
    double mcpMuonCandPhi; //<! Phi from lepton initial momentum

    /**
     * Save all MCParticle information for all neutrino induced particle
     */
    std::vector<int> mcpTrackIds;
    std::vector<int> mcpStatusCodes;
    std::vector<int> mcpPdgCodes;
    std::vector<int> mcpMothers;
    std::vector<std::string> mcpProcesses;
    std::vector<std::string> mcpEndProcesses;
    std::vector<int> mcpNDaughters;
    //std::vector<int> mcpDaughters;
    std::vector<unsigned int> mcpNumTrajPoints;
    std::vector<double> mcpVxs;
    std::vector<double> mcpVys;
    std::vector<double> mcpVzs;
    std::vector<double> mcpEndXs;
    std::vector<double> mcpEndYs;
    std::vector<double> mcpEndZs;
    std::vector<double> mcpTs;
    std::vector<double> mcpPxs;
    std::vector<double> mcpPys;
    std::vector<double> mcpPzs;
    std::vector<double> mcpEs;
    std::vector<double> mcpPs;
    std::vector<double> mcpPTs;
    std::vector<double> mcpMasses;
    std::vector<double> mcpEndPxs;
    std::vector<double> mcpEndPys;
    std::vector<double> mcpEndPzs;
    std::vector<double> mcpEndEs;

    // reco
    double nTracks;
    double muonCandLength;
    double muonCandTheta;
    double muonCandCosTheta;
    double muonCandPhi;
    double muonCandEndX;
    double muonCandEndY;
    double muonCandEndZ;
    double muonCandStartX;
    double muonCandStartY;
    double muonCandStartZ;
    double vertexX;
    double vertexY;
    double vertexZ;
    std::vector<double> tmptrackLength;
    std::vector<double> tmptrackTheta;
    std::vector<double> tmptrackCosTheta;
    std::vector<double> tmptrackPhi;
    std::vector<double> tmptrackEndX;
    std::vector<double> tmptrackEndY;
    std::vector<double> tmptrackEndZ;
    std::vector<double> tmptrackStartX;
    std::vector<double> tmptrackStartY;
    std::vector<double> tmptrackStartZ;
    std::vector<int> tmptrackMatchedMcpID;
    std::vector<double> trackLength;
    std::vector<double> trackTheta;
    std::vector<double> trackCosTheta;
    std::vector<double> trackPhi;
    std::vector<double> trackEndX;
    std::vector<double> trackEndY;
    std::vector<double> trackEndZ;
    std::vector<double> trackStartX;
    std::vector<double> trackStartY;
    std::vector<double> trackStartZ;
    std::vector<int>    trackMatchedMcpID;

};


NuMuSelection::NuMuSelection(fhicl::ParameterSet const & p)
  //:
  //  EDProducer(p)  // ,
  // More initializers here.
{

  fhicl::ParameterSet const p_fv = p.get<fhicl::ParameterSet>("FiducialVolumeSettings");

  fIsData      = p.get<bool> ("IsData");
  fMaxNShowers = p.get<int>("MaximumNShowers", 0);
  fMinNTracks  = p.get<int>("MinimumNTracks", 2);
  fMinNProtons = p.get<int>("MinimumNProtons", 1);
  fProtonEThreshold   = p.get<double>("ProtonEThreshold", 0.04);
  fElectronEThreshold = p.get<double>("ElectronEThreshold", 0.02);
  
  fSelectionLabel = p.get<std::string>("SelectionLabel", "UBXSec");
  fPIDLabel   = p.get<std::string>("ParticleIdLabel", "particleid");
  fTrackLabel = p.get<std::string>("TrackLabel", "pandoraNu::UBXSec");
  fHitLabel   = p.get<std::string>("HitLabel", "pandoraCosmicHitRemoval::UBXSec");

  fSelectionTPCObjAssn = p.get<std::string>("SelectionTPCObjAssn", "UBXSec");
  fHitTrackAssn = p.get<std::string>("HitTrackAssn", "pandoraNu::UBXSec");
  fHitTruthAssns = p.get<std::string>("HitTruthAssn", "pandoraCosmicHitRemoval::UBXSec");

  produces< std::vector<ubana::SelectionResult> >();
  produces< std::vector<ubana::TPCObject> >();
  produces< art::Assns<ubana::TPCObject, ubana::SelectionResult> >();

  // configure
  fiducialVolume.Configure(p_fv,
      geo->DetHalfHeight(),
      2.*geo->DetHalfWidth(),
      geo->DetLength());

  fiducialVolume.PrintConfig();

}

void NuMuSelection::beginJob()
{

  tree = tfs->make<TTree>("tree", "tree");
  initialiseTree(tree);

}

void NuMuSelection::produce(art::Event & e)
{

  run = e.run();
  subrun = e.subRun();
  event = e.event();

  fIsData = e.isRealData();
  isMc = !fIsData;

  selResult.SetSelectionStatus(false);
  selResult.SetFailureReason("NoUBXSec");

  clearVectors();

  /**
   * Module produces a new ubana::SelectionResult, a new ubana::TPCObject
   * and an association between the two
   */
  std::unique_ptr< std::vector<ubana::SelectionResult> > selectionCollection( new std::vector<ubana::SelectionResult> );
  std::unique_ptr< std::vector<ubana::TPCObject> > tpcObjectCollection( new std::vector<ubana::TPCObject> );
  std::unique_ptr< art::Assns<ubana::TPCObject, ubana::SelectionResult> > selTpcObjAssn( new art::Assns<ubana::TPCObject, ubana::SelectionResult>);

  // selection info
  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel(fSelectionLabel, selectionHandle);
  if (!selectionHandle.isValid()){

    std::cout << "\n[NUMUSEL] SelectionResult product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found."
      << std::endl;
    throw std::exception(); 

  }

  // get Selected NuMuCC Inclusive TPCObject from selection handle
  art::FindManyP<ubana::TPCObject> selectedTpcObjects(selectionHandle, e, fSelectionTPCObjAssn);
  art::Ptr<ubana::TPCObject> selectedTpcObject; 

  // track handle
  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);

  // hit handle
  art::Handle< std::vector< recob::Hit > > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  // hit information
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, fHitTrackAssn);

  // hit<->mcParticle associations
  art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit(hitHandle, e, fHitTruthAssns);

//  art::FindManyP<anab::ParticleID> pidFromTrack(trackHandle, e, fPIDLabel);

  /**
   * Get MCTruth, MCFlux, and GTruth information for reweighting
   */

  art::Handle< std::vector< simb::MCTruth > > mcTruthHandle;
  e.getByLabel("generator", mcTruthHandle);
  if (!mcTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
  if (mcTruthVec.size() == 0){
    std::cout << "\n[NUMUSEL] No MCTruth Information" << std::endl;
    return;
  }

  art::Handle< std::vector< simb::GTruth > > gTruthHandle;
  e.getByLabel("generator", gTruthHandle);
  if (!gTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::GTruth> > gTruthVec;
  art::fill_ptr_vector(gTruthVec, gTruthHandle);
  if (gTruthVec.size() == 0){
    std::cout << "\n[NUMUSEL] No GTruth Information" << std::endl;
    return;
  }

  art::Handle< std::vector< simb::MCFlux > > mcFluxHandle;
  e.getByLabel("generator", mcFluxHandle);
  if (!mcFluxHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCFlux> > mcFluxVec;
  art::fill_ptr_vector(mcFluxVec, mcFluxHandle);
  if (mcFluxVec.size() == 0){
    std::cout << "\n[NUMUSEL] No MCFlux Information" << std::endl;
    return;
  }

  const art::Ptr<simb::MCFlux> mcFlux = mcFluxVec.at(0);
  const art::Ptr<simb::GTruth> gTruth = gTruthVec.at(0);

  const art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);
  std::cout << "\n[NUMUSEL] PRINTING INFORMATION FOR EVENT " 
    << run << "." << subrun << "." << event << std::endl;

  /**
   * Selection 
   *
   * Demand:
   * -- Only one selected TPC object
   * -- Only one track is considered muon-like
   * -- No showers reconstructed in the shower object
   */

  recob::Track muonCandidate;
  simb::MCParticle const* mcpMuonCandidate = 0;
  recob::Vertex selectedVertex;
  size_t nSelectedShowers;
  size_t nSelectedTracks=0;

  if (selectedTpcObjects.at(0).size() == 1){ 

    isUBXSecSelected = true;

    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const std::vector<recob::Track>& selectedTracks = selectedTpcObject->GetTracks();
    selectedVertex = selectedTpcObject->GetVertex();
    nSelectedShowers = selectedTpcObject->GetNShowers();
    nSelectedTracks  = selectedTpcObject->GetNTracks();
    const ubana::TPCObjectOrigin& selectedOrigin = selectedTpcObject->GetOrigin();

    if (isMc){
      // Check origin of TPC object

      isBeamNeutrino = false;
      isMixed = false;
      isCosmic = false;

      if (selectedOrigin == ubana::kBeamNeutrino) isBeamNeutrino = true;
      else isBeamNeutrino = false;

      if (selectedOrigin == ubana::kMixed) isMixed = true;
      else isMixed = false;

      if (selectedOrigin == ubana::kCosmicRay) isCosmic = true;
      else isCosmic = false;
    }

    // assume true to begin with
    selResult.SetSelectionStatus(true);

    /*
     * Topology cuts
     */
    if ((int)selectedTracks.size() < fMinNTracks){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("nTracksLTMinNTracks");
      isTwoTrackSelected = false;
    }
    else isTwoTrackSelected = true;

    if ((int)nSelectedShowers > fMaxNShowers){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("nShowersGTMaxNShowers");
      isZeroShowerTruth = false;
    }
    else isZeroShowerTruth = true;

    /**
     * Loop tracks and check whether each track is compatibile with
     * being a muon. If more than one track is muon-like then 
     * the event does not pass.
     *
     * This is currently a stop-gap fix until working PID can be put
     * here
     */

    if (selResult.GetSelectionStatus() == true){

      int muonLikeCounter = 0;

      for (size_t i = 0; i < selectedTracks.size(); i++){

        const recob::Track& track = selectedTracks.at(i);

        std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(track.ID());

        // get associated MCParticle
        simb::MCParticle const* mcpMatchedParticle = GetAssociatedMCParticle(hits, particlesPerHit);

        tmptrackLength.push_back(track.Length());
        tmptrackTheta.push_back(track.Theta());
        tmptrackCosTheta.push_back(std::cos(track.Theta()));
        tmptrackPhi.push_back(track.Phi());
        tmptrackEndX.push_back(track.End().X());
        tmptrackEndY.push_back(track.End().Y());
        tmptrackEndZ.push_back(track.End().Z());
        tmptrackStartX.push_back(track.Start().X());
        tmptrackStartY.push_back(track.Start().Y());
        tmptrackStartZ.push_back(track.Start().Z());
        tmptrackMatchedMcpID.push_back(mcpMatchedParticle->TrackId());

        std::pair<double, double> averageDqdx = 
          pidutils.getAveragedQdX(track, hits);

        bool isMuLike = 
          pidutils.isMuonLike(averageDqdx.first, averageDqdx.second);

        if (isMuLike){ 
          muonLikeCounter++;
          muonCandidate = track;
          mcpMuonCandidate = mcpMatchedParticle;
        }

      }

      std::cout << "[NUMUSEL] >> Found " 
        << muonLikeCounter << " muon candidates! " << std::endl;

      if (muonLikeCounter == 0){
        selResult.SetSelectionStatus(false);
        selResult.SetFailureReason("NoMu");
        is0PiSelected = false;
      }
      else if (muonLikeCounter > 1){
        selResult.SetSelectionStatus(false);
        selResult.SetFailureReason("PionCandidate");
        is0PiSelected = false;
      }
      else if (muonLikeCounter == 1){
        selResult.SetSelectionStatus(true);
        selResult.SetFailureReason("Passed");
        is0PiSelected = true;
      }
    }
    tpcObjectCollection->push_back(*(selectedTpcObject.get()));


  }
  else if (selectedTpcObjects.at(0).size() == 0){
    isUBXSecSelected = false;
    isSelected = false;
    isTwoTrackSelected = false;
    isZeroShowerSelected = false;
    if (isMc){
      isBeamNeutrino = false;
      isMixed = false;
      isCosmic = false;
    }
  }
  else {
    std::cout 
      << "[NUMUSEL] >> Uh-oh, there's more than one selected TPC object, that shouldn't be right..." << std::endl;
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

  std::cout << "[NUMUSEL] >> SELECTION RESULT: " 
    << selResult.GetSelectionStatus() << std::endl;
  if (selResult.GetSelectionStatus() == false){
    std::cout << "[NUMUSEL] >> REASON: " 
      << selResult.GetFailureReason() << "\n" << std::endl;
  }

  if (selResult.GetSelectionStatus() == true){
    isSelected = true;
  }
  else isSelected = false;

  /**
   * Set variables to write to tree
   */
  if (isSelected){
    trackLength       = tmptrackLength;
    trackTheta        = tmptrackTheta;
    trackCosTheta     = tmptrackCosTheta;
    trackPhi          = tmptrackPhi;
    trackEndX         = tmptrackEndX;
    trackEndY         = tmptrackEndY;
    trackEndZ         = tmptrackEndZ;
    trackStartX       = tmptrackStartX;
    trackStartY       = tmptrackStartY;
    trackStartZ       = tmptrackStartZ;
    trackMatchedMcpID = tmptrackMatchedMcpID;
    muonCandLength    = muonCandidate.Length();
    muonCandTheta     = muonCandidate.Theta();
    muonCandCosTheta  = std::cos(muonCandidate.Theta());
    muonCandPhi       = muonCandidate.Phi();
    muonCandStartX    = muonCandidate.Vertex().X();
    muonCandStartY    = muonCandidate.Vertex().Y();
    muonCandStartZ    = muonCandidate.Vertex().Z();
    muonCandEndX      = muonCandidate.End().X();
    muonCandEndY      = muonCandidate.End().Y();
    muonCandEndZ      = muonCandidate.End().Z();
    nTracks           = (int)nSelectedTracks;
    double xyz[3];
    selectedVertex.XYZ(xyz);
    vertexX                  = xyz[0];
    vertexY                  = xyz[1];
    vertexZ                  = xyz[2];
  }
  else {
    trackLength.push_back(-999); 
    trackTheta.push_back(-999);
    trackCosTheta.push_back(-999);
    trackPhi.push_back(-999);
    trackEndX.push_back(-999);         
    trackEndY.push_back(-999);         
    trackEndZ.push_back(-999);         
    trackStartX.push_back(-999);       
    trackStartY.push_back(-999);       
    trackStartZ.push_back(-999);       
    trackMatchedMcpID.push_back(-999);
    muonCandLength           = -999;
    muonCandTheta            = -999;
    muonCandCosTheta         = -999;
    muonCandPhi              = -999;
    muonCandStartX           = -999;
    muonCandStartY           = -999;
    muonCandStartZ           = -999;
    muonCandEndX             = -999;
    muonCandEndY             = -999;
    muonCandEndZ             = -999;
    nTracks                  = -999;
    vertexX                  = -999; 
    vertexY                  = -999;
    vertexZ                  = -999;
  }

  /**
   * If simulated data, then get true neutrino interaction information
   * for efficiency calculations
   */
  if (isMc){
    const simb::MCNeutrino& mcNu      = mcTruth->GetNeutrino();
    const simb::MCParticle& mcNuP     = mcNu.Nu();
    const simb::MCParticle& mcLeptonP = mcNu.Lepton();


    mcNuVx           = (double) mcNuP.Vx();
    mcNuVy           = (double) mcNuP.Vy();
    mcNuVz           = (double) mcNuP.Vz();
    mcNuMom          = mcNuP.P();
    mcNuPx           = mcNuP.Px();
    mcNuPy           = mcNuP.Py();
    mcNuPz           = mcNuP.Pz();
    mcNuEnergy       = mcNuP.E();
    mcLeptonMom      = mcLeptonP.P();
    mcLeptonPx       = mcLeptonP.Px();
    mcLeptonPy       = mcLeptonP.Py();
    mcLeptonPz       = mcLeptonP.Pz();
    mcLeptonEnergy   = mcLeptonP.E();
    mcLeptonTheta    = mcLeptonP.Momentum().Theta();
    mcLeptonCosTheta = std::cos(mcLeptonP.Momentum().Theta());
    mcLeptonPhi      = mcLeptonP.Momentum().Phi();

    for (int i = 0; i < mcTruth->NParticles(); i++){

      const simb::MCParticle& mcParticle = mcTruth->GetParticle(i);

      mcpTrackIds.push_back(mcParticle.TrackId());
      mcpStatusCodes.push_back(mcParticle.StatusCode());
      mcpPdgCodes.push_back(mcParticle.PdgCode());
      mcpMothers.push_back(mcParticle.Mother());
      mcpProcesses.push_back(mcParticle.Process());
      mcpEndProcesses.push_back(mcParticle.EndProcess());
      mcpNDaughters.push_back(mcParticle.NumberDaughters());
      //mcpDaughters.push_back();
      mcpNumTrajPoints.push_back(mcParticle.NumberTrajectoryPoints());
      mcpVxs.push_back(mcParticle.Vx());
      mcpVys.push_back(mcParticle.Vy());
      mcpVzs.push_back(mcParticle.Vz());
      mcpEndXs.push_back(mcParticle.EndX());
      mcpEndYs.push_back(mcParticle.EndY());
      mcpEndZs.push_back(mcParticle.EndZ());
      mcpTs.push_back(mcParticle.T());
      mcpPxs.push_back(mcParticle.Px());
      mcpPys.push_back(mcParticle.Py());
      mcpPzs.push_back(mcParticle.Pz());
      mcpEs.push_back(mcParticle.E());
      mcpPs.push_back(mcParticle.P());
      mcpPTs.push_back(mcParticle.Pt());
      mcpMasses.push_back(mcParticle.Mass());
      mcpEndPxs.push_back(mcParticle.EndPx());
      mcpEndPys.push_back(mcParticle.EndPy());
      mcpEndPzs.push_back(mcParticle.EndPz());
      mcpEndEs.push_back(mcParticle.EndE());

    }

    if (isSelected){
      mcpMuonCandTrackId      = mcpMuonCandidate->TrackId();
      mcpMuonCandStatusCode   = mcpMuonCandidate->StatusCode();
      mcpMuonCandPdgCode      = mcpMuonCandidate->PdgCode();
      mcpMuonCandMother       = mcpMuonCandidate->Mother();
      mcpMuonCandProcess      = mcpMuonCandidate->Process();
      mcpMuonCandEndProcess   = mcpMuonCandidate->EndProcess();
      mcpMuonCandNDaughters   = mcpMuonCandidate->NumberDaughters();
      for (int i = 0; i < mcpMuonCandNDaughters; i++){
        mcpMuonCandDaughter.push_back(mcpMuonCandidate->Daughter(i));
      }
      mcpMuonCandNumTrajPoint = mcpMuonCandidate->NumberTrajectoryPoints();
      mcpMuonCandVx           = mcpMuonCandidate->Vx();
      mcpMuonCandVy           = mcpMuonCandidate->Vy();
      mcpMuonCandVz           = mcpMuonCandidate->Vz();
      mcpMuonCandEndX         = mcpMuonCandidate->EndX();
      mcpMuonCandEndY         = mcpMuonCandidate->EndY();
      mcpMuonCandEndZ         = mcpMuonCandidate->EndZ();
      mcpMuonCandT            = mcpMuonCandidate->T();
      mcpMuonCandPx           = mcpMuonCandidate->Px();
      mcpMuonCandPy           = mcpMuonCandidate->Py();
      mcpMuonCandPz           = mcpMuonCandidate->Pz();
      mcpMuonCandE            = mcpMuonCandidate->E();
      mcpMuonCandP            = mcpMuonCandidate->P();
      mcpMuonCandPT           = mcpMuonCandidate->Pt();
      mcpMuonCandMass         = mcpMuonCandidate->Mass();
      mcpMuonCandEndPx        = mcpMuonCandidate->EndPx();
      mcpMuonCandEndPy        = mcpMuonCandidate->EndPy();
      mcpMuonCandEndPz        = mcpMuonCandidate->EndPz();
      mcpMuonCandEndE         = mcpMuonCandidate->EndE();
      mcpMuonCandTheta        = mcpMuonCandidate->Momentum().Theta();
      mcpMuonCandCosTheta     = std::cos(mcpMuonCandidate->Momentum().Theta());
      mcpMuonCandPhi          = mcpMuonCandidate->Momentum().Phi();
    }
    else{
      mcpMuonCandTrackId      = -999;
      mcpMuonCandStatusCode   = -999;
      mcpMuonCandPdgCode      = -999;
      mcpMuonCandMother       = -999;
      mcpMuonCandProcess      = "NotSelected";
      mcpMuonCandEndProcess   = "NotSelected";
      mcpMuonCandNDaughters   = -999;
      mcpMuonCandDaughter.push_back(-999);
      mcpMuonCandNumTrajPoint = -999;
      mcpMuonCandVx           = -999;
      mcpMuonCandVy           = -999;
      mcpMuonCandVz           = -999;
      mcpMuonCandT            = -999;
      mcpMuonCandPx           = -999;
      mcpMuonCandPy           = -999;
      mcpMuonCandPz           = -999;
      mcpMuonCandE            = -999;
      mcpMuonCandP            = -999;
      mcpMuonCandPT           = -999;
      mcpMuonCandMass         = -999;
      mcpMuonCandEndPx        = -999;
      mcpMuonCandEndPy        = -999;
      mcpMuonCandEndPz        = -999;
      mcpMuonCandEndE         = -999;
      mcpMuonCandTheta        = -999;
      mcpMuonCandCosTheta     = -999;
      mcpMuonCandPhi          = -999;

    }

    mcq0Transfer = mcNuEnergy - mcLeptonEnergy;
    mcq3Transfer = 
      std::sqrt(std::pow(mcNuPx-mcLeptonPx,2) +
          std::sqrt(std::pow(mcNuPy - mcLeptonPy, 2)) +
          std::sqrt(std::pow(mcNuPz - mcLeptonPz, 2)));

    mcNuCCNC = mcTruth->GetNeutrino().CCNC();
    mcNuMode            = mcNu.Mode();
    mcNuInteractionType = mcNu.InteractionType();

    /**
     * Get information on neutrino interaction
     */
    // is true nu CC?
    if (mcNuCCNC == 0) isCC = true;
    else isCC = false;

    // is true nu_mu?
    if (mcNuP.PdgCode() == 14) isNuMu = true;
    else isNuMu = false;

    // is true nu_mubar?
    if (mcNuP.PdgCode() == -14) isNuMuBar = true;
    else isNuMuBar = false;

    // is true nu_e?
    if (mcNuP.PdgCode() == 12) isNuE = true;
    else isNuE = false;

    // is true nu_ebar?
    if (mcNuP.PdgCode() == -12) isNuEBar = true;
    else isNuEBar = false;

    // is >2 track event?
    if (isGT2TrackAboveThresholdEvent(mcTruth, fProtonEThreshold) == 1) isTwoTrackTruth = true;
    else isTwoTrackTruth = false;

    // is 0 shower event?
    if (isZeroShowerEvent(mcTruth, fElectronEThreshold) == 1) isZeroShowerTruth = true;
    else isZeroShowerTruth = false;

    // is true nu 0Pi?
    if (is0PiEvent(mcTruth)) is0PiTruth = true;
    else is0PiTruth = false;

    // is N protons above threshold greater than demand?
    if (isNPEvent(mcTruth, fProtonEThreshold, fMinNProtons) == 1) isGTNProtonsAboveThreshold = true;
    else isGTNProtonsAboveThreshold = false;

    // is true nu vertex in FV
    if (fiducialVolume.InFV(mcNuVx, mcNuVy, mcNuVz)) isTrueVtxInFV = true;
    else isTrueVtxInFV = false;

    if (isCC && isNuMu && isTrueVtxInFV) isUBXSecSignal = true;
    else isUBXSecSignal = false;

    if (isCC && is0PiTruth && isNuMu && isTrueVtxInFV && isGTNProtonsAboveThreshold) isSignal = true;
    else isSignal = false;

    if (isSignal && isBeamNeutrino && isSelected) isSelectedSignal = true;
    else isSelectedSignal = false;

    /**
     * Event is selected, put all information of interest in a tree.
     */
    if (isSelectedSignal){
      ewutil.WriteTree(e, mcFlux, mcTruth, gTruth);
    }

    std::cout << "[NUMUSEL] >> isSelected                  " << isSelected << std::endl;
    std::cout << "[NUMUSEL] >> isSignal:                   " << isSignal << std::endl;
    std::cout << "[NUMUSEL] >> isTwoTrackTruth:            " << isTwoTrackTruth << std::endl;
    std::cout << "[NUMUSEL] >> isZeroShowerTruth:          " << isZeroShowerTruth << std::endl;
    std::cout << "[NUMUSEL] >> isGTNProtonsAboveThreshold: " << isGTNProtonsAboveThreshold << std::endl;
    std::cout << "[NUMUSEL] >> isBeamNeutrino:             " << isBeamNeutrino << std::endl;
    std::cout << "[NUMUSEL] >> isCosmic:                   " << isCosmic << std::endl;
    std::cout << "[NUMUSEL] >> isMixed:                    " << isMixed << std::endl;
    std::cout << "[NUMUSEL] >> isOOFV:                     " << !isTrueVtxInFV << std::endl;
    std::cout << "[NUMUSEL] >> isNC:                       " << !isCC << std::endl;
    std::cout << "[NUMUSEL] >> isNuMuBar:                  " << isNuMuBar << std::endl;
    std::cout << "[NUMUSEL] >> isNuE:                      " << isNuE << std::endl;
    std::cout << "[NUMUSEL] >> isNuEBar:                   " << isNuEBar << std::endl;
    std::cout << "[NUMUSEL] >> is0PiTruth:                      " << !is0PiTruth << std::endl;

  }

  tree->Fill();


}

void NuMuSelection::endJob()
{
  // Implementation of optional member function here.
}

void NuMuSelection::initialiseTree(TTree *t)
{
  t->Branch("fhicl_MaxNShowers"          , &fMaxNShowers               );
  t->Branch("fhicl_MinNTracks"           , &fMinNTracks                );
  t->Branch("fhicl_MinNProtons"          , &fMinNProtons               );
  t->Branch("fhicl_ProtonEThreshold"     , &fProtonEThreshold          );
  t->Branch("fhicl_ElectronEThreshold"   , &fElectronEThreshold        );
  t->Branch("run"                        , &run                        );
  t->Branch("subrun"                     , &subrun                     );
  t->Branch("event"                      , &event                      );
  t->Branch("isSelected"                 , &isSelected                 );
  t->Branch("isCC"                       , &isCC                       );
  t->Branch("isNuMu"                     , &isNuMu                     );
  t->Branch("isNuMuBar"                  , &isNuMuBar                  );
  t->Branch("isNuE"                      , &isNuE                      );
  t->Branch("isNuEBar"                   , &isNuEBar                   );
  t->Branch("isTrueVtxInFV"              , &isTrueVtxInFV              );
  t->Branch("is0PiTruth"                 , &is0PiTruth                 );
  t->Branch("isGTNProtonsAboveThreshold" , &isGTNProtonsAboveThreshold );
  t->Branch("isUBXSecSelected"           , &isUBXSecSelected           );
  t->Branch("isUBXSecSignal"             , &isUBXSecSignal             );
  t->Branch("isTwoTrackSelected"         , &isTwoTrackSelected         );
  t->Branch("isTwoTrackTruth"            , &isTwoTrackTruth            );
  t->Branch("isZeroShowerSelected"       , &isZeroShowerSelected       );
  t->Branch("isZeroShowerTruth"          , &isZeroShowerTruth          );
  t->Branch("is0PiSelected"              , &is0PiSelected              );
  t->Branch("is0PiTruth"                 , &is0PiTruth                 );
  t->Branch("isSignal"                   , &isSignal                   );
  t->Branch("isBeamNeutrino"             , &isBeamNeutrino             );
  t->Branch("isMixed"                    , &isMixed                    );
  t->Branch("isCosmic"                   , &isCosmic                   );
  t->Branch("isSelectedSignal"           , &isSelectedSignal           );
  t->Branch("mcNuCCNC"                   , &mcNuCCNC                   );
  t->Branch("mcNuMode"                   , &mcNuMode                   );
  t->Branch("mcNuInteractionType"        , &mcNuInteractionType        );
  t->Branch("mcq0Transfer"               , &mcq0Transfer               );
  t->Branch("mcq3Transfer"               , &mcq3Transfer               );
  t->Branch("mcNuVx"                     , &mcNuVx                     );
  t->Branch("mcNuVy"                     , &mcNuVy                     );
  t->Branch("mcNuVz"                     , &mcNuVz                     );
  t->Branch("mcNuPx"                     , &mcNuPx                     );
  t->Branch("mcNuPy"                     , &mcNuPy                     );
  t->Branch("mcNuPz"                     , &mcNuPz                     );
  t->Branch("mcNuEnergy"                 , &mcNuEnergy                 );
  t->Branch("mcNuTheta"                  , &mcNuTheta                  );
  t->Branch("mcLeptonMom"                , &mcLeptonMom                );
  t->Branch("mcLeptonPx"                 , &mcLeptonPx                 );
  t->Branch("mcLeptonPy"                 , &mcLeptonPy                 );
  t->Branch("mcLeptonPz"                 , &mcLeptonPz                 );
  t->Branch("mcLeptonEnergy"             , &mcLeptonEnergy             );
  t->Branch("mcLeptonTheta"              , &mcLeptonTheta              );
  t->Branch("mcLeptonCosTheta"           , &mcLeptonCosTheta           );
  t->Branch("mcLeptonPhi"                , &mcLeptonPhi                );
  t->Branch("mcpMuonCandTrackId"         , &mcpMuonCandTrackId         );
  t->Branch("mcpMuonCandStatusCode"      , &mcpMuonCandStatusCode      );
  t->Branch("mcpMuonCandPdgCode"         , &mcpMuonCandPdgCode         );
  t->Branch("mcpMuonCandMother"          , &mcpMuonCandMother          );
  t->Branch("mcpMuonCandProcess"         , &mcpMuonCandProcess         );
  t->Branch("mcpMuonCandEndProcess"      , &mcpMuonCandEndProcess      );
  t->Branch("mcpMuonCandNDaughters"      , &mcpMuonCandNDaughters      );
  
  t->Branch("mcpMuonCandDaughter"        , "std::vector<int>"             , &mcpMuonCandDaughter);
  
  t->Branch("mcpMuonCandNumTrajPoint"    , &mcpMuonCandNumTrajPoint    );
  t->Branch("mcpMuonCandVx"              , &mcpMuonCandVx              );
  t->Branch("mcpMuonCandVy"              , &mcpMuonCandVy              );
  t->Branch("mcpMuonCandVz"              , &mcpMuonCandVz              );
  t->Branch("mcpMuonCandEndX"            , &mcpMuonCandEndX            );
  t->Branch("mcpMuonCandEndY"            , &mcpMuonCandEndY            );
  t->Branch("mcpMuonCandEndZ"            , &mcpMuonCandEndZ            );
  t->Branch("mcpMuonCandT"               , &mcpMuonCandT               );
  t->Branch("mcpMuonCandPx"              , &mcpMuonCandPx              );
  t->Branch("mcpMuonCandPy"              , &mcpMuonCandPy              );
  t->Branch("mcpMuonCandPz"              , &mcpMuonCandPz              );
  t->Branch("mcpMuonCandE"               , &mcpMuonCandE               );
  t->Branch("mcpMuonCandP"               , &mcpMuonCandP               );
  t->Branch("mcpMuonCandPT"              , &mcpMuonCandPT              );
  t->Branch("mcpMuonCandMass"            , &mcpMuonCandMass            );
  t->Branch("mcpMuonCandEndPx"           , &mcpMuonCandEndPx           );
  t->Branch("mcpMuonCandEndPy"           , &mcpMuonCandEndPy           );
  t->Branch("mcpMuonCandEndPz"           , &mcpMuonCandEndPz           );
  t->Branch("mcpMuonCandEndE"            , &mcpMuonCandEndE            );
  t->Branch("mcpMuonCandTheta"           , &mcpMuonCandTheta           );
  t->Branch("mcpMuonCandCosTheta"        , &mcpMuonCandCosTheta        );
  t->Branch("mcpMuonCandPhi"             , &mcpMuonCandPhi             );

  t->Branch("mcpTrackIds"                , "std::vector<int>"             , &mcpTrackIds      );
  t->Branch("mcpStatusCodes"             , "std::vector<int>"             , &mcpStatusCodes   );
  t->Branch("mcpPdgCodes"                , "std::vector<int>"             , &mcpPdgCodes      );
  t->Branch("mcpMothers"                 , "std::vector<int>"             , &mcpMothers       );
  t->Branch("mcpProcesses"               , "std::vector<std::string>"     , &mcpProcesses     );
  t->Branch("mcpEndProcesses"            , "std::vector<std::string>"     , &mcpEndProcesses  );
  t->Branch("mcpNDaughters"              , "std::vector<int>"             , &mcpNDaughters    );
  //t->Branch("mcpDaughters"               , "std::vector<int>"             , &mcpDaughters     );
  t->Branch("mcpNumTrajPoints"           , "std::vector<unsigned int>"    , &mcpNumTrajPoints );
  t->Branch("mcpVxs"                     , "std::vector<double>"          , &mcpVxs           );
  t->Branch("mcpVys"                     , "std::vector<double>"          , &mcpVys           );
  t->Branch("mcpVzs"                     , "std::vector<double>"          , &mcpVzs           );
  t->Branch("mcpEndXs"                   , "std::vector<double>"          , &mcpEndXs         );
  t->Branch("mcpEndYs"                   , "std::vector<double>"          , &mcpEndYs         );
  t->Branch("mcpEndZs"                   , "std::vector<double>"          , &mcpEndZs         );
  t->Branch("mcpTs"                      , "std::vector<double>"          , &mcpTs            );
  t->Branch("mcpPxs"                     , "std::vector<double>"          , &mcpPxs           );
  t->Branch("mcpPys"                     , "std::vector<double>"          , &mcpPys           );
  t->Branch("mcpPzs"                     , "std::vector<double>"          , &mcpPzs           );
  t->Branch("mcpEs"                      , "std::vector<double>"          , &mcpEs            );
  t->Branch("mcpPs"                      , "std::vector<double>"          , &mcpPs            );
  t->Branch("mcpPTs"                     , "std::vector<double>"          , &mcpPTs           );
  t->Branch("mcpMasses"                  , "std::vector<double>"          , &mcpMasses        );
  t->Branch("mcpEndPxs"                  , "std::vector<double>"          , &mcpEndPxs        );
  t->Branch("mcpEndPys"                  , "std::vector<double>"          , &mcpEndPys        );
  t->Branch("mcpEndPzs"                  , "std::vector<double>"          , &mcpEndPzs        );
  t->Branch("mcpEndEs"                   , "std::vector<double>"          , &mcpEndEs         );
  
  t->Branch("nTracks"                    , &nTracks                    );
  
  t->Branch("trackLength"                , "std::vector<double>"          , &trackLength       );
  t->Branch("trackTheta"                 , "std::vector<double>"          , &trackTheta        );
  t->Branch("trackCosTheta"              , "std::vector<double>"          , &trackCosTheta     );
  t->Branch("trackPhi"                   , "std::vector<double>"          , &trackPhi          );
  t->Branch("trackEndX"                  , "std::vector<double>"          , &trackEndX         );
  t->Branch("trackEndY"                  , "std::vector<double>"          , &trackEndY         );
  t->Branch("trackEndZ"                  , "std::vector<double>"          , &trackEndZ         );
  t->Branch("trackStartX"                , "std::vector<double>"          , &trackStartX       );
  t->Branch("trackStartY"                , "std::vector<double>"          , &trackStartY       );
  t->Branch("trackStartZ"                , "std::vector<double>"          , &trackStartZ       );
  t->Branch("trackMatchedMcpID"          , "std::vector<int>"             , &trackMatchedMcpID );
  t->Branch("muonCandLength"             , &muonCandLength             );
  t->Branch("muonCandTheta"              , &muonCandTheta              );
  t->Branch("muonCandCosTheta"           , &muonCandCosTheta           );
  t->Branch("muonCandPhi"                , &muonCandPhi                );
  t->Branch("muonCandStartX"             , &muonCandStartX             );
  t->Branch("muonCandStartY"             , &muonCandStartY             );
  t->Branch("muonCandStartZ"             , &muonCandStartZ             );
  t->Branch("muonCandEndX"               , &muonCandEndX               );
  t->Branch("muonCandEndY"               , &muonCandEndY               );
  t->Branch("muonCandEndZ"               , &muonCandEndZ               );
  t->Branch("vertexX"                    , &vertexX                    );
  t->Branch("vertexY"                    , &vertexY                    );
  t->Branch("vertexZ"                    , &vertexZ                    );
}

void NuMuSelection::clearVectors()
{
  tmptrackLength.clear(); 
  tmptrackTheta.clear();
  tmptrackCosTheta.clear();
  tmptrackPhi.clear();
  tmptrackEndX.clear();         
  tmptrackEndY.clear();         
  tmptrackEndZ.clear();         
  tmptrackStartX.clear();
  tmptrackStartY.clear();       
  tmptrackStartZ.clear();       
  tmptrackMatchedMcpID.clear();
  trackLength.clear(); 
  trackTheta.clear();
  trackCosTheta.clear();
  trackPhi.clear();
  trackEndX.clear();         
  trackEndY.clear();         
  trackEndZ.clear();         
  trackStartX.clear();
  trackStartY.clear();       
  trackStartZ.clear();       
  trackMatchedMcpID.clear();
  mcpTrackIds.clear();
  mcpStatusCodes.clear();
  mcpPdgCodes.clear();
  mcpMothers.clear();
  mcpProcesses.clear();
  mcpEndProcesses.clear();
  mcpNDaughters.clear();
  //mcpDaughters.clear();
  mcpNumTrajPoints.clear();
  mcpVxs.clear();
  mcpVys.clear();
  mcpVzs.clear();
  mcpEndXs.clear();
  mcpEndYs.clear();
  mcpEndZs.clear();
  mcpTs.clear();
  mcpPxs.clear();
  mcpPys.clear();
  mcpPzs.clear();
  mcpEs.clear();
  mcpPs.clear();
  mcpPTs.clear();
  mcpMasses.clear();
  mcpEndPxs.clear();
  mcpEndPys.clear();
  mcpEndPzs.clear();
  mcpEndEs.clear();
}


bool NuMuSelection::is0PiEvent(art::Ptr<simb::MCTruth> mcTruth)
{
  int nParticles = mcTruth->NParticles();
  for (int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);

    if (particle.Process() != "primary" || particle.StatusCode() != 1) continue;

    if (std::abs(particle.PdgCode()) == 211 || std::abs(particle.PdgCode()) == 111)
      return false;

  }

  return true;
}

bool NuMuSelection::isGT2TrackAboveThresholdEvent(art::Ptr<simb::MCTruth> mcTruth, double fProtonEThreshold){

  int nParticles = mcTruth->NParticles();
  int particleCounter = 0;

  for (int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);
    if (particle.Process() != "primary" || particle.StatusCode() !=1) continue;

    if (std::abs(particle.PdgCode()) == 13  || 
        std::abs(particle.PdgCode()) == 211 ||
        (std::abs(particle.PdgCode()) == 2212 && particle.E() > fProtonEThreshold))
      particleCounter++;

  }

  if (particleCounter >= 2) return true;
  else return false;

}

bool NuMuSelection::isZeroShowerEvent(art::Ptr<simb::MCTruth> mcTruth, double fElectronEThreshold){

  int nParticles = mcTruth->NParticles();
  int particleCounter = 0;

  for (int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);
    if (particle.Process() != "primary" || particle.StatusCode() !=1) continue;

    if (std::abs(particle.PdgCode()) == 22  || 
        (std::abs(particle.PdgCode()) == 11 && particle.E() > fElectronEThreshold))
      particleCounter++;

  }

  if (particleCounter == 0 ) return true;
  else return false;

}

bool NuMuSelection::isNPEvent(art::Ptr<simb::MCTruth> mcTruth, double fProtonEThreshold, int fMinNProtons)
{

  int nParticles = mcTruth->NParticles();
  int nProtonsAboveThreshold = 0;
  for (int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);
    if (particle.Process() != "primary" || particle.StatusCode() !=1) continue;

    if (std::abs(particle.PdgCode()) == 2212 && particle.E() > fProtonEThreshold){
      nProtonsAboveThreshold++; 
    }
  }

  if (nProtonsAboveThreshold >= fMinNProtons)
    return true;
  else return false;
}

simb::MCParticle const* NuMuSelection::GetAssociatedMCParticle(std::vector< art::Ptr< recob::Hit > > hits, art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit)
{

  std::unordered_map<int,double> trkide;
  double maxe=-1, tote=0;
  simb::MCParticle const* maxp_me = NULL; //pointer for the particle match we will calculate

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

  //loop only over our hits
  for(size_t i_h=0; i_h<hits.size(); ++i_h){

    particle_vec.clear(); match_vec.clear();
    particlesPerHit.get(hits[i_h].key(),particle_vec,match_vec);
    //the .key() gives us the index in the original collection

    //loop over particles
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
      tote += match_vec[i_p]->energy; //calculate total energy deposited
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
        maxe = trkide[ particle_vec[i_p]->TrackId() ];
        maxp_me = particle_vec[i_p];
      }
    }//end loop over particles per hit

  }

  return maxp_me;

}

DEFINE_ART_MODULE(NuMuSelection)
