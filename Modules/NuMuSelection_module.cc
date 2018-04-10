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
     * check whether event is truly CC0Pi
     */
    bool is0PiEvent(art::Ptr<simb::MCTruth> mcTruth);

    /**
     * is the event truly have 2 tracks above threshold
     */
    bool isGT2TrackAboveThresholdEvent(art::Ptr<simb::MCTruth> mcTruth, double protonEThreshold);

    /**
     * is the event truly a 0-shower event
     */
    bool isZeroShowerEvent(art::Ptr<simb::MCTruth> mcTruth, double electronEThreshold);

    /**
     * check whether event has N protons above threshold
     */
    bool isNPEvent(art::Ptr<simb::MCTruth> mcTruth, double protonEThreshold, int minNProtons);

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
    bool isData;
    bool isMc;
    int maxNShowers;
    int minNTracks;
    int minNProtons;
    double protonEThreshold;
    double electronEThreshold;

    // signal bools
    bool isSelected;
    bool isCC;
    bool is0Pi;
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
    bool isPionRejectionSelected;
    bool PionRejectionTruth;

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
    static const int kMaxParticles = 50;
    int mcpMuonCandTrackId;
    int mcpMuonCandStatusCode;
    int mcpMuonCandPdgCode;
    int mcpMuonCandMother;
    std::string mcpMuonCandProcess;
    std::string mcpMuonCandEndProcess;
    int mcpMuonCandNDaughters;
    int mcpMuonCandDaughter[kMaxParticles];
    unsigned int mcpMuonCandNumTrajPoint;
    double mcpMuonCandVx;
    double mcpMuonCandVy;
    double mcpMuonCandVz;
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
    /*    std::vector<int> mcpTrackIds;
          std::vector<int> mcpStatusCodes;
          std::vector<int> mcpPdgCodes;
          std::vector<int> mcpMothers;
          std::vector<std::string> mcpProcesses;
          std::vector<std::string> mcpEndProcesses;
          std::vector<int> mcpNDaughters;
          std::vector<int> mcpDaughters;
          std::vector<unsigned int> mcpNumTrajPoints;
          std::vector<double> mcpVxs;
          std::vector<double> mcpVys;
          std::vector<double> mcpVzs;
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
          */
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
};


NuMuSelection::NuMuSelection(fhicl::ParameterSet const & p)
  //:
  //  EDProducer(p)  // ,
  // More initializers here.
{

  isData      = p.get< bool > ("IsData");
  maxNShowers = p.get< int > ("MaximumNShowers", 0);
  minNTracks  = p.get< int > ("MinimumNTracks", 2);
  minNProtons = p.get< int > ("MinimumNProtons", 1);
  protonEThreshold   = p.get< double > ("ProtonEThreshold", 0.04);
  electronEThreshold = p.get< double > ("electronEThreshold", 0.02);

  produces< std::vector<ubana::SelectionResult> >();
  produces< std::vector<ubana::TPCObject> >();
  produces< art::Assns<ubana::TPCObject, ubana::SelectionResult> >();

  // configure
  fiducialVolume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
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

    std::cout << "\n[NUMUSEL] SelectionResult product not found. " << std::endl;
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

  // hit handle
  art::Handle< std::vector< recob::Hit > > hitHandle;
  e.getByLabel("pandoraCosmicHitRemoval::McRecoStage2", hitHandle);

  // hit information
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, "pandoraNu");

  // hit<->mcParticle associations
  art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit(hitHandle, e, "crHitRemovalTruthMatch::McRecoStage2");

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
  const art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);
  const art::Ptr<simb::GTruth> gTruth = gTruthVec.at(0);


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

    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const std::vector<recob::Track>& selectedTracks = selectedTpcObject->GetTracks();
    selectedVertex = selectedTpcObject->GetVertex();
    nSelectedShowers = selectedTpcObject->GetNShowers();
    nSelectedTracks  = selectedTpcObject->GetNTracks();
    const ubana::TPCObjectOrigin& selectedOrigin = selectedTpcObject->GetOrigin();

    if (isMc){
      // Check origin of TPC object
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
    if ((int)selectedTracks.size() < minNTracks){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("nTracksLTMinNTracks");
    }

    if ((int)nSelectedShowers > maxNShowers){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("nShowersGTMaxNShowers");
    }

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

        std::pair<double, double> averageDqdx = 
          pidutils.getAveragedQdX(track, hits);

        bool isMuLike = 
          pidutils.isMuonLike(averageDqdx.first, averageDqdx.second);

        if (isMuLike){ 
          muonLikeCounter++;
          muonCandidate = track;
          mcpMuonCandidate = GetAssociatedMCParticle(hits, particlesPerHit);
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
    }
    tpcObjectCollection->push_back(*(selectedTpcObject.get()));

  }
  else if (selectedTpcObjects.at(0).size() == 0){
    isSelected = false;
  }
  else {
    std::cout 
      << "[NUMUSEL] Uh-oh, there's more than one selected TPC object, that shouldn't be right..." << std::endl;
  }

  selectionCollection->push_back(selResult);
  util::CreateAssn(*this, 
      e, 
      *tpcObjectCollection, 
      *selectionCollection, 
      *selTpcObjAssn, 
      0, 
      selectionCollection->size());

  std::cout << "[NUMUSEL] SELECTION RESULT: " 
    << selResult.GetSelectionStatus() << std::endl;
  if (selResult.GetSelectionStatus() == false){
    std::cout << "[NUMUSEL] REASON: " 
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

    muonCandLength           = muonCandidate.Length();
    muonCandTheta            = muonCandidate.Theta();
    muonCandCosTheta         = std::cos(muonCandidate.Theta());
    muonCandPhi              = muonCandidate.Phi();
    muonCandStartX           = muonCandidate.Vertex().X();
    muonCandStartY           = muonCandidate.Vertex().Y();
    muonCandStartZ           = muonCandidate.Vertex().Z();
    muonCandEndX             = muonCandidate.End().X();
    muonCandEndY             = muonCandidate.End().Y();
    muonCandEndZ             = muonCandidate.End().Z();
    nTracks                  = (int)nSelectedTracks;
    double xyz[3];
    selectedVertex.XYZ(xyz);
    vertexX                  = xyz[0];
    vertexY                  = xyz[1];
    vertexZ                  = xyz[2];
  }
  else {
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

  e.put(std::move(selectionCollection));
  e.put(std::move(tpcObjectCollection));
  e.put(std::move(selTpcObjAssn));

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

    mcpMuonCandTrackId      = mcpMuonCandidate->TrackId();
    mcpMuonCandStatusCode   = mcpMuonCandidate->StatusCode();
    mcpMuonCandPdgCode      = mcpMuonCandidate->PdgCode();
    mcpMuonCandMother       = mcpMuonCandidate->Mother();
    mcpMuonCandProcess      = mcpMuonCandidate->Process();
    mcpMuonCandEndProcess   = mcpMuonCandidate->EndProcess();
    mcpMuonCandNDaughters   = mcpMuonCandidate->NumberDaughters();
    for (int i = 0; i < mcpMuonCandNDaughters; i++){
      mcpMuonCandDaughter[i] = mcpMuonCandidate->Daughter(i);
    }
    mcpMuonCandNumTrajPoint = mcpMuonCandidate->NumberTrajectoryPoints();
    mcpMuonCandVx           = mcpMuonCandidate->Vx();
    mcpMuonCandVy           = mcpMuonCandidate->Vy();
    mcpMuonCandVz           = mcpMuonCandidate->Vz();
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

    mcq0Transfer = mcNuEnergy - mcLeptonEnergy;
    mcq3Transfer = 
      std::sqrt(std::pow(mcNuPx-mcLeptonPx,2) +
          std::sqrt(std::pow(mcNuPy - mcLeptonPy, 2)) +
          std::sqrt(std::pow(mcNuPz - mcLeptonPz, 2)));

    mcNuMode            = mcNu.Mode();
    mcNuInteractionType = mcNu.InteractionType();

    /**
     * Get information on neutrino interaction
     */
    // is true nu CC?
    if (mcNuCCNC == 0) isCC = true;
    else isCC = false;

    // is true nu 0Pi?
    if (is0PiEvent(mcTruth)) is0Pi = true;
    else is0Pi = false;

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
    if (isGT2TrackAboveThresholdEvent(mcTruth, protonEThreshold) == 1) isTwoTrackTruth = true;
    else isTwoTrackTruth = false;

    // is 0 shower event?
    if (isZeroShowerEvent(mcTruth, electronEThreshold) == 1) isZeroShowerTruth = true;
    else isZeroShowerTruth = false;

    // is N protons above threshold greater than demand?
    if (isNPEvent(mcTruth, protonEThreshold, minNProtons) == 1) isGTNProtonsAboveThreshold = true;
    else isGTNProtonsAboveThreshold = false;

    // is true nu vertex in FV
    if (fiducialVolume.InFV(mcNuVx, mcNuVy, mcNuVz)) isTrueVtxInFV = true;
    else isTrueVtxInFV = false;

    if (isCC && is0Pi && isNuMu && isTrueVtxInFV && isGTNProtonsAboveThreshold) isSignal = true;
    else isSignal = false;

    if (isSignal && isBeamNeutrino && isSelected) isSelectedSignal = true;
    else isSelectedSignal = false;

    /**
     * Event is selected, put all information of interest in a tree.
     */
    if (isSelectedSignal){
      std::cout << "[NUMUSEL] Found and selected signal event!" << std::endl;
      ewutil.WriteTree(e, mcFlux, mcTruth, gTruth);
    }
    else if (isSelected){

      std::cout << "[NUMUSEL] Selected background event: " << std::endl;
      std::cout << "isCosmic: " << isCosmic << std::endl;
      std::cout << "isMixed:  " << isMixed << std::endl;
      std::cout << "isOOFV:   " << !isTrueVtxInFV << std::endl;
      std::cout << "isNC:     " << !isCC << std::endl;
      std::cout << "isNuMuBar:" << isNuMuBar << std::endl;
      std::cout << "isNuE:    " << isNuE << std::endl;
      std::cout << "isNuEBar: " << isNuEBar << std::endl;
      std::cout << "is0Pi:    " << !is0Pi << std::endl;
    }

  }

  tree->Fill();
}

void NuMuSelection::endJob()
{
  // Implementation of optional member function here.
}

void NuMuSelection::initialiseTree(TTree *t)
{
  t->Branch("fhicl_maxNShowers"          , &maxNShowers                );
  t->Branch("fhicl_minNTracks"           , &minNTracks                 );
  t->Branch("fhicl_minNProtons"          , &minNProtons                );
  t->Branch("fhicl_protonEThreshold"     , &protonEThreshold           );
  t->Branch("run"                        , &run                        );
  t->Branch("subrun"                     , &subrun                     );
  t->Branch("event"                      , &event                      );
  t->Branch("isSelected"                 , &isSelected                 );
  t->Branch("isCC"                       , &isCC                       );
  t->Branch("is0Pi"                      , &is0Pi                      );
  t->Branch("isNuMu"                     , &isNuMu                     );
  t->Branch("isNuMuBar"                  , &isNuMuBar                  );
  t->Branch("isNuE"                      , &isNuE                      );
  t->Branch("isNuEBar"                   , &isNuEBar                   );
  t->Branch("isTrueVtxInFV"              , &isTrueVtxInFV              );
  t->Branch("isGTNProtonsAboveThreshold" , &isGTNProtonsAboveThreshold );
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
  t->Branch("mcpMuonCandDaughter"        , mcpMuonCandDaughter         , "mcpMuonCandDaughter/I");
  t->Branch("mcpMuonCandNumTrajPoint"    , &mcpMuonCandNumTrajPoint    );
  t->Branch("mcpMuonCandVx"              , &mcpMuonCandVx              );
  t->Branch("mcpMuonCandVy"              , &mcpMuonCandVy              );
  t->Branch("mcpMuonCandVz"              , &mcpMuonCandVz              );
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
  t->Branch("nTracks"                    , &nTracks                    );
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

bool NuMuSelection::isGT2TrackAboveThresholdEvent(art::Ptr<simb::MCTruth> mcTruth, double protonEThreshold){

  int nParticles = mcTruth->NParticles();
  int particleCounter = 0;

  for (int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);
    if (particle.Process() != "primary" || particle.StatusCode() !=1) continue;

    if (std::abs(particle.PdgCode()) == 13  || 
        std::abs(particle.PdgCode()) == 211 ||
        (std::abs(particle.PdgCode()) == 2212 && particle.E() > protonEThreshold))
      particleCounter++;

  }

  if (particleCounter >= 2) return true;
  else return false;

}

bool NuMuSelection::isZeroShowerEvent(art::Ptr<simb::MCTruth> mcTruth, double electronEThreshold){

  int nParticles = mcTruth->NParticles();
  int particleCounter = 0;

  for (int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);
    if (particle.Process() != "primary" || particle.StatusCode() !=1) continue;

    if (std::abs(particle.PdgCode()) == 22  || 
        (std::abs(particle.PdgCode()) == 11 && particle.E() > electronEThreshold))
      particleCounter++;

  }

  if (particleCounter == 0 ) return true;
  else return false;

}

bool NuMuSelection::isNPEvent(art::Ptr<simb::MCTruth> mcTruth, double protonEThreshold, int minNProtons)
{

  int nParticles = mcTruth->NParticles();
  int nProtonsAboveThreshold = 0;
  for (int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);
    if (particle.Process() != "primary" || particle.StatusCode() !=1) continue;

    if (std::abs(particle.PdgCode()) == 2212 && particle.E() > protonEThreshold){
      nProtonsAboveThreshold++; 
    }
  }

  if (nProtonsAboveThreshold >= minNProtons)
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
