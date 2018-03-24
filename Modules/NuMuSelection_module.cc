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

    /**
     * initialises branches in output tree
     */
    void initialiseTree(TTree *t);

    /**
     * check whether event is truly CC0Pi
     */
    bool is0PiEvent(art::Ptr<simb::MCTruth> mcTruth);

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

    // vars
    int fRun;
    int fSubRun;
    int fEvent;

    // truth
    double mcNuVx;
    double mcNuVy;
    double mcNuVz;
    int mcNuCCNC;
    int mcNuMode;
    int mcNuInteractionType;
    double mcNuEnergy;
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
    double mcq3Transfer;
    double mcq0Transfer;

    // reco
    double nTracks;
    double candMuonLength;
    double candMuonTheta;
    double candMuonCosTheta;
    double candMuonPhi;
    double candMuonEndX;
    double candMuonEndY;
    double candMuonEndZ;
    double candMuonStartX;
    double candMuonStartY;
    double candMuonStartZ;
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

  // hit information
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, "pandoraNu");

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
    << fRun << "." << fSubRun << "." << fEvent << std::endl;

  /**
   * Selection 
   *
   * Demand:
   * -- Only one selected TPC object
   * -- Only one track is considered muon-like
   * -- No showers reconstructed in the shower object
   */
 
  recob::Track muonCandidate;
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

    // Check origin of TPC object
    if (selectedOrigin == ubana::kBeamNeutrino) isBeamNeutrino = true;
    else isBeamNeutrino = false;

    if (selectedOrigin == ubana::kMixed) isMixed = true;
    else isMixed = false;

    if (selectedOrigin == ubana::kCosmicRay) isCosmic = true;
    else isCosmic = false;

    int muonLikeCounter = 0;

    /**
     * Loop tracks and check whether each track is compatibile with
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
        pidutils.getAveragedQdX(track, hits);

      bool isMuLike = 
        pidutils.isMuonLike(averageDqdx.first, averageDqdx.second);

      if (isMuLike){ 
        muonLikeCounter++;
        muonCandidate = track;
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

    /*
     * Topology cuts
     */
    if ((int)nSelectedShowers > maxNShowers){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("nShowersGTMaxNShowers");
    }

    if ((int)selectedTracks.size() < minNTracks){
      selResult.SetSelectionStatus(false);
      selResult.SetFailureReason("nTracksLTMinNTracks");
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

  e.put(std::move(selectionCollection));
  e.put(std::move(tpcObjectCollection));
  e.put(std::move(selTpcObjAssn));

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
    candMuonLength           = muonCandidate.Length();
    candMuonTheta            = muonCandidate.Theta();
    candMuonCosTheta         = std::cos(muonCandidate.Theta());
    candMuonPhi              = muonCandidate.Phi();
    candMuonStartX           = muonCandidate.Vertex().X();
    candMuonStartY           = muonCandidate.Vertex().Y();
    candMuonStartZ           = muonCandidate.Vertex().Z();
    candMuonEndX             = muonCandidate.End().X();
    candMuonEndY             = muonCandidate.End().Y();
    candMuonEndZ             = muonCandidate.End().Z();
    nTracks                  = (int)nSelectedTracks;
    double xyz[3];
    selectedVertex.XYZ(xyz);
    vertexX                  = xyz[0];
    vertexY                  = xyz[1];
    vertexZ                  = xyz[2];
  }
  else {
    candMuonLength           = -999;
    candMuonTheta            = -999;
    candMuonCosTheta         = -999;
    candMuonPhi              = -999;
    candMuonStartX           = -999;
    candMuonStartY           = -999;
    candMuonStartZ           = -999;
    candMuonEndX             = -999;
    candMuonEndY             = -999;
    candMuonEndZ             = -999;
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
      const simb::MCNeutrino& mcNu = mcTruth->GetNeutrino();
      const simb::MCParticle& mcNuP = mcNu.Nu();
      const simb::MCParticle& mcLeptonP = mcNu.Lepton();
    
      mcNuVx = (double) mcNuP.Vx();
      mcNuVy = (double) mcNuP.Vy();
      mcNuVz = (double) mcNuP.Vz();
      mcNuPx = mcNuP.Px();
      mcNuPy = mcNuP.Py();
      mcNuPz = mcNuP.Pz();
      mcNuEnergy = mcNuP.E();
      mcLeptonMom = mcLeptonP.P();
      mcLeptonPx = mcLeptonP.Px();
      mcLeptonPy = mcLeptonP.Py();
      mcLeptonPz = mcLeptonP.Pz();
      mcLeptonTheta = mcLeptonP.Momentum().Theta();
      mcLeptonCosTheta = std::cos(mcLeptonP.Momentum().Theta());
      mcLeptonPhi = mcLeptonP.Momentum().Phi();

      mcq0Transfer = mcNuEnergy - mcLeptonEnergy;
      mcq3Transfer = 
        std::sqrt(std::pow(mcNuPx-mcLeptonPx,2) +
            std::sqrt(std::pow(mcNuPy - mcLeptonPy, 2)) +
            std::sqrt(std::pow(mcNuPz - mcLeptonPz, 2)));

      mcNuMode = mcNu.Mode();
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

    // is true nu vertex in FV
    if (fiducialVolume.InFV(mcNuVx, mcNuVy, mcNuVz)) isTrueVtxInFV = true;
    else isTrueVtxInFV = false;
  
    if (isCC && is0Pi && isNuMu && isTrueVtxInFV) isSignal = true;
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
      std::cout << "isNPi:    " << !is0Pi << std::endl;
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
  t->Branch("isSelected", &isSelected);
  t->Branch("isCC", &isCC);
  t->Branch("is0Pi", &is0Pi);
  t->Branch("isNuMu", &isNuMu);
  t->Branch("isNuMuBar", &isNuMuBar);
  t->Branch("isNuE", &isNuE);
  t->Branch("isNuEBar", &isNuEBar);
  t->Branch("isTrueVtxInFV", &isTrueVtxInFV);
  t->Branch("isSignal", &isSignal);
  t->Branch("isBeamNeutrino", &isBeamNeutrino);
  t->Branch("isMixed", &isMixed);
  t->Branch("isCosmic", &isCosmic);
  t->Branch("isSelectedSignal", &isSelectedSignal);
  t->Branch("mcNuCCNC", &mcNuCCNC);
  t->Branch("mcNuMode", &mcNuMode);
  t->Branch("mcNuInteractionType", &mcNuInteractionType);
  t->Branch("mcNuVx", &mcNuVx);
  t->Branch("mcNuVy", &mcNuVy);
  t->Branch("mcNuVz", &mcNuVz);
  t->Branch("mcNuPx", &mcNuPx);
  t->Branch("mcNuPy", &mcNuPy);
  t->Branch("mcNuPz", &mcNuPz);
  t->Branch("mcNuEnergy", &mcNuEnergy);
  t->Branch("mcNuTheta", &mcNuTheta);
  t->Branch("mcLeptonMom", &mcLeptonMom);
  t->Branch("mcLeptonPx", &mcLeptonPx);
  t->Branch("mcLeptonPy", &mcLeptonPy);
  t->Branch("mcLeptonPz", &mcLeptonPz);
  t->Branch("mcLeptonEnergy", &mcLeptonEnergy);
  t->Branch("mcLeptonTheta", &mcLeptonTheta);
  t->Branch("mcLeptonCosTheta", &mcLeptonCosTheta);
  t->Branch("mcLeptonPhi", &mcLeptonPhi);
  t->Branch("nTracks", &nTracks);
  t->Branch("candMuonLength", &candMuonLength);
  t->Branch("candMuonTheta", &candMuonTheta);
  t->Branch("candMuonCosTheta", & candMuonCosTheta);
  t->Branch("candMuonPhi", &candMuonPhi);
  t->Branch("candMuonStartX", &candMuonStartX);
  t->Branch("candMuonStartY", &candMuonStartY);
  t->Branch("candMuonStartZ", &candMuonStartZ);
  t->Branch("candMuonEndX", &candMuonEndX);
  t->Branch("candMuonEndY", &candMuonEndY);
  t->Branch("candMuonEndZ", &candMuonEndZ);
  t->Branch("vertexX", &vertexX);
  t->Branch("vertexY", &vertexY);
  t->Branch("vertexZ", &vertexZ);
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

DEFINE_ART_MODULE(NuMuSelection)
