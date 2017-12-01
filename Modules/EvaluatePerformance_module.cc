////////////////////////////////////////////////////////////////////////
// Class:       EvaluatePerformance
// Plugin Type: analyzer (art v2_05_00)
// File:        EvaluatePerformance_module.cc
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
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft Includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"

// ROOT Includes
#include "TTree.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"

// UBXSec Includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"

// FMWK Includes

class EvaluatePerformance;


class EvaluatePerformance : public art::EDAnalyzer {
  public:
    explicit EvaluatePerformance(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    EvaluatePerformance(EvaluatePerformance const &) = delete;
    EvaluatePerformance(EvaluatePerformance &&) = delete;
    EvaluatePerformance & operator = (EvaluatePerformance const &) = delete;
    EvaluatePerformance & operator = (EvaluatePerformance &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    bool isCC0PiEvent(art::Ptr<simb::MCTruth> mcTruth);

  private:

    art::ServiceHandle< geo::Geometry > geo;
    art::ServiceHandle< art::TFileService > tfs;
    TTree* selectionEfficiency;
    ubana::FiducialVolume fiducialVolume;  

    // fcl vars
    std::string fSelectionLabel;

    // vars
    
    int run;
    int subRun;
    int event;

    bool isEventPassed;
    bool isCC0Pi;
    bool isBeamNeutrino;
    bool isCosmic;
    bool isMixed;
    bool isCC;
    bool isMuonNeutrino;
    bool isMuonAntiNeutrino;
    bool isElectronNeutrino;
    bool isElectronAntineutrino;
    bool isTrueVtxInFV;
    bool isSignal;
    int mcNuCCNC;

    double vx, vy, vz;

    double mcNuEnergy;
    double mcNuMom;
    double mcNuPx, mcNuPy, mcNuPz;
    double mcLeptonEnergy;
    double mcLeptonMom;
    double mcLeptonPx, mcLeptonPy, mcLeptonPz;
    double mcLeptonTheta;
    double mcLeptonPhi;

    double momentumTransfer;
    double energyTransfer;

    double muonCandidateLength;

    std::vector<bool> purityVector = {0,0,0,0,0,0,0,0};

    // Efficiency histograms
    TEfficiency* mcNuEnergyEff;
    TEfficiency* mcLeptonMomEff;
    TEfficiency* mcLeptonThetaEff;
    TEfficiency* mcLeptonPhiEff;

    TEfficiency* mcNuEnergyCC0PiEff;
    TEfficiency* mcLeptonMomCC0PiEff;
    TEfficiency* mcLeptonThetaCC0PiEff;
    TEfficiency* mcLeptonPhiCC0PiEff;

    // Purity histograms
    TEfficiency* mcNuEnergyPur;
    TEfficiency* mcLeptonMomPur;
    TEfficiency* mcLeptonThetaPur;
    TEfficiency* mcLeptonPhiPur;

    TEfficiency* mcNuEnergyCC0PiPur;
    TEfficiency* mcLeptonMomCC0PiPur;
    TEfficiency* mcLeptonThetaCC0PiPur;
    TEfficiency* mcLeptonPhiCC0PiPur;

    // length histograms
    TH1D* selectedTrackLengthTrueCC;
    TH1D* selectedTrackLengthTrueCC0Pi;
    TH1D* selectedTrackLengthCosmic;
    TH1D* selectedTrackLengthMixed;
    TH1D* selectedTrackLengthOutOfFV;
    TH1D* selectedTrackLengthAntiNuMu;
    TH1D* selectedTrackLengthNuE;
    TH1D* selectedTrackLengthAntiNuE;
    TH1D* selectedTrackLengthNC;
    TH1D* selectedTrackLengthOther;

    // truth histograms
    TH2D* trueq3q0;
};


EvaluatePerformance::EvaluatePerformance(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
{

  // configure
  fiducialVolume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
      geo->DetHalfHeight(),
      2.*geo->DetHalfWidth(),
      geo->DetLength());

  fiducialVolume.PrintConfig();

  // get fcl params
  fSelectionLabel = p.get<std::string> ("SelectionLabel");
}

void EvaluatePerformance::analyze(art::Event const & e)
{
  run = e.run();
  subRun = e.subRun();
  event = e.event();

  std::cout << "\n>| PRINTING INFORMATION ABOUT EVENT " << run << "." << subRun << "." << event << std::endl;;

  // MC truth
  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
  e.getByLabel("generator", mcTruthHandle);
  if (!mcTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
  if (mcTruthVec.size() == 0){ 
    std::cout << ">|>| No Neutrinos " << std::endl;
    return;
  }

  // selection info
  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel(fSelectionLabel, selectionHandle);
  if (!selectionHandle.isValid()){

    std::cout << ">|>| SelectionResult product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found."
      << std::endl;
    throw std::exception(); 

  }

  // get Selected TPCObject from selection handle
  art::FindManyP<ubana::TPCObject> selectedTpcObjects(selectionHandle, e, fSelectionLabel);
  art::Ptr<ubana::TPCObject> selectedTpcObject; 

  // get objects from handle
  art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);
  const simb::MCNeutrino& mcNu = mcTruth->GetNeutrino();
  const simb::MCParticle& mcNuP = mcNu.Nu();
  const simb::MCParticle& mcLeptonP = mcNu.Lepton();

  double vx = (double)mcNuP.Vx();
  double vy = (double)mcNuP.Vy();
  double vz = (double)mcNuP.Vz();
  mcNuCCNC  = mcNu.CCNC();

  // variables for making efficiency plots
  mcNuEnergy = mcNuP.E();
  mcNuPx = mcNuP.Px();
  mcNuPy = mcNuP.Py();
  mcNuPz = mcNuP.Pz();
  mcLeptonEnergy = mcLeptonP.E();
  mcLeptonMom = mcLeptonP.P();
  mcLeptonPx = mcLeptonP.Px();
  mcLeptonPy = mcLeptonP.Py();
  mcLeptonPz = mcLeptonP.Pz();
  mcLeptonTheta = std::cos(mcLeptonP.Momentum().Theta());
  mcLeptonPhi = mcLeptonP.Momentum().Phi();

  // -------------------------------------------------------------------------
  // Truth information
  // --- "is there a true numu CC interaction in the FV for this event"
  // -------------------------------------------------------------------------

  // is true nu CC?
  if (mcNuCCNC == 0) isCC = true;
  else isCC = false;

  // is true nu CC0Pi?
  if (isCC0PiEvent(mcTruth)) isCC0Pi = true;
  else isCC0Pi = false;

  // is true nu muon nu interaction?
  if (mcNuP.PdgCode() == 14) isMuonNeutrino = true;
  else isMuonNeutrino = false;

  // else check if muon antineutrino...
  if (mcNuP.PdgCode() == -14) isMuonAntiNeutrino = true;
  else isMuonAntiNeutrino = false;

  // ... or electron neutrino
  if (mcNuP.PdgCode() == 12) isElectronNeutrino = true;
  else isElectronNeutrino = false;

  // ... or electron antineutrino
  if (mcNuP.PdgCode() == -12) isElectronAntineutrino = true;
  else isElectronAntineutrino = false;

  // is true nu interaction in FV
  if (fiducialVolume.InFV(vx, vy, vz)) isTrueVtxInFV = true;
  else isTrueVtxInFV = false;

  if (isCC && isCC0Pi && isMuonNeutrino && isTrueVtxInFV) isSignal = true;
  else isSignal = false;

  //
  // Calculate energy and momentum transfer
  //

  // -------------------------------------------------------------------------
  // Reconstructed information
  // --- "did we select an object in this event"
  // --- "if we did, is it a neutrino interaction"
  // -------------------------------------------------------------------------


  if (selectedTpcObjects.at(0).size() == 1){
    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const ubana::TPCObjectOrigin& selectedOrigin = selectedTpcObject->GetOrigin();
    const std::vector<recob::Track>& selectedTracks = selectedTpcObject->GetTracks();

    // get candidate muon vars
    muonCandidateLength = 0;
    for (size_t i = 0; i < selectedTracks.size(); i++){

      const recob::Track& track = selectedTracks.at(i);
      if (track.Length() > muonCandidateLength){

        muonCandidateLength = track.Length();

      }

    }

    // is selected TPC object a true beam neutrino interaction
    isBeamNeutrino = false;
    isMixed = false;
    isCosmic = false;

    if (selectedOrigin == ubana::kBeamNeutrino) isBeamNeutrino = true;

    if (selectedOrigin == ubana::kMixed) isMixed = true;

    if (selectedOrigin == ubana::kCosmicRay) isCosmic = true;


  }

  // check if selected and fill efficiency histogram
  for (auto const& selectionStatus : (*selectionHandle)){

    if (selectionStatus.GetSelectionStatus() == true){

      std::cout << ">|>| Event Selected." << std::endl;
      if (isBeamNeutrino && isSignal){
        // event passes

        std::cout << ">|>|>| Event is true numuCC from beam in FV" << std::endl;

        isEventPassed = true;

        // efficiency
        mcNuEnergyEff->Fill(isEventPassed, mcNuEnergy);
        mcLeptonMomEff->Fill(isEventPassed, mcLeptonMom);
        mcLeptonThetaEff->Fill(isEventPassed, mcLeptonTheta);
        mcLeptonPhiEff->Fill(isEventPassed, mcLeptonPhi);

        // purity
        mcNuEnergyPur->Fill(true, mcNuEnergy);
        mcLeptonMomPur->Fill(true, mcLeptonMom);
        mcLeptonThetaPur->Fill(true, mcLeptonTheta);
        mcLeptonPhiPur->Fill(true, mcLeptonPhi);

        selectedTrackLengthTrueCC->Fill(muonCandidateLength);
        if (isCC0Pi){

          // efficiency
          mcNuEnergyCC0PiEff->Fill(isEventPassed, mcNuEnergy);
          mcLeptonMomCC0PiEff->Fill(isEventPassed, mcLeptonMom);
          mcLeptonThetaCC0PiEff->Fill(isEventPassed, mcLeptonTheta);
          mcLeptonPhiCC0PiEff->Fill(isEventPassed, mcLeptonPhi);

          // purity
          mcNuEnergyCC0PiPur->Fill(true, mcNuEnergy);
          mcLeptonMomCC0PiPur->Fill(true, mcLeptonMom);
          mcLeptonThetaCC0PiPur->Fill(true, mcLeptonTheta);
          mcLeptonPhiCC0PiPur->Fill(true, mcLeptonPhi);

          selectedTrackLengthTrueCC0Pi->Fill(muonCandidateLength);

        }
      
        energyTransfer = mcNuEnergy - mcLeptonEnergy;
        momentumTransfer = std::sqrt(std::pow(mcNuPx - mcLeptonPx,2) + std::pow(mcNuPy - mcLeptonPy,2) + std::pow(mcNuPz - mcLeptonPz,2));

        trueq3q0->Fill(momentumTransfer, energyTransfer);

      }
      else{ 
        std::cout << ">|>|>| Event is not true numuCC from beam in FV" << std::endl;

        purityVector.at(0) = isCosmic;
        purityVector.at(1) = isMixed;
        purityVector.at(2) = !isTrueVtxInFV;
        purityVector.at(3) = isMuonNeutrino;
        purityVector.at(4) = isMuonAntiNeutrino;
        purityVector.at(5) = isElectronNeutrino;
        purityVector.at(6) = isElectronAntineutrino;
        purityVector.at(7) = !isCC;

        std::cout << ">|>|>| FAILURE REASON(S)" << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(0) << " isCosmic" << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(1) << " isMixed" << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(2) << " isTrueVtxOutOfFV" << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(3) << " isMuonNeutrino " << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(4) << " isMuonAntiNeutrino" << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(5) << " isElectronNeutrino" << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(6) << " isElectronAntiNeutrino" << std::endl;
        std::cout << ">|>|>|>| " << purityVector.at(7) << " isNC" << std::endl;

        // purity 

        mcNuEnergyPur->Fill(false, mcNuEnergy);
        mcLeptonMomPur->Fill(false, mcLeptonMom);
        mcLeptonThetaPur->Fill(false, mcLeptonTheta);
        mcLeptonPhiPur->Fill(false, mcLeptonPhi);

        mcNuEnergyCC0PiPur->Fill(false, mcNuEnergy);
        mcLeptonMomCC0PiPur->Fill(false, mcLeptonMom);
        mcLeptonThetaCC0PiPur->Fill(false, mcLeptonTheta);
        mcLeptonPhiCC0PiPur->Fill(false, mcLeptonPhi);

        std::cout << "muonCandLen " << muonCandidateLength << std::endl;

        if (isCosmic){ selectedTrackLengthCosmic->Fill(muonCandidateLength); }
        else if (isMixed){ selectedTrackLengthMixed->Fill(muonCandidateLength); }
        else if (!isTrueVtxInFV){ selectedTrackLengthOutOfFV->Fill(muonCandidateLength); }
        else if (isMuonNeutrino){ selectedTrackLengthNC->Fill(muonCandidateLength); }
        else if (isMuonAntiNeutrino){ selectedTrackLengthAntiNuMu->Fill(muonCandidateLength); }
        else if (isElectronNeutrino){ selectedTrackLengthNuE->Fill(muonCandidateLength); }
        else if (isElectronAntineutrino){ selectedTrackLengthAntiNuE->Fill(muonCandidateLength); }
        else {selectedTrackLengthOther->Fill(muonCandidateLength);}

      }

    }
    else {

      // event does not pass

      std::cout << ">|>| Not Selected." << std::endl;

      isEventPassed = false;
      if (isSignal){

        std::cout << ">|>|>| Event is true numuCC from beam in FV" << std::endl;
        mcNuEnergyEff->Fill(isEventPassed, mcNuEnergy);
        mcLeptonMomEff->Fill(isEventPassed, mcLeptonMom);
        mcLeptonThetaEff->Fill(isEventPassed, mcLeptonTheta);
        mcLeptonPhiEff->Fill(isEventPassed, mcLeptonPhi);

        if (isCC0Pi){

          mcNuEnergyCC0PiEff->Fill(isEventPassed, mcNuEnergy);
          mcLeptonMomCC0PiEff->Fill(isEventPassed, mcLeptonMom);
          mcLeptonThetaCC0PiEff->Fill(isEventPassed, mcLeptonTheta);
          mcLeptonPhiCC0PiEff->Fill(isEventPassed, mcLeptonPhi);


        }

      }
      else std::cout << ">|>|>| Event is not true numuCC from beam in FV" << std::endl;

    }
  }

  selectionEfficiency->Fill();
}

void EvaluatePerformance::beginJob()
{
  // Implementation of optional member function here.
  selectionEfficiency = tfs->make<TTree>("selectionEfficiency", "selectionEfficiency");
  selectionEfficiency->Branch("purityVector", &purityVector, "purityVector");  

  // efficiencies

  mcNuEnergyEff = tfs->make<TEfficiency>("mcNuEnergyEff", ";E_{#nu}^{true}; #epsilon", 15, 0, 3);
  mcLeptonMomEff = tfs->make<TEfficiency>("mcLeptonMomEff", ";P_{l}^{true}; #epsilon", 15, 0, 2);
  mcLeptonThetaEff = tfs->make<TEfficiency>("mcLeptonThetaEff", ";cos(#theta_{l}^{true}); #epsilon;", 10, -1, 1);
  mcLeptonPhiEff = tfs->make<TEfficiency>("mcLeptonPhiEff", ";#phi_{l}^{true}; #epsilon", 15, -3, 3);
  mcNuEnergyCC0PiEff = tfs->make<TEfficiency>("mcNuEnergyCC0PiEff", ";E_{#nu}^{true}; #epsilon", 15, 0, 3);
  mcLeptonMomCC0PiEff = tfs->make<TEfficiency>("mcLeptonMomCC0PiEff", ";P_{l}^{true}; #epsilon", 15, 0, 2);
  mcLeptonThetaCC0PiEff = tfs->make<TEfficiency>("mcLeptonThetaCC0PiEff", ";cos(#theta_{l}^{true}); #epsilon;", 10, -1, 1);
  mcLeptonPhiCC0PiEff = tfs->make<TEfficiency>("mcLeptonPhiCC0PiEff", ";#phi_{l}^{true};#epsilon", 15, -3, 3);

  // purity

  mcNuEnergyPur = tfs->make<TEfficiency>("mcNuEnergyPur", ";E_{#nu}^{true};p", 15, 0, 3);
  mcLeptonMomPur = tfs->make<TEfficiency>("mcLeptonMomPur", ";P_{l}^{true};p", 15, 0, 2);
  mcLeptonThetaPur = tfs->make<TEfficiency>("mcLeptonThetaPur", ";cos(#theta_{l}^{true});p", 10, -1, 1);
  mcLeptonPhiPur = tfs->make<TEfficiency>("mcLeptonPhiPur", ";#phi_{l}^{true};p", 15, -3, 3);

  mcNuEnergyCC0PiPur = tfs->make<TEfficiency>("mcNuEnergyCC0PiPur", ";E_{#nu}^{true};p", 15, 0, 3);
  mcLeptonMomCC0PiPur = tfs->make<TEfficiency>("mcLeptonMomCC0PiPur", ";P_{l}^{true};p", 15, 0, 2);
  mcLeptonThetaCC0PiPur = tfs->make<TEfficiency>("mcLeptonThetaCC0PiPur", ";cos(#theta_{l}^{true});p", 10, -1, 1);
  mcLeptonPhiCC0PiPur = tfs->make<TEfficiency>("mcLeptonPhiCC0PiPur", ";#phi_{l}^{true};p", 15, -3, 3);

  // length histograms
  selectedTrackLengthTrueCC = tfs->make<TH1D>("selectedTrackLengthTrueCC", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthTrueCC0Pi = tfs->make<TH1D>("selectedTrackLengthTrueCC0Pi", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthCosmic = tfs->make<TH1D>("selectedTrackLengthCosmic", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthMixed = tfs->make<TH1D>("selectedTrackLengthMixed", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthOutOfFV = tfs->make<TH1D>("selectedTrackLengthOutOfFV", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthAntiNuMu = tfs->make<TH1D>("selectedTrackLengthAntiNuMu", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthNuE = tfs->make<TH1D>("selectedTrackLengthNuE", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthAntiNuE = tfs->make<TH1D>("selectedTrackLengthAntiNuE", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthNC = tfs->make<TH1D>("selectedTrackLengthNC", ";Candidate muon length; # tracks", 30, 0, 700);
  selectedTrackLengthOther = tfs->make<TH1D>("selectedTrackLengthOther", ";Candidate muon length; # tracks", 30, 0, 700);

  // truth histograms
  trueq3q0 = tfs->make<TH2D>("trueq3q0", ";True q3 (GeV);True q0 (GeV)", 50, 0, 1.2, 50, 0, 1.2);

}

void EvaluatePerformance::endJob()
{


}

bool EvaluatePerformance::isCC0PiEvent(art::Ptr<simb::MCTruth> mcTruth) {

  int nParticles = mcTruth->NParticles();
  for(int i = 0; i < nParticles; i++){

    const simb::MCParticle& particle = mcTruth->GetParticle(i);

    if (particle.Process() != "primary" || particle.StatusCode() != 1) continue;    

    if (std::abs(particle.PdgCode()) == 211 || std::abs(particle.PdgCode()) == 111){

      return false;

    }
  }

  return true;
}

DEFINE_ART_MODULE(EvaluatePerformance)
