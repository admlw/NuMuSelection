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

  private:

    art::ServiceHandle< geo::Geometry > geo;
    art::ServiceHandle< art::TFileService > tfs;
    TTree* selectionEfficiency;
    ubana::FiducialVolume fiducialVolume;  

    // vars
    bool isEventPassed;
    bool isCC0pi;
    bool isBeamNeutrino;
    bool isCosmic;
    bool isMixed;
    bool isRecoInFV;
    bool isCC;
    bool isMuonNeutrino;
    bool isMuonAntineutrino;
    bool isElectronNeutrino;
    bool isElectronAntineutrino;
    bool isInFV;
    int mcNuCCNC;

    double vx, vy, vz;

    double mcNuEnergy;
    double mcLeptonMom;
    double mcLeptonTheta;
    double mcLeptonPhi;

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

}

void EvaluatePerformance::analyze(art::Event const & e)
{

  // MC truth
  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
  e.getByLabel("generator", mcTruthHandle);
  if (!mcTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);
  if (mcTruthVec.size() == 0){ 
    std::cout << " >> No Neutrinos " << std::endl;
    return;
  }

  // selection info
  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel("UBXSec", selectionHandle);
  if (!selectionHandle.isValid()){

    std::cout << " >> SelectionResult product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found."
      << std::endl;
    throw std::exception(); 

  }

  // get Selected TPCObject from selection handle
  art::FindManyP<ubana::TPCObject> selectedTpcObjects(selectionHandle, e, "UBXSec");
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
  mcLeptonMom = mcLeptonP.P();
  mcLeptonTheta = std::cos(mcLeptonP.Momentum().Theta());
  mcLeptonPhi = mcLeptonP.Momentum().Phi();

  // ------------------
  // Truth-based cuts
  // ------------------

  // is true nu CC?
  if (mcNuCCNC != 0){

    std::cout << " >> Event is NC " << std::endl;
    isCC = false;

  }
  else isCC = true;

  // is true nu muon nu interaction?
  if (mcNuP.PdgCode() != 14){

    std::cout << " >> not a nu_mu interaction!" << std::endl;
    isMuonNeutrino = false;

  }
  else isMuonNeutrino = true;

  if (mcNuP.PdgCode() == -14) isMuonAntineutrino = true;
  else isMuonAntineutrino = false;

  if (mcNuP.PdgCode() == 12) isElectronNeutrino = true;
  else isElectronNeutrino = false;

  if (mcNuP.PdgCode() == -12) isElectronAntineutrino = true;
  else isElectronAntineutrino = false;

  // is true nu interaction in FV
  if (!fiducialVolume.InFV(vx, vy, vz)){

    std::cout << " >> true interaction vertex " << vx << ", " << vy << ", " << vz
      << "is out of TPC! " << std::endl;
    isInFV = false;

  }
  else isInFV = true;

  // ---------------------
  // Reco-level cuts
  // ---------------------

  if (selectedTpcObjects.at(0).size() != 0){
    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const recob::Vertex& selectedVertex = selectedTpcObject->GetVertex();
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

    // get reconstructed vertex position
    // if true vertex is in the fiducial volume but reconstructed is out, this goes in the denominator
    double xyz[3] = {0.0,0.0,0.0};
    selectedVertex.XYZ(xyz);
    if (fiducialVolume.InFV(xyz[0], xyz[1], xyz[2]) == true){

      std::cout << " >> Reconstruced object is in TPC" << std::endl;
      isRecoInFV = true;

    }
    else{
      std::cout << " >> Reconstructed object is out of TPC!" << std::endl;
      isRecoInFV = false;
    }


    // get selected TPCObject  origin
    if ( selectedOrigin == ubana::kBeamNeutrino){

      std::cout << " >> Found a beam neutrino! " << std::endl; 
      isBeamNeutrino = true;
      isCosmic = false;
      isMixed = false;

    }
    else if (selectedOrigin == ubana::kCosmicRay) {

      std::cout << " >> Not a beam neutrino!" << std::endl;
      isBeamNeutrino = false;
      isCosmic = true;
      isMixed = false;

    }
    else if (selectedOrigin == ubana::kMixed){

      std::cout << " >> Not a beam neutrino!" << std::endl;
      isBeamNeutrino = false;
      isCosmic = false;
      isMixed = true;

    }
    else isBeamNeutrino = false;
  }

  // check if selected and fill efficiency histogram
  for (auto const& selectionStatus : (*selectionHandle)){

    if (isCC && isMuonNeutrino && isInFV){
      std::cout << " >> Found a true CC Event " << std::endl;

      int nParticles = mcTruth->NParticles();
      for(int i = 0; i < nParticles; i++){

        const simb::MCParticle& particle = mcTruth->GetParticle(i);
        if (particle.Process() != "primary" || particle.StatusCode() != 1) continue;    

        if (std::abs(particle.PdgCode()) == 211 || std::abs(particle.PdgCode()) == 111){

          isCC0pi = false;

        }
        else isCC0pi = true;

      }
    }

    if (selectionStatus.GetSelectionStatus() == true){

      if (isBeamNeutrino == true && isRecoInFV == true && isCC == true && isMuonNeutrino == true && isInFV == true){
        // event passes

        std::cout << " -->> Selected." << std::endl;
        isEventPassed = true;

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

        if (isCC0pi == true){

          mcNuEnergyCC0PiEff->Fill(isEventPassed, mcNuEnergy);
          mcLeptonMomCC0PiEff->Fill(isEventPassed, mcLeptonMom);
          mcLeptonThetaCC0PiEff->Fill(isEventPassed, mcLeptonTheta);
          mcLeptonPhiCC0PiEff->Fill(isEventPassed, mcLeptonPhi);
        
          mcNuEnergyCC0PiPur->Fill(true, mcNuEnergy);
          mcLeptonMomCC0PiPur->Fill(true, mcLeptonMom);
          mcLeptonThetaCC0PiPur->Fill(true, mcLeptonTheta);
          mcLeptonPhiCC0PiPur->Fill(true, mcLeptonPhi);

        }

      }
      else{ 
        std::cout << "selected event is not signal!" << std::endl;

        purityVector.at(0) = isCosmic;
        purityVector.at(1) = isMixed;
        purityVector.at(2) = !(isInFV && !isRecoInFV);
        purityVector.at(3) = isMuonAntineutrino;
        purityVector.at(4) = isElectronNeutrino;
        purityVector.at(5) = isElectronAntineutrino;
        purityVector.at(6) = !isCC;

        // purity 
        
        mcNuEnergyPur->Fill(false, mcNuEnergy);
        mcLeptonMomPur->Fill(false, mcLeptonMom);
        mcLeptonThetaPur->Fill(false, mcLeptonTheta);
        mcLeptonPhiPur->Fill(false, mcLeptonPhi);

        mcNuEnergyCC0PiPur->Fill(false, mcNuEnergy);
        mcLeptonMomCC0PiPur->Fill(false, mcLeptonMom);
        mcLeptonThetaCC0PiPur->Fill(false, mcLeptonTheta);
        mcLeptonPhiCC0PiPur->Fill(false, mcLeptonPhi);

        for ( size_t i = 0; i < purityVector.size(); i++ ){

          std::cout << " >> " << purityVector.at(i) << std::endl;;

        }

        if (isCosmic){ selectedTrackLengthCosmic->Fill(muonCandidateLength); continue;}
        if (isMixed){ selectedTrackLengthMixed->Fill(muonCandidateLength); continue;}
        if (!(isInFV && !isRecoInFV)){ selectedTrackLengthOutOfFV->Fill(muonCandidateLength); continue;}
        if (isMuonAntineutrino){ selectedTrackLengthAntiNuMu->Fill(muonCandidateLength); continue;}
        if (isElectronNeutrino){ selectedTrackLengthNuE->Fill(muonCandidateLength); continue;}
        if (isElectronAntineutrino){ selectedTrackLengthAntiNuE->Fill(muonCandidateLength); continue;}
        if (isCC){ selectedTrackLengthNC->Fill(muonCandidateLength); continue;}

      }

    }
    else {

      // event does not pass
      isEventPassed = false;
      if (isBeamNeutrino == true && isRecoInFV == true && isCC == true && isMuonNeutrino == true && isInFV == true){

        mcNuEnergyEff->Fill(isEventPassed, mcNuEnergy);
        mcLeptonMomEff->Fill(isEventPassed, mcLeptonMom);
        mcLeptonThetaEff->Fill(isEventPassed, mcLeptonTheta);
        mcLeptonPhiEff->Fill(isEventPassed, mcLeptonPhi);

        std::cout << " -->>Not Selected." << std::endl;

        if (isCC0pi == true){
          mcNuEnergyCC0PiEff->Fill(isEventPassed, mcNuEnergy);
          mcLeptonMomCC0PiEff->Fill(isEventPassed, mcLeptonMom);
          mcLeptonThetaCC0PiEff->Fill(isEventPassed, mcLeptonTheta);
          mcLeptonPhiCC0PiEff->Fill(isEventPassed, mcLeptonPhi);
        }
      }
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

}

void EvaluatePerformance::endJob()
{


}

DEFINE_ART_MODULE(EvaluatePerformance)
