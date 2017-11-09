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
    bool isEventPassed = false;
    bool isCC0pi = true;
    bool isBeamNeutrino;
    bool isRecoInTPC;
    bool isCC;
    bool isMuonNeutrino;
    bool isInFV;
    int mcNuCCNC;

    double vx, vy, vz;

    double mcNuEnergy;
    double mcLeptonMom;
    double mcLeptonTheta;
    double mcLeptonPhi;

    // Efficiency histograms
    TEfficiency* mcNuEnergyEff;
    TEfficiency* mcLeptonMomEff;
    TEfficiency* mcLeptonThetaEff;
    TEfficiency* mcLeptonPhiEff;

    TEfficiency* mcNuEnergyCC0PiEff;
    TEfficiency* mcLeptonMomCC0PiEff;
    TEfficiency* mcLeptonThetaCC0PiEff;
    TEfficiency* mcLeptonPhiCC0PiEff;

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


  // is true nu interaction in FV
  if (!fiducialVolume.InFV(vx, vy, vz)){

    std::cout << "true interaction vertex " << vx << ", " << vy << ", " << vz
      << "is out of TPC! " << std::endl;
    isInFV = false;

  }
  else isInFV = true;

  // ---------------------
  // Reco-level cuts
  // ---------------------
  std::cout << selectedTpcObjects.at(0).size() << std::endl;
  if (selectedTpcObjects.at(0).size() != 0){
    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const recob::Vertex& selectedVertex = selectedTpcObject->GetVertex();
    const ubana::TPCObjectOrigin& selectedOrigin = selectedTpcObject->GetOrigin();

    // get reconstructed vertex position
    // if true vertex is in the fiducial volume but reconstructed is out, this goes in the denominator
    double xyz[3] = {0.0,0.0,0.0};
    selectedVertex.XYZ(xyz);
    if (fiducialVolume.InFV(xyz[0], xyz[1], xyz[2]) == true){

      std::cout << "Reconstruced object is in TPC" << std::endl;
      isRecoInTPC = true;

    }
    else{
      std::cout << "Reconstructed object is out of TPC!" << std::endl;
      isRecoInTPC = false;
    }


    // get selected TPCObject  origin
    if ( selectedOrigin == ubana::kBeamNeutrino){

      std::cout << "Found a beam neutrino! " << std::endl; 
      isBeamNeutrino = true;

    }
    else{

      std::cout << "Not a beam neutrino!" << std::endl;
      isBeamNeutrino = false;

    }
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

      if (selectionStatus.GetSelectionStatus() == true && isBeamNeutrino == true && isRecoInTPC == true){

        // event passes
        isEventPassed = true;

        mcNuEnergyEff->Fill(isEventPassed, mcNuEnergy);
        mcLeptonMomEff->Fill(isEventPassed, mcLeptonMom);
        mcLeptonThetaEff->Fill(isEventPassed, mcLeptonTheta);
        mcLeptonPhiEff->Fill(isEventPassed, mcLeptonPhi);
        std::cout << " -->> Selected." << std::endl;

        if (isCC0pi == true){
          mcNuEnergyCC0PiEff->Fill(isEventPassed, mcNuEnergy);
          mcLeptonMomCC0PiEff->Fill(isEventPassed, mcLeptonMom);
          mcLeptonThetaCC0PiEff->Fill(isEventPassed, mcLeptonTheta);
          mcLeptonPhiCC0PiEff->Fill(isEventPassed, mcLeptonPhi);
        }
      }
      else {
        // event does not pass
        isEventPassed = false;

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
}

void EvaluatePerformance::beginJob()
{
  // Implementation of optional member function here.
  selectionEfficiency = tfs->make<TTree>("selectionEfficiency", "selectionEfficiency");
  mcNuEnergyEff = tfs->make<TEfficiency>("mcNuEnergyEff", ";E_{#nu}^{true}; #epsilon", 15, 0, 3);
  mcLeptonMomEff = tfs->make<TEfficiency>("mcLeptonMomEff", ";P_{l}^{true}; #epsilon", 15, 0, 2);
  mcLeptonThetaEff = tfs->make<TEfficiency>("mcLeptonThetaEff", ";#theta_{l}^{true}; #epsilon;", 10, -1, 1);
  mcLeptonPhiEff = tfs->make<TEfficiency>("mcLeptonPhiEff", ";#phi_{l}^{true}; #epsilon", 15, -3, 3);
  mcNuEnergyCC0PiEff = tfs->make<TEfficiency>("mcNuEnergyCC0PiEff", ";E_{#nu}^{true}; #epsilon", 15, 0, 3);
  mcLeptonMomCC0PiEff = tfs->make<TEfficiency>("mcLeptonMomCC0PiEff", ";P_{l}^{true}; #epsilon", 15, 0, 2);
  mcLeptonThetaCC0PiEff = tfs->make<TEfficiency>("mcLeptonThetaCC0PiEff", ";#theta_{l}^{true}; #epsilon;", 10, -1, 1);
  mcLeptonPhiCC0PiEff = tfs->make<TEfficiency>("mcLeptonPhiCC0PiEff", ";#phi_{l}^{true};#epsilon", 15, -3, 3);

}

void EvaluatePerformance::endJob()
{


}

DEFINE_ART_MODULE(EvaluatePerformance)
