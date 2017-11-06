////////////////////////////////////////////////////////////////////////
// Class:       ProduceEfficiencies
// Plugin Type: analyzer (art v2_05_00)
// File:        ProduceEfficiencies_module.cc
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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT Includes
#include "TTree.h"
#include "TEfficiency.h"

// UBXSec Includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"

// FMWK Includes

class ProduceEfficiencies;


class ProduceEfficiencies : public art::EDAnalyzer {
  public:
    explicit ProduceEfficiencies(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ProduceEfficiencies(ProduceEfficiencies const &) = delete;
    ProduceEfficiencies(ProduceEfficiencies &&) = delete;
    ProduceEfficiencies & operator = (ProduceEfficiencies const &) = delete;
    ProduceEfficiencies & operator = (ProduceEfficiencies &&) = delete;

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
    bool isCC0pi = true;
    double mcNuEnergy;
    int mcNuCCNC;

    double vx, vy, vz;

    // Efficiency histograms
    TEfficiency* mcNuEnergyEff;
    TEfficiency* mcNuEnergyCC0PiEff;
};


ProduceEfficiencies::ProduceEfficiencies(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
{


}

void ProduceEfficiencies::analyze(art::Event const & e)
{

  // get MC neutrino
  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
  e.getByLabel("generator", mcTruthHandle);
  if (!mcTruthHandle.isValid()) return;
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);

  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel("UBXSec", selectionHandle);
  if (!selectionHandle.isValid()){

    std::cout << " >> SelectionResult product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found."
      << std::endl;
    throw std::exception(); 

  }

  if (mcTruthVec.size() == 0){ 
    std::cout << " >> No Neutrinos " << std::endl;
    return;
  }

  art::Ptr<simb::MCTruth> mcTruth = mcTruthVec.at(0);
  const simb::MCNeutrino& mcNu = mcTruth->GetNeutrino();
  const simb::MCParticle& mcNuP = mcNu.Nu();

  vx = mcNuP.Vx();
  vy = mcNuP.Vy();
  vz = mcNuP.Vz();
  mcNuCCNC   = mcNu.CCNC();
  
  mcNuEnergy = mcNuP.E();

  // get true numu CC interactions in FV
  if (mcTruth->Origin() != simb::kBeamNeutrino){

    std::cout << " Not a beam neutrino! " << std::endl; 
    return;

  }

  if (){

    std::cout << "interaction vertex " << vx << ", " << vy << ", " << vz
      << "is out of TPC! " << std::endl;
    return;
  
  }

  if (mcNuCCNC != 0){
  
    std::cout << " >> Event is NC " << std::endl;
    return;
  
  }
  if (mcNuP.PdgCode() != 14){

    std::cout << " >> not a nu_mu interaction!" << std::endl;
    return;

  }

  // check if selected and fill efficiency histogram
  for (auto const& selectionStatus : (*selectionHandle)){

    std::cout << " >> Found a true CC Event " << std::endl;
    std::cout << " >> Printing daughter particles: " << std::endl; 
   
    int nParticles = mcTruth->NParticles();
    for(int i = 0; i < nParticles; i++){
      
        const simb::MCParticle& particle = mcTruth->GetParticle(i);
        if (particle.Process() != "primary" || particle.StatusCode() != 1) continue;    
        
        if (std::abs(particle.PdgCode()) == 211 || std::abs(particle.PdgCode()) == 111){

          isCC0pi = false;

        }
        else isCC0pi = true;

    }

    if (selectionStatus.GetSelectionStatus()){

      // event passes
      isEventPassed = true;

      mcNuEnergyEff->Fill(isEventPassed, mcNuEnergy);
      std::cout << " -->> Selected." << std::endl;

      if (isCC0pi == true)
        mcNuEnergyCC0PiEff->Fill(isEventPassed, mcNuEnergy);

    }
    else {
      // event does not pass
      isEventPassed = false;

      mcNuEnergyEff->Fill(isEventPassed, mcNuEnergy);
      std::cout << " -->>Not Selected." << std::endl;
     
      if (isCC0pi == true)
        mcNuEnergyCC0PiEff->Fill(isEventPassed, mcNuEnergy);

    }

  }

}

void ProduceEfficiencies::beginJob()
{
  // Implementation of optional member function here.
  selectionEfficiency = tfs->make<TTree>("selectionEfficiency", "selectionEfficiency");
  mcNuEnergyEff = tfs->make<TEfficiency>("mcNuEnergyEff", ";E_{#nu}^{true}; #epsilon", 15, 0, 3);
  mcNuEnergyCC0PiEff = tfs->make<TEfficiency>("mcNuEnergyCC0PiEff", ";E_{#nu}^{true}; #epsilon", 15, 0, 3);

}

void ProduceEfficiencies::endJob()
{


}

DEFINE_ART_MODULE(ProduceEfficiencies)
