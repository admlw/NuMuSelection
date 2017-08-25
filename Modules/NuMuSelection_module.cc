////////////////////////////////////////////////////////////////////////
// Class:       NuMuSelection
// Plugin Type: analyzer (art v2_05_00)
// File:        NuMuSelection_module.cc
//
// Generated at Thu Aug 24 21:55:04 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TTree.h"

class NuMuSelection;


class NuMuSelection : public art::EDAnalyzer {
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
        void analyze(art::Event const & e) override;

        // Selected optional functions.
        void reconfigure(fhicl::ParameterSet const & p) override;
        void beginJob() override;
        void endJob() override;

    private:

        // fcl pars
        bool isData;
        bool isGetTruthTree;
        std::string fTruthLabel;
        std::string fParticleLabel;

        bool isMc;

        int fRun;
        int fSubRun;
        int fEvent;
        
        TTree *truthTree;
};


NuMuSelection::NuMuSelection(fhicl::ParameterSet const & p)
    :
        EDAnalyzer(p)  // ,
        // More initializers here.
{}

void NuMuSelection::reconfigure(fhicl::ParameterSet const & p)
{

    isData = p.get< bool > ("IsData");
    isGetTruthTree = p.get< bool > ("ProduceTruthTree");
    fTruthLabel = p.get< std::string > ("TruthLabel");
    fParticleLabel = p.get< std::string > ("MCParticleLabel");

}

void NuMuSelection::beginJob()
{

    art::ServiceHandle< art::TFileService > tfs;

    truthTree = tfs->make<TTree>("truthTree", "truthTree");

}

void NuMuSelection::analyze(art::Event const & e)
{
    fRun = e.run();
    fSubRun = e.subRun();
    fEvent = e.event();

    isData = e.isRealData();
    isMc = !isData;
    if (isData)
        std::cout << ">| Running over PHYSICS DATA" << std::endl;
    else
        std::cout << ">| Running over SIMULATED DATA" << std::endl;

    if (isGetTruthTree){

        std::cout << "\n>| Producing Truth Tree" << std::endl;

        art::Handle< std::vector<simb::MCTruth> > truthHandle; // you can't handle the truth
        e.getByLabel(fTruthLabel, truthHandle);
        
        for (auto const& truth : *(truthHandle)){

            auto const& trueOrigin = truth.Origin(); // simb::kBeamNeutrino or simb::kCosmicRay

            if (trueOrigin == simb::kBeamNeutrino){

                auto const& nu = truth.GetNeutrino();
                int nuPdg = nu.Nu().PdgCode(); 
                int nuCcNc = nu.CCNC();

                std::cout << ">>|------------- " << std::endl;
                std::cout << ">>| Found true neutrino with pdg code " << nuPdg << std::endl;
                std::cout << ">>| Interaction type? (CC = 0, NC = 1)" << nuCcNc << std::endl;

            }

        }


    }

    truthTree->Fill();
}

void NuMuSelection::endJob()
{
    // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NuMuSelection)
