////////////////////////////////////////////////////////////////////////
// Class:       SampleInformation
// Plugin Type: analyzer (art v2_05_00)
// File:        SampleInformation_module.cc
//
// Generated at Thu Oct 12 12:22:45 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// basic includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// other art includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"


class SampleInformation : public art::EDAnalyzer {
public:
  explicit SampleInformation(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SampleInformation(SampleInformation const &) = delete;
  SampleInformation(SampleInformation &&) = delete;
  SampleInformation & operator = (SampleInformation const &) = delete;
  SampleInformation & operator = (SampleInformation &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // services
  art::ServiceHandle< art::TFileService > tfs;

  // file
  TFile *file;

  // bools
  bool isPrintEventInfo = true;
  bool isPrintSummaryInfo = true;

  // mctruth origin information
  TH1D* mcTruthOrigin;
  std::vector<int> mcTruthOriginVector = {0,0,0,0,0};

};


SampleInformation::SampleInformation(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void SampleInformation::beginJob()
{

  file = tfs->make<TFile>("SampleInfo.root", "RECREATE");
  
  // mctruth origin information
  mcTruthOrigin = tfs->make<TH1D>("mcTruthOrigin", ";mcTruthOrigin;", 5, 0, 5);

}


void SampleInformation::analyze(art::Event const & e)
{

  if (e.isRealData()) return;

  // get MCTruth Information
  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
  e.getByLabel("generator", mcTruthHandle);

  int run = (int)e.run();
  int subRun = (int)e.subRun();
  int event = (int)e.event();

  if (isPrintEventInfo == true)
    std::cout << "------- RUN.SUBRUN.EVENT " 
      << run << "." 
      << subRun << "." 
      << event << " -------" << std::endl;

  for ( auto const& thisMcTruth : (*mcTruthHandle)){

    simb::Origin_t mctOrigin = thisMcTruth.Origin();
    if (mctOrigin == simb::kBeamNeutrino)
      mcTruthOriginVector.at(0)++;
    else if (mctOrigin == simb::kCosmicRay)
      mcTruthOriginVector.at(1)++;
    else if (mctOrigin == simb::kSuperNovaNeutrino)
      mcTruthOriginVector.at(2)++;
    else if (mctOrigin == simb::kSingleParticle)
      mcTruthOriginVector.at(3)++;
    else if (mctOrigin == simb::kUnknown)
      mcTruthOriginVector.at(4)++;
  
    if (isPrintEventInfo == true){
      
      std::cout << "MCTruth Origin: " << mctOrigin << std::endl; 

    }

  }

}
void SampleInformation::endJob()
{

  TFile& file = tfs->file();
  file.cd();

  for (int i = 0; i < mcTruthOrigin->GetNbinsX(); i++)
    mcTruthOrigin->SetBinContent(i, mcTruthOriginVector.at(i));

  if (isPrintSummaryInfo == true){

    std::cout << "\n\n\n\n>>>>>>> SUMMARY INFORMATION <<<<<<<" << std::endl;
    std::cout << "Sample contains truth information inheriting from: "
      << "\n\tBeam Neutrinos     : " << mcTruthOriginVector.at(0)
      << "\n\tCosmic Rays        : " << mcTruthOriginVector.at(1)
      << "\n\tSupernova Neutrinos: " << mcTruthOriginVector.at(2)
      << "\n\tSingle Particles   : " << mcTruthOriginVector.at(3)
      << "\n\tUnkwnown           : " << mcTruthOriginVector.at(4)
      << std::endl;

  }

  // write files
  mcTruthOrigin->Write();

}

DEFINE_ART_MODULE(SampleInformation)
