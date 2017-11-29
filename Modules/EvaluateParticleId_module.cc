////////////////////////////////////////////////////////////////////////
// Class:       EvaluateParticleId
// Plugin Type: analyzer (art v2_05_00)
// File:        EvaluateParticleId_module.cc
//
// Generated at Mon Nov 13 20:18:00 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// base includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

// local includes
#include "uboone/NuMuSelection/Algos/particleIdUtility.h"

class EvaluateParticleId;


class EvaluateParticleId : public art::EDAnalyzer {
public:
  explicit EvaluateParticleId(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EvaluateParticleId(EvaluateParticleId const &) = delete;
  EvaluateParticleId(EvaluateParticleId &&) = delete;
  EvaluateParticleId & operator = (EvaluateParticleId const &) = delete;
  EvaluateParticleId & operator = (EvaluateParticleId &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  pidutil::particleIdUtility pidutils;

  art::ServiceHandle< art::TFileService > tfs;

  // histograms
  TH1D* hAveragedQdX;
  TH2D* hAveragedQdXLength;
  TH2D* hAveragedQdXTheta;

  std::pair<double, double> averagedQdX;

};


EvaluateParticleId::EvaluateParticleId(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void EvaluateParticleId::analyze(art::Event const & e)
{

  // get event information
  int fEvent = e.event();
  int fRun = e.run();
  int fSubRun = e.subRun();

  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel("pandoraNu", trackHandle);

  art::Handle< std::vector< raw::RawDigit > > rawDHandle;
  e.getByLabel("wcNoiseFilter", rawDHandle);
  std::vector< art::Ptr< raw::RawDigit > > rawDVec;
  art::fill_ptr_vector(rawDVec, rawDHandle);

  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, "pandoraNu");

  if (trackHandle.product()->size() != 1) return;

  for (auto const& track : (*trackHandle)){

    // get associated hits
    std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(track.ID());

    averagedQdX = pidutils.getAveragedQdX(track, hits, rawDVec, false); 
    
    double averagedQdXMean = averagedQdX.first;
    //TH1D* averagedQdXhisto = averagedQdX.second;

    if (averagedQdXMean < 500)
      std::cout << fRun << "." << fSubRun << "." << fEvent << std::endl;

    hAveragedQdX->Fill(averagedQdXMean);
    hAveragedQdXLength->Fill(averagedQdXMean, track.Length());
    hAveragedQdXTheta->Fill(averagedQdXMean, track.Theta());

    TFile* chargeDistributions = new TFile("chargeDistros.root", "UPDATE");
    chargeDistributions->cd();
    //averagedQdXhisto->Write();
    //averagedQdXhisto->Delete();
    chargeDistributions->Close();

  }

}

void EvaluateParticleId::beginJob()
{

  hAveragedQdX = tfs->make<TH1D>("hAveragedQdX", ";average dQdX;", 100, 0, 2000);
  hAveragedQdXLength = tfs->make<TH2D>("hAveragedQdXLength"," ;average hit charge (ADC x Ticks); Length (cm)", 200, 0, 10000, 100, 0, 200);
  hAveragedQdXTheta = tfs->make<TH2D>("hAveragedQdXTheta", ";average hit charge (ADC x Ticks); Theta (degrees)", 200, 0, 10000, 100, 0, 3.15);

}

void EvaluateParticleId::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(EvaluateParticleId)
