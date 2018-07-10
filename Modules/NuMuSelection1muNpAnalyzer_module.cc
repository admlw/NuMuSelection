////////////////////////////////////////////////////////////////////////
// Class:       NuMuSelection1muNpAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        NuMuSelection1muNpAnalyzer_module.cc
//
// Generated at Mon Jul  9 10:35:20 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// This is a rewriting of the NuMuSelection producer module as an analyzer 
// module. 
//
// With the introduction of the new PID algorithms it became clear that the 
// ability to change the PID cuts on the fly would be extremely useful without
// having the full weight of LArSoft.
//
// Hindsight is 20/20.

// Base includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ART includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"        
#include "lardataobj/RecoBase/Hit.h"           
#include "lardata/Utilities/AssociationUtil.h" 
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h" 

// ROOT includes
#include "TTree.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"        
#include "uboone/UBXSec/DataTypes/TPCObject.h"              
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"        

// EventWeight includes                                     
#include "uboone/EventWeight/EventWeightTreeUtility.h"      

// local Includes                                           
#include "uboone/NuMuSelection/Algos/particleIdUtility.h"   


class NuMuSelection1muNpAnalyzer;


class NuMuSelection1muNpAnalyzer : public art::EDAnalyzer {
  public:
    explicit NuMuSelection1muNpAnalyzer(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    NuMuSelection1muNpAnalyzer(NuMuSelection1muNpAnalyzer const &) = delete;
    NuMuSelection1muNpAnalyzer(NuMuSelection1muNpAnalyzer &&) = delete;
    NuMuSelection1muNpAnalyzer & operator = (NuMuSelection1muNpAnalyzer const &) = delete;
    NuMuSelection1muNpAnalyzer & operator = (NuMuSelection1muNpAnalyzer &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endSubRun(art::SubRun const & sr) override;

    // User Functions

    /**
     * Initialises all tree variables for ana_tree
     */
    void initialiseAnalysisTree(TTree*, bool fIsData);

    /**
     * Resizes vectors to size 0 at the start of an event
     */
    void resizeVectors();

    std::pair< const simb::MCParticle*, double > GetAssociatedMCParticle(std::vector< art::Ptr<recob::Hit> > hits, art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit);

  private:

    // initialise services
    art::ServiceHandle< art::TFileService > tfs;

    // initialise classes
    ubana::SelectionResult selResult;
    ubana::FiducialVolume fiducialVolume;
    uboone::EWTreeUtil ewutil;

    // fhicl parameters
    std::string fSelectionLabel;
    std::string fHitLabel;
    std::string fTrackLabel;
    std::string fPIDLabel;
    std::string fSelectionTPCObjAssn;
    std::string fHitTrackAssn;
    std::string fHitTruthAssn;
    bool fIsData;

    // tree and tree variables
    TTree *ana_tree;
    int run;
    int subrun;
    int event;
    bool isData;
    bool isSimulation;
    bool isUBXSecSelected;
    bool isBeamNeutrino;
    bool isMixed;
    bool isCosmic;
    int nSelectedTracks;
    int nSelectedShowers;
    std::vector<double>* noBragg_fwd_mip = nullptr;
    std::vector<double>* bragg_fwd_mu    = nullptr;
    std::vector<double>* bragg_fwd_p     = nullptr;
    std::vector<double>* bragg_fwd_pi    = nullptr;
    std::vector<double>* bragg_fwd_k     = nullptr;
    std::vector<double>* bragg_bwd_mu    = nullptr;
    std::vector<double>* bragg_bwd_p     = nullptr;
    std::vector<double>* bragg_bwd_pi    = nullptr;
    std::vector<double>* bragg_bwd_k     = nullptr;
    std::vector<double>* track_length    = nullptr;
    std::vector<double>* track_theta     = nullptr;
    std::vector<double>* track_costheta  = nullptr;
    std::vector<double>* track_phi       = nullptr;
    std::vector<double>* track_startx    = nullptr;
    std::vector<double>* track_starty    = nullptr;
    std::vector<double>* track_startz    = nullptr;
    std::vector<double>* track_endx      = nullptr;
    std::vector<double>* track_endy      = nullptr;
    std::vector<double>* track_endz      = nullptr;
    double true_purity;
    double true_PDG;

};


NuMuSelection1muNpAnalyzer::NuMuSelection1muNpAnalyzer(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
{

  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");
  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolumeSettings");

  fSelectionLabel      = p_labels.get<std::string>("SelectionLabel", "UBXSec");
  fHitLabel            = p_labels.get<std::string>("HitLabel", "pandoraCosmicHitRemoval::UBXSec");
  fTrackLabel          = p_labels.get<std::string>("TrackLabel", "pandoraNu::UBXSec");
  fPIDLabel            = p_labels.get<std::string>("ParticleIdLabel");

  fSelectionTPCObjAssn = p_labels.get<std::string>("SelectionTPCObjAssn", "UBXSec");
  fHitTrackAssn        = p_labels.get<std::string>("hitTrackAssn", "pandoraNu::UBXSec");
  fHitTruthAssn        = p_labels.get<std::string>("hitTruthAssn", "pandoraCosmicHitRemoval::UBXSec");

  fIsData = p.get<bool>("IsData", "false");

}

void NuMuSelection1muNpAnalyzer::analyze(art::Event const & e)
{

  run    = e.run();
  subrun = e.subRun();
  event  = e.event();
  isData = e.isRealData();
  isSimulation = !isData;

  // get selection information for this event
  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel(fSelectionLabel, selectionHandle);

  if (!selectionHandle.isValid()){

    std::cout << "[NuMuSelection] ubana::SelectionResult product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "[NuMuSelection] ubana::SelectionResult product not found." << std::endl;
    throw std::exception();

  }

  art::FindManyP<ubana::TPCObject> selectedTPCObjects(selectionHandle, e, fSelectionTPCObjAssn);
  art::Ptr<ubana::TPCObject> selectedTPCObject;

  // other handles to information we're going to need
  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);

  art::Handle< std::vector< recob::Hit > > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  // the assns we need
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, fHitTrackAssn);
  art::FindManyP< anab::ParticleID > pidFromTrack(trackHandle, e, fPIDLabel);


  // finally, get  MCTruth, MCFlux, and GTruth information for reweighting                      

  art::Handle< std::vector< simb::MCTruth > > mcTruthHandle;                          
  std::vector< art::Ptr<simb::MCTruth> > mcTruthVec;                                  

  art::Handle< std::vector< simb::MCFlux > > mcFluxHandle;                            
  std::vector< art::Ptr<simb::MCFlux> > mcFluxVec;                                    

  art::Handle< std::vector< simb::GTruth > > gTruthHandle;                            
  std::vector< art::Ptr<simb::GTruth> > gTruthVec;                                    


  if (isSimulation){                                                                     
    e.getByLabel("generator", mcTruthHandle);                                         
    if (!mcTruthHandle.isValid()) return;                                             
    art::fill_ptr_vector(mcTruthVec, mcTruthHandle);                                  
    if (mcTruthVec.size() == 0){                                                      
      std::cout << "\n[NuMuSelection] No MCTruth Information" << std::endl;                 
      return;                                                                         
    }                                                                                 

    e.getByLabel("generator", gTruthHandle);                                          
    if (!gTruthHandle.isValid()) return;                                              
    art::fill_ptr_vector(gTruthVec, gTruthHandle);                                    
    if (gTruthVec.size() == 0){                                                       
      std::cout << "\n[NuMuSelection] No GTruth Information" << std::endl;                  
      return;                                                                         
    }                                                                                 

    e.getByLabel("generator", mcFluxHandle);                                          
    if (!mcFluxHandle.isValid()) return;                                              
    art::fill_ptr_vector(mcFluxVec, mcFluxHandle);                                    
    if (mcFluxVec.size() == 0){                                                       
      std::cout << "\n[NuMuSelection] No MCFlux Information" << std::endl;                  
      return;                                                                         
    }                                                                                 
  }                                                                                   

  art::Ptr<simb::MCFlux>  mcFlux;
  art::Ptr<simb::GTruth>  gTruth;
  art::Ptr<simb::MCTruth> mcTruth;

  if (isSimulation){
    mcFlux  = mcFluxVec.at(0);
    gTruth  = gTruthVec.at(0);
    mcTruth = mcTruthVec.at(0);
  }

  std::cout << "[NuMuSelection] --- run.subrun.event: " 
    << run << "." << subrun << "." << event << std::endl;

  // initialise variables which need scope
  recob::Vertex selectedVertex;
  nSelectedShowers = 0;
  nSelectedTracks = 0;
  resizeVectors();

  int nSelectedTPCObjects = selectedTPCObjects.at(0).size();

  // we expect a single TPC object
  if (nSelectedTPCObjects == 1){
    std::cout << "[NuMuSelection] Event is UBXSec Selected." << std::endl;
    isUBXSecSelected = true;

    // get selected TPC object and associated information
    selectedTPCObject = selectedTPCObjects.at(0).at(0);
    const std::vector<recob::Track>& selectedTracks = selectedTPCObject->GetTracks();
    selectedVertex   = selectedTPCObject->GetVertex();
    nSelectedShowers = (int)selectedTPCObject->GetNShowers();
    nSelectedTracks  = (int)selectedTPCObject->GetNTracks();
    const ubana::TPCObjectOrigin& selectedOrigin = selectedTPCObject->GetOrigin();

    if (isSimulation){

      isBeamNeutrino = false;
      isMixed        = false;
      isCosmic       = false;

      if (selectedOrigin == ubana::kBeamNeutrino) isBeamNeutrino = true;
      else isBeamNeutrino = false;

      if (selectedOrigin == ubana::kMixed) isMixed = true;
      else isMixed = false;

      if (selectedOrigin == ubana::kCosmicRay) isCosmic = true;
      else isCosmic = false;

    }

    // loop tracks in TPC object and get information
    for (size_t i = 0; i < selectedTracks.size(); i++){

      const recob::Track& track = selectedTracks.at(i);

      std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(track.ID());
      std::vector< art::Ptr< anab::ParticleID > > pids = pidFromTrack.at(track.ID());

      if (pids.size() == 0){

        std::cout << "[NuMuSelection] No track<->PID association found for track with ID " << track.ID() << ". Throwing exception. " << std::endl;
        throw std::exception();

      }

      // get PID information
      std::vector< anab::sParticleIDAlgScores > algScoresVec = pids.at(0)->ParticleIDAlgScores();

      for (size_t i_alg = 0; i_alg < algScoresVec.size(); i_alg++){

        anab::sParticleIDAlgScores algScore = algScoresVec.at(i_alg);

        if(algScore.fAlgName == "BraggPeakLLH" && algScore.fPlaneID.Plane == 2){                                   
          if (anab::kVariableType(algScore.fVariableType) == anab::kLogL_fwd){     

            if (algScore.fAssumedPdg == 13  ) bragg_fwd_mu->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 2212) bragg_fwd_p->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 211 ) bragg_fwd_pi->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 321 ) bragg_fwd_k->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 0   ) noBragg_fwd_mip->push_back(algScore.fValue);   

          }                                                                        
          else if (anab::kVariableType(algScore.fVariableType) == anab::kLogL_bwd){

            if (algScore.fAssumedPdg == 13  ) bragg_bwd_mu->push_back(algScore.fValue);      
            if (algScore.fAssumedPdg == 2212) bragg_bwd_p->push_back(algScore.fValue);      
            if (algScore.fAssumedPdg == 211 ) bragg_bwd_pi->push_back(algScore.fValue);      
            if (algScore.fAssumedPdg == 321 ) bragg_bwd_k->push_back(algScore.fValue);      

          }                                                                        

        }

      }

      track_length->push_back(track.Length());
      track_theta->push_back(track.Theta());
      track_costheta->push_back(std::cos(track.Theta()));
      track_phi->push_back(track.Phi());
      track_startx->push_back(track.Start().X());
      track_starty->push_back(track.Start().Y());
      track_startz->push_back(track.Start().Z());
      track_endx->push_back(track.End().X());
      track_endy->push_back(track.End().Y());
      track_endz->push_back(track.End().Z());

      if (isSimulation){
        art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit(hitHandle, e, fHitTruthAssn);
        std::pair< const simb::MCParticle*, double > mcpInformation = GetAssociatedMCParticle(hits, particlesPerHit);

        const simb::MCParticle* mcpMatchedParticle = mcpInformation.first;
        true_purity = mcpInformation.second;
        
        true_PDG = mcpMatchedParticle->PdgCode();

      }
    }
  }
  else if (nSelectedTPCObjects == 0){
    // this is probably fine

    std::cout << "[NuMuSelection] Event has " << nSelectedTPCObjects << std::endl;
    isUBXSecSelected = false;
    return;

  }
  // this could indicate a problem
  else if (nSelectedTPCObjects > 1){

    std::cout << "[NuMuSelection] Event has " << nSelectedTPCObjects << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "[NuMuSelection] Event has " << nSelectedTPCObjects << std::endl;
    throw std::exception();

  }

  ana_tree->Fill();

}

void NuMuSelection1muNpAnalyzer::beginJob()
{

  ana_tree = tfs->make<TTree>("analysis_tree", "analysis tree");
  initialiseAnalysisTree(ana_tree, fIsData);

}

void NuMuSelection1muNpAnalyzer::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void NuMuSelection1muNpAnalyzer::initialiseAnalysisTree(TTree *tree, bool fIsData){

  tree->Branch("run"              , &run);
  tree->Branch("subrun"           , &subrun);
  tree->Branch("event"            , &event);
  tree->Branch("isData"           , &isData);
  tree->Branch("isSimulation"     , &isSimulation);
  tree->Branch("isUBXSecSelected" , &isUBXSecSelected);
  tree->Branch("noBragg_fwd_mip"  , &noBragg_fwd_mip);
  tree->Branch("bragg_fwd_mu"     , "std::vector<double>" , &bragg_fwd_mu);
  tree->Branch("bragg_fwd_p"      , "std::vector<double>" , &bragg_fwd_p );
  tree->Branch("bragg_fwd_pi"     , "std::vector<double>" , &bragg_fwd_pi);
  tree->Branch("bragg_fwd_k"      , "std::vector<double>" , &bragg_fwd_k );
  tree->Branch("bragg_bwd_mu"     , "std::vector<double>" , &bragg_bwd_mu);
  tree->Branch("bragg_bwd_p"      , "std::vector<double>" , &bragg_bwd_p );
  tree->Branch("bragg_bwd_pi"     , "std::vector<double>" , &bragg_bwd_pi);
  tree->Branch("bragg_bwd_k"      , "std::vector<double>" , &bragg_bwd_k );
  tree->Branch("track_length"     , "std::vector<double>" , &track_length   );
  tree->Branch("track_theta"      , "std::vector<double>" , &track_theta    );
  tree->Branch("track_costheta"   , "std::vector<double>" , &track_costheta );
  tree->Branch("track_phi"        , "std::vector<double>" , &track_phi      );
  tree->Branch("track_startx"     , "std::vector<double>" , &track_startx   );
  tree->Branch("track_starty"     , "std::vector<double>" , &track_starty   );
  tree->Branch("track_startz"     , "std::vector<double>" , &track_startz   );
  tree->Branch("track_endx"       , "std::vector<double>" , &track_endx     );
  tree->Branch("track_endy"       , "std::vector<double>" , &track_endy     );
  tree->Branch("track_endz"       , "std::vector<double>" , &track_endz     );
  tree->Branch("true_purity"      , &true_purity);
  tree->Branch("true_PDG"         , &true_PDG);

}

void NuMuSelection1muNpAnalyzer::resizeVectors(){

    noBragg_fwd_mip->resize(0);
    bragg_fwd_mu->resize(0);
    bragg_fwd_p->resize(0);
    bragg_fwd_pi->resize(0);
    bragg_fwd_k->resize(0);
    bragg_bwd_mu->resize(0);
    bragg_bwd_p->resize(0);
    bragg_bwd_pi->resize(0);
    bragg_bwd_k->resize(0);
    track_length->resize(0);
    track_theta->resize(0);
    track_costheta->resize(0);
    track_phi->resize(0);
    track_startx->resize(0);
    track_starty->resize(0);
    track_startz->resize(0);
    track_endx->resize(0);
    track_endy->resize(0);
    track_endz->resize(0);

}

std::pair< const simb::MCParticle*, double > NuMuSelection1muNpAnalyzer::GetAssociatedMCParticle(std::vector< art::Ptr< recob::Hit > > hits, art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit)
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

  std::pair< const simb::MCParticle* , double > returner;
  returner.first  = maxp_me;
  returner.second = maxe/tote;

  return returner;

}

DEFINE_ART_MODULE(NuMuSelection1muNpAnalyzer)
