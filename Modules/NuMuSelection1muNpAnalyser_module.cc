////////////////////////////////////////////////////////////////////////
// Class:       NuMuSelection1muNpAnalyser
// Plugin Type: analyzer (art v2_05_01)
// File:        NuMuSelection1muNpAnalyser_module.cc
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
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardata/Utilities/AssociationUtil.h" 
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h" 
#include "larcore/Geometry/Geometry.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// ROOT includes
#include "TTree.h"

// ubxsec includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"        
#include "uboone/UBXSec/DataTypes/TPCObject.h"              
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"        

// ParticleID includes
#include "uboone/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

// EventWeight includes                                     
#include "uboone/EventWeight/EventWeightTreeUtility.h"      

// local Includes                                           
#include "uboone/NuMuSelection/Algos/mapBuilderUtility.h"   


class NuMuSelection1muNpAnalyser;


class NuMuSelection1muNpAnalyser : public art::EDAnalyzer {
  public:
    explicit NuMuSelection1muNpAnalyser(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    NuMuSelection1muNpAnalyser(NuMuSelection1muNpAnalyser const &) = delete;
    NuMuSelection1muNpAnalyser(NuMuSelection1muNpAnalyser &&) = delete;
    NuMuSelection1muNpAnalyser & operator = (NuMuSelection1muNpAnalyser const &) = delete;
    NuMuSelection1muNpAnalyser & operator = (NuMuSelection1muNpAnalyser &&) = delete;

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
    void resizeVectors(bool isData);

    /**
     * puts dummy variables in variables for writing to tree
     */
    void emplaceDummyVars();

    std::pair< const simb::MCParticle*, double > GetAssociatedMCParticle(std::vector< art::Ptr<recob::Hit> > hits, art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit);

  private:

    // initialise services
    art::ServiceHandle< art::TFileService > tfs;
    art::ServiceHandle< geo::Geometry > geo;

    // initialise classes
    ubana::FiducialVolume _fiducialVolume;
    uboone::EWTreeUtil _ewutil;
    numusel::MapBuilderUtil _mbutil;
    trkf::TrackMomentumCalculator _trkmom;

    // fhicl parameters
    std::string fSelectionLabel;
    std::string fHitLabel;
    std::string fTrackLabel;
    std::string fMCSLabel;
    std::string fTrackSmearedCalorimetryAssn;
    std::string fTrackUnsmearedCalorimetryAssn;
    std::string fTrackParticleIdAssn;
    std::string fSmearedCalorimetryAssn;
    std::string fUnsmearedCalorimetryAssn;
    std::string fSelectionTPCObjAssn;
    std::string fHitTrackAssn;
    std::string fHitTruthAssn;
    std::string fTruthMcpAssn;
    bool fIsData;

    // subrun counting tree and variables
    TTree *sr_tree;
    int sr_run;
    int sr_subrun;
    double sr_begintime;
    double sr_endtime;
    double sr_pot;

    // tree and tree variables
    TTree *ana_tree;
    int run;
    int subrun;
    int event;
    bool isData;
    bool isSimulation;
    bool isUBXSecSelected;

    // -- ubxsec
    bool isBeamNeutrino = false;
    bool isMixed = false;
    bool isCosmic = false;
    bool isInFV = false;
    int nSelectedTracks = -999;
    int nSelectedShowers = -999;
    int nSelectedPfparticles = -999;

    // -- particleid
    std::vector<double>* noBragg_fwd_mip  = nullptr;
    std::vector<double>* bragg_fwd_mu     = nullptr;
    std::vector<double>* bragg_fwd_p      = nullptr;
    std::vector<double>* bragg_fwd_pi     = nullptr;
    std::vector<double>* bragg_fwd_k      = nullptr;
    std::vector<double>* bragg_bwd_mu     = nullptr;
    std::vector<double>* bragg_bwd_p      = nullptr;
    std::vector<double>* bragg_bwd_pi     = nullptr;
    std::vector<double>* bragg_bwd_k      = nullptr;

    std::vector<double>* track_dep_energy = nullptr;

    // -- reco track
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
    std::vector<double>* track_chi2      = nullptr;
    std::vector<int>* track_ndof         = nullptr;
    std::vector<int>* track_ntrajpoints  = nullptr;
    std::vector<bool>* track_isContained = nullptr;
    std::vector< std::vector<double> >* track_dedxperhit_smeared = nullptr;
    std::vector< std::vector<double> >* track_dedxperhit_unsmeared = nullptr;
    std::vector< std::vector<double> >* track_resrangeperhit = nullptr;

    // -- reco track mcs 
    // n.b. the particle hypothesis is always a muon here
    std::vector<double>* track_mcs_muassmp_fwd = nullptr;
    std::vector<double>* track_mcs_muassmp_bwd = nullptr;
    std::vector<double>* track_mcs_muassmp_energy_fwd = nullptr;
    std::vector<double>* track_mcs_muassmp_energy_bwd = nullptr;
    std::vector<double>* track_mcs_muassmp_fwd_uncertainty = nullptr;
    std::vector<double>* track_mcs_muassmp_bwd_uncertainty = nullptr;
    std::vector<int>* track_mcs_muassmp_particlehypothesis = nullptr;
    std::vector<double>* track_mcs_muassmp_fwd_loglikelihood = nullptr;
    std::vector<double>* track_mcs_muassmp_bwd_loglikelihood = nullptr;
    std::vector<double>* track_range_mom_muassumption = nullptr;
    std::vector<double>* track_range_mom_passumption = nullptr;
    std::vector<double>* track_range_energy_muassumption = nullptr;
    std::vector<double>* track_range_energy_passumption = nullptr;

    // -- reco vertex
    double vertex_x = -999;
    double vertex_y = -999;
    double vertex_z = -999;

    // -- truth information
    int true_truth_nparticles = 0;
    std::vector<int>* true_genie_pdg                         = nullptr;
    std::vector<int>* true_genie_statuscode                  = nullptr;
    std::vector<int>* true_genie_trackid                     = nullptr;
    std::vector<int>* true_genie_motherid                    = nullptr;
    std::vector<std::string>* true_genie_process             = nullptr;
    std::vector<std::string>* true_genie_endprocess          = nullptr;
    std::vector<int>* true_genie_numdaugters                 = nullptr;
    std::vector< std::vector<int> >* true_genie_daughterids  = nullptr;
    std::vector<double>* true_genie_startx                   = nullptr;
    std::vector<double>* true_genie_starty                   = nullptr;
    std::vector<double>* true_genie_startz                   = nullptr;
    std::vector<double>* true_genie_startt                   = nullptr; 
    std::vector<double>* true_genie_endx                     = nullptr;
    std::vector<double>* true_genie_endy                     = nullptr;
    std::vector<double>* true_genie_endz                     = nullptr;
    std::vector<double>* true_genie_endt                     = nullptr;
    std::vector<double>* true_genie_mass                     = nullptr;
    std::vector<double>* true_genie_starte                   = nullptr;
    std::vector<double>* true_genie_startp                   = nullptr;
    std::vector<double>* true_genie_startpx                  = nullptr;
    std::vector<double>* true_genie_startpy                  = nullptr;
    std::vector<double>* true_genie_startpz                  = nullptr;
    std::vector<double>* true_genie_ende                     = nullptr;
    std::vector<double>* true_genie_endpx                    = nullptr;
    std::vector<double>* true_genie_endpy                    = nullptr;
    std::vector<double>* true_genie_endpz                    = nullptr;
    std::vector<int>* true_genie_rescatter                   = nullptr;
    std::vector<int>* true_mcp_pdg                         = nullptr;
    std::vector<int>* true_mcp_statuscode                  = nullptr;
    std::vector<int>* true_mcp_trackid                     = nullptr;
    std::vector<int>* true_mcp_motherid                    = nullptr;
    std::vector<std::string>* true_mcp_process             = nullptr;
    std::vector<std::string>* true_mcp_endprocess          = nullptr;
    std::vector<int>* true_mcp_numdaugters                 = nullptr;
    std::vector< std::vector<int> >* true_mcp_daughterids  = nullptr;
    std::vector<double>* true_mcp_startx                   = nullptr;
    std::vector<double>* true_mcp_starty                   = nullptr;
    std::vector<double>* true_mcp_startz                   = nullptr;
    std::vector<double>* true_mcp_startt                   = nullptr; 
    std::vector<double>* true_mcp_endx                     = nullptr;
    std::vector<double>* true_mcp_endy                     = nullptr;
    std::vector<double>* true_mcp_endz                     = nullptr;
    std::vector<double>* true_mcp_endt                     = nullptr;
    std::vector<double>* true_mcp_mass                     = nullptr;
    std::vector<double>* true_mcp_starte                   = nullptr;
    std::vector<double>* true_mcp_startp                   = nullptr;
    std::vector<double>* true_mcp_startpx                  = nullptr;
    std::vector<double>* true_mcp_startpy                  = nullptr;
    std::vector<double>* true_mcp_startpz                  = nullptr;
    std::vector<double>* true_mcp_ende                     = nullptr;
    std::vector<double>* true_mcp_endpx                    = nullptr;
    std::vector<double>* true_mcp_endpy                    = nullptr;
    std::vector<double>* true_mcp_endpz                    = nullptr;
    std::vector<int>* true_mcp_rescatter                   = nullptr;
    int true_nu_ccnc = 0;
    int true_nu_mode = 0;
    int true_nu_interactiontype = 0;
    int true_nu_target = 0;
    int true_nu_hitnuc = 0;
    int true_nu_hitquark = 0;
    double true_nu_w = 0;
    double true_nu_x = 0;
    double true_nu_y = 0;
    double true_nu_qsqr = 0;
    double true_nu_pt = 0;
    double true_nu_theta = 0;

    // -- mcparticle
    std::vector<double>* true_match_purity                   = nullptr;
    std::vector<int>* true_match_pdg                         = nullptr;
    std::vector<int>* true_match_statuscode                  = nullptr;
    std::vector<int>* true_match_trackid                     = nullptr;
    std::vector<int>* true_match_motherid                    = nullptr;
    std::vector<std::string>* true_match_process             = nullptr;
    std::vector<std::string>* true_match_endprocess          = nullptr;
    std::vector<int>* true_match_numdaugters                 = nullptr;
    std::vector< std::vector<int> >* true_match_daughterids  = nullptr;
    std::vector<double>* true_match_startx                   = nullptr;
    std::vector<double>* true_match_starty                   = nullptr;
    std::vector<double>* true_match_startz                   = nullptr;
    std::vector<double>* true_match_startt                   = nullptr;
    std::vector<double>* true_match_endx                     = nullptr;
    std::vector<double>* true_match_endy                     = nullptr;
    std::vector<double>* true_match_endz                     = nullptr;
    std::vector<double>* true_match_endt                     = nullptr;
    std::vector<double>* true_match_mass                     = nullptr;
    std::vector<double>* true_match_starte                   = nullptr;
    std::vector<double>* true_match_startp                   = nullptr;
    std::vector<double>* true_match_startpx                  = nullptr;
    std::vector<double>* true_match_startpy                  = nullptr;
    std::vector<double>* true_match_startpz                  = nullptr;
    std::vector<double>* true_match_ende                     = nullptr;
    std::vector<double>* true_match_endpx                    = nullptr;
    std::vector<double>* true_match_endpy                    = nullptr;
    std::vector<double>* true_match_endpz                    = nullptr;
    std::vector<int>* true_match_rescatter                   = nullptr;


};


NuMuSelection1muNpAnalyser::NuMuSelection1muNpAnalyser(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
{

  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");
  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolumeSettings");

  fSelectionLabel      = p_labels.get<std::string>("SelectionLabel", "UBXSec");
  fHitLabel            = p_labels.get<std::string>("HitLabel", "pandoraCosmicHitRemoval::UBXSec");
  fTrackLabel          = p_labels.get<std::string>("TrackLabel", "pandoraNu::UBXSec");
  fMCSLabel            = p_labels.get<std::string>("MCSLabel", "pandoraNuMCSMu");

  fTrackSmearedCalorimetryAssn   = p_labels.get<std::string>("TrackSmearedCaloAssn", "simcalibration::numusel");
  fTrackUnsmearedCalorimetryAssn = p_labels.get<std::string>("TrackUnsmearedCaloAssn", "pidcalibration::numusel");
  fTrackParticleIdAssn    = p_labels.get<std::string>("TrackParticleIdAssn");
  fSelectionTPCObjAssn    = p_labels.get<std::string>("SelectionTPCObjAssn", "UBXSec");
  fHitTrackAssn           = p_labels.get<std::string>("HitTrackAssn", "pandoraNu::UBXSec");
  fHitTruthAssn           = p_labels.get<std::string>("HitTruthAssn", "pandoraCosmicHitRemoval::UBXSec");
  fTruthMcpAssn           = p_labels.get<std::string>("TruthMcpAssn", "largeant");

  fIsData = p.get<bool>("IsData", "false");

  // configure
  _fiducialVolume.Configure(p_fv,
      geo->DetHalfHeight(),
      2.*geo->DetHalfWidth(),
      geo->DetLength());

  _fiducialVolume.PrintConfig();

}

void NuMuSelection1muNpAnalyser::analyze(art::Event const & e)
{

  run    = e.run();
  subrun = e.subRun();
  event  = e.event();
  isData = e.isRealData();
  isSimulation = !isData;

  std::cout << "[NuMuSelection] --- run.subrun.event: " 
    << run << "." << subrun << "." << event << std::endl;

  // resize vectors here to avoid any over-running memory issues
  NuMuSelection1muNpAnalyser::resizeVectors(isData);

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
  std::vector< art::Ptr<recob::Track> > trackPtrVector;
  art::fill_ptr_vector(trackPtrVector, trackHandle);

  art::Handle< std::vector< recob::Hit > > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  art::Handle< std::vector< recob::MCSFitResult > > mcsHandle;
  e.getByLabel(fMCSLabel, mcsHandle);
  std::vector< art::Ptr<recob::MCSFitResult> > mcsFitPtrVector;
  art::fill_ptr_vector(mcsFitPtrVector, mcsHandle);

  // there is no association between tracks and the MCSResult stored in the file, 
  // so going to build a map between each track ID and an art::Ptr to the 
  // recob::MCSFitResult for quick finding
  std::map<int, art::Ptr<recob::MCSFitResult> > trackIdMcsFitMap = _mbutil.buildTrackIdMcsResultMap(trackPtrVector, mcsFitPtrVector);

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

  // the assns we need
  art::FindManyP< recob::Hit > hitsFromTrack(trackHandle, e, fHitTrackAssn);
  art::FindManyP< anab::ParticleID > pidFromTrack(trackHandle, e, fTrackParticleIdAssn);
  art::FindManyP< anab::Calorimetry > unsmearedCaloFromTrack(trackHandle, e, fTrackUnsmearedCalorimetryAssn);
  art::FindManyP< anab::Calorimetry > smearedCaloFromTrack(trackHandle, e, fTrackSmearedCalorimetryAssn);

  art::Ptr<simb::MCFlux>  mcFlux;
  art::Ptr<simb::GTruth>  gTruth;
  art::Ptr<simb::MCTruth> mcTruth;

  if (isSimulation){

    art::FindManyP< simb::MCParticle > mcparticleFromTruth(mcTruthHandle, e, fTruthMcpAssn);

    mcFlux  = mcFluxVec.at(0);
    gTruth  = gTruthVec.at(0);
    mcTruth = mcTruthVec.at(0);

    true_truth_nparticles = mcTruth->NParticles();
    const simb::MCNeutrino mcNeutrino = mcTruth->GetNeutrino();
    true_nu_ccnc = mcNeutrino.CCNC();
    true_nu_mode = mcNeutrino.Mode();
    true_nu_interactiontype = mcNeutrino.InteractionType();
    true_nu_target = mcNeutrino.Target();
    true_nu_hitnuc = mcNeutrino.HitNuc();
    true_nu_hitquark = mcNeutrino.HitQuark();
    true_nu_w = mcNeutrino.W();
    true_nu_x = mcNeutrino.X();
    true_nu_y = mcNeutrino.Y();
    true_nu_qsqr = mcNeutrino.QSqr();
    true_nu_pt = mcNeutrino.Pt();
    true_nu_theta = mcNeutrino.Theta();


    for (int i = 0; i < true_truth_nparticles; i++){

      simb::MCParticle mcParticle = mcTruth->GetParticle(i);

      true_genie_pdg->push_back(mcParticle.PdgCode());
      true_genie_statuscode->push_back(mcParticle.StatusCode());
      true_genie_trackid->push_back(mcParticle.TrackId());
      true_genie_motherid->push_back(mcParticle.Mother());
      true_genie_process->push_back(mcParticle.Process());
      true_genie_endprocess->push_back(mcParticle.EndProcess());
      true_genie_numdaugters->push_back(mcParticle.NumberDaughters());

      std::vector<int> daughterIds;
      daughterIds.resize(0);
      for (int j = 0; j < mcParticle.NumberDaughters(); j++){
        daughterIds.push_back(mcParticle.Daughter(j));
      }
      true_genie_daughterids->push_back(daughterIds);

      true_genie_startx->push_back(mcParticle.Vx());
      true_genie_starty->push_back(mcParticle.Vy());
      true_genie_startz->push_back(mcParticle.Vz());
      true_genie_startt->push_back(mcParticle.T());
      true_genie_endx->push_back(mcParticle.EndX());
      true_genie_endy->push_back(mcParticle.EndY());
      true_genie_endz->push_back(mcParticle.EndZ());
      true_genie_endt->push_back(mcParticle.EndT());
      true_genie_mass->push_back(mcParticle.Mass());
      true_genie_starte->push_back(mcParticle.E());
      true_genie_startp->push_back(mcParticle.P());
      true_genie_startpx->push_back(mcParticle.Px());
      true_genie_startpy->push_back(mcParticle.Py());
      true_genie_startpz->push_back(mcParticle.Pz());
      true_genie_ende->push_back(mcParticle.EndE());
      true_genie_endpx->push_back(mcParticle.EndPx());
      true_genie_endpy->push_back(mcParticle.EndPy());
      true_genie_endpz->push_back(mcParticle.EndPz());
      true_genie_rescatter->push_back(mcParticle.Rescatter());

    }

    // MCParticle.at(0) is always the neutrino 
    std::cout << "[NuMuSelection] True Neutrino PDG:  " << true_genie_pdg->at(0) << std::endl;
    std::cout << "[NuMuSelection] True Neutrino CCNC: " << true_nu_ccnc << std::endl;
    std::cout << "[NuMuSelection] True Neutrino Mode: " << true_nu_mode << std::endl;
    std::cout << "[NuMuSelection] True Neutrino Interaction Type: " << true_nu_interactiontype << std::endl;

    // now get the MCParticles from the MCTruth<->MCParticle Assn
    std::vector< art::Ptr<simb::MCParticle> > mcparticles = mcparticleFromTruth.at(mcTruth.key());

    for (size_t i = 0; i < mcparticles.size(); i++){

      art::Ptr< simb::MCParticle > mcParticle = mcparticles.at(i);

      true_mcp_pdg->push_back(mcParticle->PdgCode());
      true_mcp_statuscode->push_back(mcParticle->StatusCode());
      true_mcp_trackid->push_back(mcParticle->TrackId());
      true_mcp_motherid->push_back(mcParticle->Mother());
      true_mcp_process->push_back(mcParticle->Process());
      true_mcp_endprocess->push_back(mcParticle->EndProcess());
      true_mcp_numdaugters->push_back(mcParticle->NumberDaughters());

      std::vector<int> daughterIds;
      daughterIds.resize(0);
      for (int j = 0; j < mcParticle->NumberDaughters(); j++){
        daughterIds.push_back(mcParticle->Daughter(j));
      }
      true_mcp_daughterids->push_back(daughterIds);

      true_mcp_startx->push_back(mcParticle->Vx());
      true_mcp_starty->push_back(mcParticle->Vy());
      true_mcp_startz->push_back(mcParticle->Vz());
      true_mcp_startt->push_back(mcParticle->T());
      true_mcp_endx->push_back(mcParticle->EndX());
      true_mcp_endy->push_back(mcParticle->EndY());
      true_mcp_endz->push_back(mcParticle->EndZ());
      true_mcp_endt->push_back(mcParticle->EndT());
      true_mcp_mass->push_back(mcParticle->Mass());
      true_mcp_starte->push_back(mcParticle->E());
      true_mcp_startp->push_back(mcParticle->P());
      true_mcp_startpx->push_back(mcParticle->Px());
      true_mcp_startpy->push_back(mcParticle->Py());
      true_mcp_startpz->push_back(mcParticle->Pz());
      true_mcp_ende->push_back(mcParticle->EndE());
      true_mcp_endpx->push_back(mcParticle->EndPx());
      true_mcp_endpy->push_back(mcParticle->EndPy());
      true_mcp_endpz->push_back(mcParticle->EndPz());
      true_mcp_rescatter->push_back(mcParticle->Rescatter());

    }

    isInFV = _fiducialVolume.InFV(true_genie_startx->at(0), true_genie_starty->at(0), true_genie_startz->at(0));
    std::cout << "[NuMuSelection] isInFV: " << isInFV << std::endl;


  }

  // initialise variables which need scope
  recob::Vertex selectedVertex;
  nSelectedShowers = 0;
  nSelectedTracks = 0;
  nSelectedPfparticles = 0;

  int nSelectedTPCObjects = selectedTPCObjects.at(0).size();

  // we expect a single TPC object
  if (nSelectedTPCObjects == 1){
    std::cout << "[NuMuSelection] Event is UBXSec Selected." << std::endl;
    isUBXSecSelected = true;

    // get selected TPC object and associated information
    selectedTPCObject = selectedTPCObjects.at(0).at(0);
    const std::vector<recob::Track>& selectedTracks = selectedTPCObject->GetTracks();
    selectedVertex   = selectedTPCObject->GetVertex();
    double xyz[3];
    selectedVertex.XYZ(xyz);
    vertex_x = xyz[0];
    vertex_y = xyz[1];
    vertex_z = xyz[2];

    const ubana::TPCObjectOrigin& selectedOrigin = selectedTPCObject->GetOrigin();
    
    selectedTPCObject->GetMultiplicity(nSelectedPfparticles,nSelectedTracks,nSelectedShowers);

    std::cout << "[NuMuSelection] Vertex position: " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << std::endl; 
    std::cout << "[NuMuSelection] Number of tracks: " << nSelectedTracks << std::endl;
    std::cout << "[NuMuSelection] Number of showers: " << nSelectedShowers << std::endl;

    if (isSimulation){

      _ewutil.WriteTree(e, mcFlux, mcTruth, gTruth);

      isBeamNeutrino = false;
      isMixed        = false;
      isCosmic       = false;

      if (selectedOrigin == ubana::kBeamNeutrino) isBeamNeutrino = true;
      else isBeamNeutrino = false;

      if (selectedOrigin == ubana::kMixed) isMixed = true;
      else isMixed = false;

      if (selectedOrigin == ubana::kCosmicRay) isCosmic = true;
      else isCosmic = false;

      std::cout << "[NuMuSelection] isBeamNeutrino: " << isBeamNeutrino << std::endl;
      std::cout << "[NuMuSelection] isMixed: " << isMixed << std::endl;
      std::cout << "[NuMuSelection] isCosmic: " << isCosmic << std::endl;

    }

    // loop tracks in TPC object and get information
    for (size_t i = 0; i < selectedTracks.size(); i++){

      const recob::Track& track = selectedTracks.at(i);

      std::vector< art::Ptr< recob::Hit > > hits = hitsFromTrack.at(track.ID());
      std::vector< art::Ptr< anab::ParticleID > > pids = pidFromTrack.at(track.ID());
      std::vector< art::Ptr< anab::Calorimetry > > unsmearedCalos = unsmearedCaloFromTrack.at(track.ID());
      std::vector< art::Ptr< anab::Calorimetry > > smearedCalos   = smearedCaloFromTrack.at(track.ID());

      // check whether track is contained in FV or not
      bool start_isContained = _fiducialVolume.InFV(track.Start().X(), track.Start().Y(), track.Start().Z());
      bool end_isContained   = _fiducialVolume.InFV(track.End().X(), track.End().Y(), track.End().Z());

      if (start_isContained && end_isContained)
        track_isContained->push_back(true);
      else track_isContained->push_back(false);

      // also get MCS fit result from MCS map we built earlier
      art::Ptr<recob::MCSFitResult> mcsFitResult;

      for (auto const& x : trackIdMcsFitMap){

        if (track.ID() == x.first){

          mcsFitResult = x.second;
          break;

        }

      }

      // Note MCS assumption is muon-like
      track_mcs_muassmp_fwd->push_back(mcsFitResult->fwdMomentum());
      track_mcs_muassmp_bwd->push_back(mcsFitResult->bwdMomentum());
      track_mcs_muassmp_energy_fwd->push_back((std::sqrt(std::pow(mcsFitResult->fwdMomentum(),2)+std::pow(105.7,2)) - 105.7)*1000.);
      track_mcs_muassmp_energy_bwd->push_back((std::sqrt(std::pow(mcsFitResult->bwdMomentum(),2)+std::pow(105.7,2)) - 105.7)*1000.);
      track_mcs_muassmp_fwd_uncertainty->push_back(mcsFitResult->fwdMomUncertainty());
      track_mcs_muassmp_bwd_uncertainty->push_back(mcsFitResult->bwdMomUncertainty());
      track_mcs_muassmp_fwd_loglikelihood->push_back(mcsFitResult->fwdLogLikelihood());
      track_mcs_muassmp_bwd_loglikelihood->push_back(mcsFitResult->bwdLogLikelihood());
      track_mcs_muassmp_particlehypothesis->push_back(mcsFitResult->particleIdHyp());

      // get momentum and energy from range (for contained particles)
      track_range_mom_muassumption->push_back(_trkmom.GetTrackMomentum(track.Length(),13)*1000.);
      track_range_mom_passumption->push_back(_trkmom.GetTrackMomentum(track.Length(),2212)*1000.);
      track_range_energy_muassumption->push_back((std::sqrt(std::pow(_trkmom.GetTrackMomentum(track.Length(),13)*1000.,2)+std::pow(105.7,2))-105.7)* 1000.);
      track_range_energy_passumption ->push_back((std::sqrt(std::pow(_trkmom.GetTrackMomentum(track.Length(),2212)*1000.,2)+std::pow(938.272,2))-938.272)*1000.);

      if (pids.size() == 0){

        std::cout << "[NuMuSelection] No track<->PID association found for track with ID " << track.ID() << ". Throwing exception. " << std::endl;
        throw std::exception();

      }

      // get PID information
      std::vector< anab::sParticleIDAlgScores > algScoresVec = pids.at(0)->ParticleIDAlgScores();

      for (size_t i_alg = 0; i_alg < algScoresVec.size(); i_alg++){

        anab::sParticleIDAlgScores algScore = algScoresVec.at(i_alg);

        if(algScore.fAlgName == "BraggPeakLLH" && UBPID::uB_getSinglePlane(algScore.fPlaneID) == 2){                                   
          if (anab::kVariableType(algScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(algScore.fTrackDir) == anab::kForward){     

            if (algScore.fAssumedPdg == 13  ) bragg_fwd_mu->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 2212) bragg_fwd_p->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 211 ) bragg_fwd_pi->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 321 ) bragg_fwd_k->push_back(algScore.fValue);   
            if (algScore.fAssumedPdg == 0   ) noBragg_fwd_mip->push_back(algScore.fValue);   

          }                                                                        
          else if (anab::kVariableType(algScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(algScore.fTrackDir) == anab::kBackward){

            if (algScore.fAssumedPdg == 13  ) bragg_bwd_mu->push_back(algScore.fValue);      
            if (algScore.fAssumedPdg == 2212) bragg_bwd_p->push_back(algScore.fValue);      
            if (algScore.fAssumedPdg == 211 ) bragg_bwd_pi->push_back(algScore.fValue);      
            if (algScore.fAssumedPdg == 321 ) bragg_bwd_k->push_back(algScore.fValue);      

          }                                                                        

        }

        if (algScore.fAlgName == "DepEvsRangeE"
            && anab::kVariableType(algScore.fVariableType) == anab::kEdeposited
            && UBPID::uB_getSinglePlane(algScore.fPlaneID) == 2){

            track_dep_energy->push_back(algScore.fValue / 1000.);


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
      track_chi2->push_back(track.Chi2());
      track_ndof->push_back(track.Ndof());
      track_ntrajpoints->push_back(track.NumberTrajectoryPoints());

      for (size_t j = 0; j < smearedCalos.size(); j++){

        if (smearedCalos.at(j)->PlaneID().Plane == 2){

          track_dedxperhit_unsmeared->push_back(unsmearedCalos.at(j)->dEdx());
          track_dedxperhit_smeared->push_back(smearedCalos.at(j)->dEdx());
          track_resrangeperhit->push_back(smearedCalos.at(j)->ResidualRange());

        }
        else continue;

      }

      if (isSimulation){
        art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit(hitHandle, e, fHitTruthAssn);
        std::pair< const simb::MCParticle*, double > mcpInformation = GetAssociatedMCParticle(hits, particlesPerHit);

        const simb::MCParticle* mcpMatchedParticle = mcpInformation.first;

        true_match_purity->push_back(mcpInformation.second);
        true_match_pdg->push_back(mcpMatchedParticle->PdgCode());
        true_match_statuscode->push_back(mcpMatchedParticle->StatusCode());
        true_match_trackid->push_back(mcpMatchedParticle->TrackId());
        true_match_motherid->push_back(mcpMatchedParticle->Mother());

        // need something smart to deal with these two
        true_match_process->push_back(mcpMatchedParticle->Process());
        true_match_endprocess->push_back(mcpMatchedParticle->EndProcess());

        true_match_numdaugters->push_back(mcpMatchedParticle->NumberDaughters());

        std::vector<int> daughterids = {0};
        for (int j = 0; j < mcpMatchedParticle->NumberDaughters(); j++){

          daughterids.push_back(mcpMatchedParticle->Daughter(j));  

        }

        true_match_daughterids->push_back(daughterids);

        true_match_startx->push_back(mcpMatchedParticle->Vx());
        true_match_starty->push_back(mcpMatchedParticle->Vy());
        true_match_startz->push_back(mcpMatchedParticle->Vz());
        true_match_startt->push_back(mcpMatchedParticle->T());
        true_match_endx->push_back(mcpMatchedParticle->EndX());
        true_match_endy->push_back(mcpMatchedParticle->EndY());
        true_match_endz->push_back(mcpMatchedParticle->EndZ());
        true_match_endt->push_back(mcpMatchedParticle->EndT());
        true_match_mass->push_back(mcpMatchedParticle->Mass());
        true_match_starte->push_back(mcpMatchedParticle->E());
        true_match_startp->push_back(mcpMatchedParticle->P());
        true_match_startpx->push_back(mcpMatchedParticle->Px());
        true_match_startpy->push_back(mcpMatchedParticle->Py());
        true_match_startpz->push_back(mcpMatchedParticle->Pz());
        true_match_ende->push_back(mcpMatchedParticle->EndE());
        true_match_endpx->push_back(mcpMatchedParticle->EndPx());
        true_match_endpy->push_back(mcpMatchedParticle->EndPy());
        true_match_endpz->push_back(mcpMatchedParticle->EndPz());
        true_match_rescatter->push_back(mcpMatchedParticle->Rescatter());

      }
    }

  }
  else if (nSelectedTPCObjects == 0){
    // this is probably fine

    std::cout << "[NuMuSelection] Event is not UBXSec selected" << std::endl;
    isUBXSecSelected = false;
    emplaceDummyVars();

  }
  // this could indicate a problem
  else if (nSelectedTPCObjects > 1){

    std::cout << "[NuMuSelection] Event has " << nSelectedTPCObjects << " TPC Objects! There should only be one! Something is very wrong! Aaaaahh!" <<  std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "[NuMuSelection] Event has " << nSelectedTPCObjects << " TPC Objects when there should only be one. Something is very wrong!" << std::endl;
    throw std::exception();

  }

  std::cout << "[NuMuSelection] Filling Tree... " << std::endl;
  ana_tree->Fill();

}

void NuMuSelection1muNpAnalyser::beginJob()
{

  ana_tree = tfs->make<TTree>("analysis_tree", "analysis tree");
  initialiseAnalysisTree(ana_tree, fIsData);

  sr_tree = tfs->make<TTree>("sr_tree", "subrun POT-utility tree");
  sr_tree->Branch("sr_run", &sr_run);
  sr_tree->Branch("sr_subrun", &sr_subrun);
  sr_tree->Branch("sr_begintime", &sr_begintime);
  sr_tree->Branch("sr_endtime", &sr_endtime);
  sr_tree->Branch("sr_pot", &sr_pot);

}

void NuMuSelection1muNpAnalyser::endSubRun(art::SubRun const & sr)
{

  sr_run = sr.run();
  sr_subrun = sr.subRun();
  sr_begintime = sr.beginTime().value();
  sr_endtime = sr.endTime().value();

  art::Handle<sumdata::POTSummary> potsum_h;
  if (!fIsData){
    if (sr.getByLabel("generator", potsum_h)){
      sr_pot = potsum_h->totpot;
    }
    else sr_pot = 0;
  }

  if (fIsData){
    if (sr.getByLabel("beamdata", potsum_h)){
      sr_pot = potsum_h->totpot;
    }
    else sr_pot = 0;
  }

  sr_tree->Fill();

}

void NuMuSelection1muNpAnalyser::initialiseAnalysisTree(TTree *tree, bool fIsData){

  tree->Branch("run"                          , &run                                              );
  tree->Branch("subrun"                       , &subrun                                           );
  tree->Branch("event"                        , &event                                            );
  tree->Branch("isData"                       , &isData                                           );
  tree->Branch("isSimulation"                 , &isSimulation                                     );
  tree->Branch("isUBXSecSelected"             , &isUBXSecSelected                                 );
  tree->Branch("nSelectedTracks"              , &nSelectedTracks                                  );
  tree->Branch("nSelectedShowers"             , &nSelectedShowers                                 );
  tree->Branch("nSelectedPfparticles"         , &nSelectedPfparticles);
  tree->Branch("noBragg_fwd_mip"              , "std::vector<double>"                                , &noBragg_fwd_mip);
  tree->Branch("bragg_fwd_mu"                 , "std::vector<double>"                                , &bragg_fwd_mu   );
  tree->Branch("bragg_fwd_p"                  , "std::vector<double>"                                , &bragg_fwd_p    );
  tree->Branch("bragg_fwd_pi"                 , "std::vector<double>"                                , &bragg_fwd_pi   );
  tree->Branch("bragg_fwd_k"                  , "std::vector<double>"                                , &bragg_fwd_k    );
  tree->Branch("bragg_bwd_mu"                 , "std::vector<double>"                                , &bragg_bwd_mu   );
  tree->Branch("bragg_bwd_p"                  , "std::vector<double>"                                , &bragg_bwd_p    );
  tree->Branch("bragg_bwd_pi"                 , "std::vector<double>"                                , &bragg_bwd_pi   );
  tree->Branch("bragg_bwd_k"                  , "std::vector<double>"                                , &bragg_bwd_k    );
  tree->Branch("track_isContained"            , "std::vector<bool>"                                  , &track_isContained);
  tree->Branch("track_dep_energy"             , "std::vector<double>"                                , &track_dep_energy  );
  tree->Branch("track_length"                 , "std::vector<double>"                                , &track_length   );
  tree->Branch("track_theta"                  , "std::vector<double>"                                , &track_theta    );
  tree->Branch("track_costheta"               , "std::vector<double>"                                , &track_costheta );
  tree->Branch("track_phi"                    , "std::vector<double>"                                , &track_phi      );
  tree->Branch("track_startx"                 , "std::vector<double>"                                , &track_startx   );
  tree->Branch("track_starty"                 , "std::vector<double>"                                , &track_starty   );
  tree->Branch("track_startz"                 , "std::vector<double>"                                , &track_startz   );
  tree->Branch("track_endx"                   , "std::vector<double>"                                , &track_endx     );
  tree->Branch("track_endy"                   , "std::vector<double>"                                , &track_endy     );
  tree->Branch("track_endz"                   , "std::vector<double>"                                , &track_endz     );
  tree->Branch("track_chi2"                   , "std::vector<double>"                                , &track_chi2        );
  tree->Branch("track_ndof"                   , "std::vector<int>"                                   , &track_ndof        );
  tree->Branch("track_ntrajpoints"            , "std::vector<int>"                                   , &track_ntrajpoints );
  tree->Branch("track_dedxperhit_smeared"     , "std::vector< std::vector<double> >"                 , &track_dedxperhit_smeared);
  tree->Branch("track_dedxperhit_unsmeared"   , "std::vector< std::vector<double> >"                 , &track_dedxperhit_unsmeared);
  tree->Branch("track_resrangeperhit"         , "std::vector< std::vector<double> >"                 , & track_resrangeperhit);
  tree->Branch("vertex_x"                     , &vertex_x);
  tree->Branch("vertex_y"                     , &vertex_y);
  tree->Branch("vertex_z"                     , &vertex_z);
  tree->Branch("track_mcs_muassmp_fwd"                , "std::vector<double>"                                , &track_mcs_muassmp_fwd);
  tree->Branch("track_mcs_muassmp_bwd"                , "std::vector<double>"                                , &track_mcs_muassmp_bwd);
  tree->Branch("track_mcs_muassmp_energy_fwd"         , "std::vector<double>"                                , &track_mcs_muassmp_energy_fwd);
  tree->Branch("track_mcs_muassmp_energy_bwd"         , "std::vector<double>"                                , &track_mcs_muassmp_energy_bwd);
  tree->Branch("track_mcs_muassmp_fwd_uncertainty"    , "std::vector<double>"                                , &track_mcs_muassmp_fwd_uncertainty);
  tree->Branch("track_mcs_muassmp_bwd_uncertainty"    , "std::vector<double>"                                , &track_mcs_muassmp_bwd_uncertainty);
  tree->Branch("track_mcs_muassmp_particlehypothesis" , "std::vector<int>"                                   , &track_mcs_muassmp_particlehypothesis);
  tree->Branch("track_mcs_muassmp_fwd_loglikelihood"  , "std::vector<double>"                                , &track_mcs_muassmp_fwd_loglikelihood);
  tree->Branch("track_mcs_muassmp_bwd_loglikelihood"  , "std::vector<double>"                                , &track_mcs_muassmp_bwd_loglikelihood);
  tree->Branch("track_range_mom_muassumption"         , "std::vector<double>"                                , &track_range_mom_muassumption);
  tree->Branch("track_range_mom_passumption"          , "std::vector<double>"                                , &track_range_mom_passumption);
  tree->Branch("track_range_energy_muassumption"      , "std::vector<double>"                                , &track_range_energy_muassumption);
  tree->Branch("track_range_energy_passumption"       , "std::vector<double>"                                , &track_range_energy_passumption);

  if (!fIsData){
    tree->Branch("isBeamNeutrino"            , &isBeamNeutrino);
    tree->Branch("isCosmic"                  , &isCosmic);
    tree->Branch("isMixed"                   , &isMixed);
    tree->Branch("isInFV"                    , &isInFV);
    tree->Branch("true_truth_nparticles"     , &true_truth_nparticles);
    tree->Branch("true_genie_pdg"            , "std::vector<int>"              , &true_genie_pdg         );
    tree->Branch("true_genie_statuscode "    , "std::vector<int>"              , &true_genie_statuscode  );
    tree->Branch("true_genie_trackid"        , "std::vector<int>"              , &true_genie_trackid     );
    tree->Branch("true_genie_motherid"       , "std::vector<int>"              , &true_genie_motherid    );
    tree->Branch("true_genie_process"        , "std::vector<std::string>"      , &true_genie_process     );
    tree->Branch("true_genie_endprocess"     , "std::vector<std::string>"      , &true_genie_endprocess  );
    tree->Branch("true_genie_numdaugters"    , "std::vector<int>"              , &true_genie_numdaugters );
    tree->Branch("true_genie_daughterids"    , "std::vector<std::vector<int>>" , &true_genie_daughterids );
    tree->Branch("true_genie_startx"         , "std::vector<double>"           , &true_genie_startx      );
    tree->Branch("true_genie_starty"         , "std::vector<double>"           , &true_genie_starty      );
    tree->Branch("true_genie_startz"         , "std::vector<double>"           , &true_genie_startz      );
    tree->Branch("true_genie_startt"         , "std::vector<double>"           , &true_genie_startt      );
    tree->Branch("true_genie_endx"           , "std::vector<double>"           , &true_genie_endx        );
    tree->Branch("true_genie_endy"           , "std::vector<double>"           , &true_genie_endy        );
    tree->Branch("true_genie_endz"           , "std::vector<double>"           , &true_genie_endz        );
    tree->Branch("true_genie_endt"           , "std::vector<double>"           , &true_genie_endt        );
    tree->Branch("true_genie_mass"           , "std::vector<double>"           , &true_genie_mass        );
    tree->Branch("true_genie_starte"         , "std::vector<double>"           , &true_genie_starte      );
    tree->Branch("true_genie_startp"         , "std::vector<double>"           , &true_genie_startp      );
    tree->Branch("true_genie_startpx"        , "std::vector<double>"           , &true_genie_startpx     );
    tree->Branch("true_genie_startpy"        , "std::vector<double>"           , &true_genie_startpy     );
    tree->Branch("true_genie_startpz"        , "std::vector<double>"           , &true_genie_startpz     );
    tree->Branch("true_genie_ende"           , "std::vector<double>"           , &true_genie_ende        );
    tree->Branch("true_genie_endpx"          , "std::vector<double>"           , &true_genie_endpx       );
    tree->Branch("true_genie_endpy"          , "std::vector<double>"           , &true_genie_endpy       );
    tree->Branch("true_genie_endpz"          , "std::vector<double>"           , &true_genie_endpz       );
    tree->Branch("true_genie_rescatter"      , "std::vector<int>"              , &true_genie_rescatter   );
    tree->Branch("true_mcp_pdg"              , "std::vector<int>"              , &true_mcp_pdg         );
    tree->Branch("true_mcp_statuscode "      , "std::vector<int>"              , &true_mcp_statuscode  );
    tree->Branch("true_mcp_trackid"          , "std::vector<int>"              , &true_mcp_trackid     );
    tree->Branch("true_mcp_motherid"         , "std::vector<int>"              , &true_mcp_motherid    );
    tree->Branch("true_mcp_process"          , "std::vector<std::string>"      , &true_mcp_process     );
    tree->Branch("true_mcp_endprocess"       , "std::vector<std::string>"      , &true_mcp_endprocess  );
    tree->Branch("true_mcp_numdaugters"      , "std::vector<int>"              , &true_mcp_numdaugters );
    tree->Branch("true_mcp_daughterids"      , "std::vector<std::vector<int>>" , &true_mcp_daughterids );
    tree->Branch("true_mcp_startx"           , "std::vector<double>"           , &true_mcp_startx      );
    tree->Branch("true_mcp_starty"           , "std::vector<double>"           , &true_mcp_starty      );
    tree->Branch("true_mcp_startz"           , "std::vector<double>"           , &true_mcp_startz      );
    tree->Branch("true_mcp_startt"           , "std::vector<double>"           , &true_mcp_startt      );
    tree->Branch("true_mcp_endx"             , "std::vector<double>"           , &true_mcp_endx        );
    tree->Branch("true_mcp_endy"             , "std::vector<double>"           , &true_mcp_endy        );
    tree->Branch("true_mcp_endz"             , "std::vector<double>"           , &true_mcp_endz        );
    tree->Branch("true_mcp_endt"             , "std::vector<double>"           , &true_mcp_endt        );
    tree->Branch("true_mcp_mass"             , "std::vector<double>"           , &true_mcp_mass        );
    tree->Branch("true_mcp_starte"           , "std::vector<double>"           , &true_mcp_starte      );
    tree->Branch("true_mcp_startp"           , "std::vector<double>"           , &true_mcp_startp      );
    tree->Branch("true_mcp_startpx"          , "std::vector<double>"           , &true_mcp_startpx     );
    tree->Branch("true_mcp_startpy"          , "std::vector<double>"           , &true_mcp_startpy     );
    tree->Branch("true_mcp_startpz"          , "std::vector<double>"           , &true_mcp_startpz     );
    tree->Branch("true_mcp_ende"             , "std::vector<double>"           , &true_mcp_ende        );
    tree->Branch("true_mcp_endpx"            , "std::vector<double>"           , &true_mcp_endpx       );
    tree->Branch("true_mcp_endpy"            , "std::vector<double>"           , &true_mcp_endpy       );
    tree->Branch("true_mcp_endpz"            , "std::vector<double>"           , &true_mcp_endpz       );
    tree->Branch("true_mcp_rescatter"        , "std::vector<int>"              , &true_mcp_rescatter   );
    tree->Branch("true_nu_ccnc"              , &true_nu_ccnc);
    tree->Branch("true_nu_mode"              , &true_nu_mode);
    tree->Branch("true_nu_interactiontype"   , &true_nu_interactiontype);
    tree->Branch("true_nu_target"            , &true_nu_target);
    tree->Branch("true_nu_hitnuc"            , &true_nu_hitnuc);
    tree->Branch("true_nu_hitquark"          , &true_nu_hitquark);
    tree->Branch("true_nu_w"                 , &true_nu_w);
    tree->Branch("true_nu_x"                 , &true_nu_x);
    tree->Branch("true_nu_y"                 , &true_nu_y);
    tree->Branch("true_nu_qsqr"              , &true_nu_qsqr);
    tree->Branch("true_nu_pt"                , &true_nu_pt);
    tree->Branch("true_nu_theta"             , &true_nu_theta);
    tree->Branch("true_match_purity"         , "std::vector<double>"           , &true_match_purity      );
    tree->Branch("true_match_pdg"            , "std::vector<int>"              , &true_match_pdg         );
    tree->Branch("true_match_statuscode "    , "std::vector<int>"              , &true_match_statuscode  );
    tree->Branch("true_match_trackid"        , "std::vector<int>"              , &true_match_trackid     );
    tree->Branch("true_match_motherid"       , "std::vector<int>"              , &true_match_motherid    );
    tree->Branch("true_match_process"        , "std::vector<std::string>"      , &true_match_process     );
    tree->Branch("true_match_endprocess"     , "std::vector<std::string>"      , &true_match_endprocess  );
    tree->Branch("true_match_numdaugters"    , "std::vector<int>"              , &true_match_numdaugters );
    tree->Branch("true_match_daughterids"    , "std::vector<std::vector<int>>" , &true_match_daughterids );
    tree->Branch("true_match_startx"         , "std::vector<double>"           , &true_match_startx      );
    tree->Branch("true_match_starty"         , "std::vector<double>"           , &true_match_starty      );
    tree->Branch("true_match_startz"         , "std::vector<double>"           , &true_match_startz      );
    tree->Branch("true_match_startt"         , "std::vector<double>"           , &true_match_startt      );
    tree->Branch("true_match_endx"           , "std::vector<double>"           , &true_match_endx        );
    tree->Branch("true_match_endy"           , "std::vector<double>"           , &true_match_endy        );
    tree->Branch("true_match_endz"           , "std::vector<double>"           , &true_match_endz        );
    tree->Branch("true_match_endt"           , "std::vector<double>"           , &true_match_endt        );
    tree->Branch("true_match_mass"           , "std::vector<double>"           , &true_match_mass        );
    tree->Branch("true_match_starte"         , "std::vector<double>"           , &true_match_starte      );
    tree->Branch("true_match_startp"         , "std::vector<double>"           , &true_match_startp      );
    tree->Branch("true_match_startpx"        , "std::vector<double>"           , &true_match_startpx     );
    tree->Branch("true_match_startpy"        , "std::vector<double>"           , &true_match_startpy     );
    tree->Branch("true_match_startpz"        , "std::vector<double>"           , &true_match_startpz     );
    tree->Branch("true_match_ende"           , "std::vector<double>"           , &true_match_ende        );
    tree->Branch("true_match_endpx"          , "std::vector<double>"           , &true_match_endpx       );
    tree->Branch("true_match_endpy"          , "std::vector<double>"           , &true_match_endpy       );
    tree->Branch("true_match_endpz"          , "std::vector<double>"           , &true_match_endpz       );
    tree->Branch("true_match_rescatter"      , "std::vector<int>"              , &true_match_rescatter   );
  }
}

void NuMuSelection1muNpAnalyser::resizeVectors(bool isData){

  noBragg_fwd_mip->resize(0);
  bragg_fwd_mu->resize(0);
  bragg_fwd_p->resize(0);
  bragg_fwd_pi->resize(0);
  bragg_fwd_k->resize(0);
  bragg_bwd_mu->resize(0);
  bragg_bwd_p->resize(0);
  bragg_bwd_pi->resize(0);
  bragg_bwd_k->resize(0);
  track_isContained->resize(0);
  track_dep_energy->resize(0);
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
  track_chi2->resize(0);
  track_ndof->resize(0);
  track_ntrajpoints->resize(0);
  track_dedxperhit_smeared->resize(0);
  track_dedxperhit_unsmeared->resize(0);
  track_resrangeperhit->resize(0);
  track_mcs_muassmp_fwd->resize(0);
  track_mcs_muassmp_bwd->resize(0);
  track_mcs_muassmp_energy_fwd->resize(0);
  track_mcs_muassmp_energy_bwd->resize(0);
  track_mcs_muassmp_fwd_uncertainty->resize(0);
  track_mcs_muassmp_bwd_uncertainty->resize(0);
  track_mcs_muassmp_particlehypothesis->resize(0);
  track_mcs_muassmp_fwd_loglikelihood->resize(0);
  track_mcs_muassmp_bwd_loglikelihood->resize(0);
  track_range_mom_muassumption->resize(0);
  track_range_mom_passumption->resize(0);
  track_range_energy_muassumption->resize(0);
  track_range_energy_passumption->resize(0);
  if (!isData){
    true_genie_pdg->resize(0);
    true_genie_statuscode->resize(0);
    true_genie_trackid->resize(0);
    true_genie_motherid->resize(0);
    true_genie_process->resize(0);
    true_genie_endprocess->resize(0);
    true_genie_numdaugters->resize(0);
    true_genie_daughterids->resize(0);
    true_genie_startx->resize(0);
    true_genie_starty->resize(0);
    true_genie_startz->resize(0);
    true_genie_startt->resize(0);
    true_genie_endx->resize(0);
    true_genie_endy->resize(0);
    true_genie_endz->resize(0);
    true_genie_endt->resize(0);
    true_genie_mass->resize(0);
    true_genie_starte->resize(0);
    true_genie_startp->resize(0);
    true_genie_startpx->resize(0);
    true_genie_startpy->resize(0);
    true_genie_startpz->resize(0);
    true_genie_ende->resize(0);
    true_genie_endpx->resize(0);
    true_genie_endpy->resize(0);
    true_genie_endpz->resize(0);
    true_genie_rescatter->resize(0);
    true_mcp_pdg->resize(0);
    true_mcp_statuscode->resize(0);
    true_mcp_trackid->resize(0);
    true_mcp_motherid->resize(0);
    true_mcp_process->resize(0);
    true_mcp_endprocess->resize(0);
    true_mcp_numdaugters->resize(0);
    true_mcp_daughterids->resize(0);
    true_mcp_startx->resize(0);
    true_mcp_starty->resize(0);
    true_mcp_startz->resize(0);
    true_mcp_startt->resize(0);
    true_mcp_endx->resize(0);
    true_mcp_endy->resize(0);
    true_mcp_endz->resize(0);
    true_mcp_endt->resize(0);
    true_mcp_mass->resize(0);
    true_mcp_starte->resize(0);
    true_mcp_startp->resize(0);
    true_mcp_startpx->resize(0);
    true_mcp_startpy->resize(0);
    true_mcp_startpz->resize(0);
    true_mcp_ende->resize(0);
    true_mcp_endpx->resize(0);
    true_mcp_endpy->resize(0);
    true_mcp_endpz->resize(0);
    true_mcp_rescatter->resize(0);
    true_match_purity->resize(0);
    true_match_pdg->resize(0);
    true_match_statuscode->resize(0);
    true_match_trackid->resize(0);
    true_match_motherid->resize(0);
    true_match_process->resize(0);
    true_match_endprocess->resize(0);
    true_match_numdaugters->resize(0);
    true_match_daughterids->resize(0);
    true_match_startx->resize(0);
    true_match_starty->resize(0);
    true_match_startz->resize(0);
    true_match_startt->resize(0);
    true_match_endx->resize(0);
    true_match_endy->resize(0);
    true_match_endz->resize(0);
    true_match_endt->resize(0);
    true_match_mass->resize(0);
    true_match_starte->resize(0);
    true_match_startp->resize(0);
    true_match_startpx->resize(0);
    true_match_startpy->resize(0);
    true_match_startpz->resize(0);
    true_match_ende->resize(0);
    true_match_endpx->resize(0);
    true_match_endpy->resize(0);
    true_match_endpz->resize(0);
    true_match_rescatter->resize(0);
  }

}

void NuMuSelection1muNpAnalyser::emplaceDummyVars(){

  noBragg_fwd_mip->resize(1, -999);
  bragg_fwd_mu->resize(1, -999);
  bragg_fwd_p->resize(1, -999);
  bragg_fwd_pi->resize(1, -999);
  bragg_fwd_k->resize(1, -999);
  bragg_bwd_mu->resize(1, -999);
  bragg_bwd_p->resize(1, -999);
  bragg_bwd_pi->resize(1, -999);
  bragg_bwd_k->resize(1, -999);
  track_isContained->resize(1,-999);
  track_dep_energy->resize(1,-999);
  track_length->resize(1, -999);
  track_theta->resize(1, -999);
  track_costheta->resize(1, -999);
  track_phi->resize(1, -999);
  track_startx->resize(1, -999);
  track_starty->resize(1, -999);
  track_startz->resize(1, -999);
  track_endx->resize(1, -999);
  track_endy->resize(1, -999);
  track_endz->resize(1, -999);
  track_chi2->resize(1, -999);
  track_ndof->resize(1, -999);
  track_ntrajpoints->resize(1, -999);
  track_dedxperhit_smeared->resize(1, {{-999.}});
  track_dedxperhit_unsmeared->resize(1, {{-999.}});
  track_resrangeperhit->resize(1, {{-999.}});
  track_mcs_muassmp_fwd->resize(1, -999);
  track_mcs_muassmp_bwd->resize(1, -999);
  track_mcs_muassmp_energy_fwd->resize(1, -999);
  track_mcs_muassmp_energy_bwd->resize(1, -999);
  track_mcs_muassmp_fwd_uncertainty->resize(1, -999);
  track_mcs_muassmp_bwd_uncertainty->resize(1, -999);
  track_mcs_muassmp_particlehypothesis->resize(1, -999);
  track_mcs_muassmp_fwd_loglikelihood->resize(1, -999);
  track_mcs_muassmp_bwd_loglikelihood->resize(1, -999);
  track_range_mom_muassumption->resize(1,-999);
  track_range_mom_passumption->resize(1,-999);
  track_range_energy_muassumption->resize(1,-999);
  track_range_energy_passumption->resize(1,-999);
  true_match_purity->resize(1, -999);
  true_match_pdg->resize(1, -999);
  true_match_statuscode->resize(1, -999);
  true_match_trackid->resize(1, -999);
  true_match_motherid->resize(1, -999);
  true_match_process->resize(1, "dummy");
  true_match_endprocess->resize(1, "dummy");
  true_match_numdaugters->resize(1, -999);
  true_match_daughterids->resize(1, {{-999}});
  true_match_startx->resize(1, -999);
  true_match_starty->resize(1, -999);
  true_match_startz->resize(1, -999);
  true_match_startt->resize(1, -999);
  true_match_endx->resize(1, -999);
  true_match_endy->resize(1, -999);
  true_match_endz->resize(1, -999);
  true_match_endt->resize(1, -999);
  true_match_mass->resize(1, -999);
  true_match_starte->resize(1, -999);
  true_match_startp->resize(1, -999);
  true_match_startpx->resize(1, -999);
  true_match_startpy->resize(1, -999);
  true_match_startpz->resize(1, -999);
  true_match_ende->resize(1, -999);
  true_match_endpx->resize(1, -999);
  true_match_endpy->resize(1, -999);
  true_match_endpz->resize(1, -999);
  true_match_rescatter->resize(1, -999);
  nSelectedTracks = -999;
  nSelectedShowers = -999;
  nSelectedPfparticles = -999;
  isBeamNeutrino = false;
  isCosmic = false;
  isMixed = false;
  vertex_x = -999;
  vertex_y = -999;
  vertex_z = -999;


}

std::pair< const simb::MCParticle*, double > NuMuSelection1muNpAnalyser::GetAssociatedMCParticle(std::vector< art::Ptr< recob::Hit > > hits, art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit)
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

DEFINE_ART_MODULE(NuMuSelection1muNpAnalyser)
