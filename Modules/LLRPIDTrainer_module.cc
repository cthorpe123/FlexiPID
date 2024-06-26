////////////////////////////////////////////////////////////////////////
// Class:       LLRPIDTrainer
// Plugin Type: analyzer (art v3_03_01)
// File:        LLRPIDTrainer_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

// C++ STL includes
#include <vector>
#include <string>

// larsoft/uboonecode includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// root includes
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

namespace FlexiPID {
  class LLRPIDTrainer;
}

class FlexiPID::LLRPIDTrainer : public art::EDAnalyzer {
  public:
    explicit LLRPIDTrainer(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    LLRPIDTrainer(LLRPIDTrainer const&) = delete;
    LLRPIDTrainer(LLRPIDTrainer&&) = delete;
    LLRPIDTrainer& operator=(LLRPIDTrainer const&) = delete;
    LLRPIDTrainer& operator=(LLRPIDTrainer&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

    void beginSubRun(const art::SubRun& sr);
    void endSubRun(const art::SubRun& sr);

  private:

    // Output trees
    TTree * OutputTree;
    int t_run,t_subrun,t_event;
    std::vector<int> t_TrackTruePDG;
    std::vector<double> t_TrackTrueMomentum;
    std::vector<double> t_TrackTruthPurity;
    std::vector<float> t_TrackLength;
    std::vector<float> t_TrackAngleTrueReco; // Angle between true and reco track directions    
    std::vector<std::vector<float>> t_ResidualRange_Plane0;
    std::vector<std::vector<float>> t_dEdx_Plane0;
    std::vector<std::vector<float>> t_Pitch_Plane0;
    std::vector<std::vector<float>> t_ResidualRange_Plane1;
    std::vector<std::vector<float>> t_dEdx_Plane1;
    std::vector<std::vector<float>> t_Pitch_Plane1;
    std::vector<std::vector<float>> t_ResidualRange_Plane2;
    std::vector<std::vector<float>> t_dEdx_Plane2;
    std::vector<std::vector<float>> t_Pitch_Plane2;

    // Fhicl parameters
    const std::string f_PFParticleModuleLabel;
    const std::string f_TrackModuleLabel;
    const std::string f_HitModuleLabel;
    const std::string f_CaloModuleLabel;
    const std::string f_TrackHitAssnLabel;
    const std::string f_HitTruthAssnLabel; 

};

FlexiPID::LLRPIDTrainer::LLRPIDTrainer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  f_PFParticleModuleLabel(p.get<std::string>("PFParticleModuleLabel")),
  f_TrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
  f_HitModuleLabel(p.get<std::string>("HitModuleLabel")),
  f_CaloModuleLabel(p.get<std::string>("CaloModuleLabel")),
  f_TrackHitAssnLabel(p.get<std::string>("TrackHitAssnLabel")),
  f_HitTruthAssnLabel(p.get<std::string>("HitTruthAssnLabel"))
{

}

void FlexiPID::LLRPIDTrainer::analyze(art::Event const& e)
{

  // Begin by resetting everything

  t_run = e.run();
  t_subrun = e.subRun();
  t_event = e.event();

  t_TrackTruePDG.clear();
  t_TrackTrueMomentum.clear();
  t_TrackTruthPurity.clear();
  t_TrackLength.clear();
  t_TrackAngleTrueReco.clear();
  t_ResidualRange_Plane0.clear();
  t_dEdx_Plane0.clear();
  t_Pitch_Plane0.clear();
  t_ResidualRange_Plane1.clear();
  t_dEdx_Plane1.clear();
  t_Pitch_Plane1.clear();
  t_ResidualRange_Plane2.clear();
  t_dEdx_Plane2.clear();
  t_Pitch_Plane2.clear();

  // Load all the various data products from the event

  art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
  std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;
  if(!e.getByLabel(f_PFParticleModuleLabel,Handle_PFParticle)) 
    throw cet::exception("LLRPIDTrainer") << "No PFParticle Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);

  art::Handle<std::vector<recob::Track>> Handle_Tracks;
  std::vector<art::Ptr<recob::Track>> Vect_Tracks;
  if(!e.getByLabel(f_TrackModuleLabel,Handle_Tracks))  
    throw cet::exception("LLRPIDTrainer") << "No Track data product!" << std::endl;
  art::fill_ptr_vector(Vect_Tracks,Handle_Tracks);

  art::Handle<std::vector<anab::Calorimetry>> Handle_Calorimetry;
  std::vector<art::Ptr<anab::Calorimetry>> Vect_Calorimetry;
  if(!e.getByLabel(f_CaloModuleLabel,Handle_Calorimetry))  
    throw cet::exception("LLRPIDTrainer") << "No Calorimetry data product!" << std::endl;
  art::fill_ptr_vector(Vect_Calorimetry,Handle_Calorimetry);

  art::Handle<std::vector<recob::Hit>> Handle_Hit;
  std::vector<art::Ptr<recob::Hit>> Vect_Hit;
  if(!e.getByLabel(f_HitModuleLabel,Handle_Hit)) 
    throw cet::exception("LLRPIDTrainer") << "No Hit Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_Hit,Handle_Hit);


  art::FindManyP<recob::Track> Assoc_PFParticleTrack(Vect_PFParticle,e,f_TrackModuleLabel); 
  art::FindManyP<recob::Hit> Assoc_TrackHit(Vect_Tracks,e,f_TrackHitAssnLabel);
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> ParticlesPerHit(Handle_Hit,e,f_HitTruthAssnLabel);
  art::FindManyP<anab::Calorimetry> Assoc_TrackCalo(Vect_Tracks,e,f_CaloModuleLabel); 

  // Check if there is reco'd neutrino 
  size_t neutrinoID = 99999;
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->IsPrimary() && (abs(pfp->PdgCode()) == 14 || abs(pfp->PdgCode()) == 12)){
      neutrinoID = pfp->Self();
    }
  }
  if(neutrinoID == 99999) return;

  // Loop over PFParticles
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->Parent() != neutrinoID) continue; 

    // Grab the associated track
    std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack.at(pfp.key());
    if(pfpTracks.size() != 1) continue;
    art::Ptr<recob::Track> track = pfpTracks.at(0);

    // Truth match the track
    std::vector<art::Ptr<recob::Hit>> hits = Assoc_TrackHit.at(track.key());
    std::unordered_map<int,double>  trkide;
    int maxhits=-1;
    simb::MCParticle const* matchedParticle = NULL;
    std::vector<simb::MCParticle const*> particleVec;
    std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

    for(size_t i_hit=0;i_hit<hits.size();++i_hit){
      particleVec.clear();
      matchVec.clear();
      ParticlesPerHit.get(hits[i_hit].key(),particleVec,matchVec);
      for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){
        trkide[particleVec[i_particle]->TrackId()]++; 
        if(trkide[particleVec[i_particle]->TrackId()] > maxhits){
          maxhits = trkide[particleVec[i_particle]->TrackId()];
          matchedParticle = particleVec[i_particle];
        }
      }
    }
    if(matchedParticle == NULL) continue;

    t_TrackTruePDG.push_back(matchedParticle->PdgCode());
    t_TrackTrueMomentum.push_back(TVector3(matchedParticle->Px(),matchedParticle->Py(),matchedParticle->Pz()).Mag());
    t_TrackTruthPurity.push_back((double)maxhits/hits.size());
    t_TrackLength.push_back(track->Length());

    // Calculate the angle between the true and reconstructed directions - check if particle is reco'd in the
    // correct direction later
    TVector3 true_dir(matchedParticle->Px(),matchedParticle->Py(),matchedParticle->Pz());
    TVector3 reco_dir(track->StartDirection().X(),track->StartDirection().Y(),track->StartDirection().Z());
    t_TrackAngleTrueReco.push_back(180/3.141*true_dir.Angle(reco_dir)); 

    t_ResidualRange_Plane0.push_back(std::vector<float>()); 
    t_dEdx_Plane0.push_back(std::vector<float>()); 
    t_Pitch_Plane0.push_back(std::vector<float>()); 

    t_ResidualRange_Plane1.push_back(std::vector<float>()); 
    t_dEdx_Plane1.push_back(std::vector<float>()); 
    t_Pitch_Plane1.push_back(std::vector<float>()); 

    t_ResidualRange_Plane2.push_back(std::vector<float>()); 
    t_dEdx_Plane2.push_back(std::vector<float>()); 
    t_Pitch_Plane2.push_back(std::vector<float>()); 

    // Get the calorimetry data
    std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack = Assoc_TrackCalo.at(track.key());
    for(art::Ptr<anab::Calorimetry> calo : caloFromTrack){
      int plane = calo->PlaneID().Plane;
      if(plane != 0 && plane != 1 && plane != 2) continue;        

      if(plane == 0){
        t_ResidualRange_Plane0.back() = calo->ResidualRange();
        t_dEdx_Plane0.back() = calo->dEdx();
        t_Pitch_Plane0.back() = calo->TrkPitchVec();
      }
      if(plane == 1){
        t_ResidualRange_Plane1.back() = calo->ResidualRange();
        t_dEdx_Plane1.back() = calo->dEdx();
        t_Pitch_Plane1.back() = calo->TrkPitchVec();
      }
      if(plane == 2){
        t_ResidualRange_Plane2.back() = calo->ResidualRange();
        t_dEdx_Plane2.back() = calo->dEdx();
        t_Pitch_Plane2.back() = calo->TrkPitchVec();
      }
    }

  }

  OutputTree->Fill();

}

void FlexiPID::LLRPIDTrainer::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;

  OutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");

  OutputTree->Branch("run",&t_run);
  OutputTree->Branch("subrun",&t_subrun);
  OutputTree->Branch("event",&t_event);

  OutputTree->Branch("TrackTruePDG",&t_TrackTruePDG);
  OutputTree->Branch("TrackTrueMomentum",&t_TrackTrueMomentum);
  OutputTree->Branch("TrackTruthPurity",&t_TrackTruthPurity);
  OutputTree->Branch("TrackLength",&t_TrackLength);
  OutputTree->Branch("TrackAngleTrueReco",&t_TrackAngleTrueReco);

  OutputTree->Branch("ResidualRange_Plane0",&t_ResidualRange_Plane0);
  OutputTree->Branch("dEdx_Plane0",&t_dEdx_Plane0);
  OutputTree->Branch("Pitch_Plane0",&t_Pitch_Plane0);

  OutputTree->Branch("ResidualRange_Plane1",&t_ResidualRange_Plane1);
  OutputTree->Branch("dEdx_Plane1",&t_dEdx_Plane1);
  OutputTree->Branch("Pitch_Plane1",&t_Pitch_Plane1);

  OutputTree->Branch("ResidualRange_Plane2",&t_ResidualRange_Plane2);
  OutputTree->Branch("dEdx_Plane2",&t_dEdx_Plane2);
  OutputTree->Branch("Pitch_Plane2",&t_Pitch_Plane2);

}

void FlexiPID::LLRPIDTrainer::endJob()
{
}

void FlexiPID::LLRPIDTrainer::beginSubRun(const art::SubRun& sr)
{
}

void FlexiPID::LLRPIDTrainer::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(FlexiPID::LLRPIDTrainer)
