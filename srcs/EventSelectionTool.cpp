#include "../include/EventSelectionTool.h"
#include <iostream>
#include "TLeaf.h"
#include "TBranch.h"

namespace selection{
 
  void EventSelectionTool::LoadEventList(const std::string &file_name, EventList &event_list){
  
    TFile f(file_name.c_str());
    TTree *t_event    = (TTree*) f.Get("event_tree");
    TTree *t_particle = (TTree*) f.Get("particle_tree");
    TTree *t_track    = (TTree*) f.Get("recotrack_tree");
    TTree *t_shower   = (TTree*) f.Get("recoshower_tree");

    UniqueEventIdList unique_event_list;
    EventSelectionTool::GetUniqueEventList(t_event, unique_event_list);

    unsigned int n_events = t_event->GetEntries();

    for(unsigned int i = 0; i < unique_event_list.size(); ++i){
    
      ParticleList mcparticles;
      ParticleList recoparticles;
      TrackList    tracks;
      ShowerList   showers;

      TBranch *b_event_id = t_event->GetBranch("event_id");
      TBranch *b_time_now = t_event->GetBranch("time_now");
      TBranch *b_r_vertex = t_event->GetBranch("r_vertex");
      TBranch *b_t_vertex = t_event->GetBranch("t_vertex");
      TBranch *b_t_nuance = t_event->GetBranch("t_interaction");
      TBranch *b_t_iscc   = t_event->GetBranch("t_iscc");
      
      double temp_r_vertex[3];
      double temp_t_vertex[3];
      unsigned int nuance = std::numeric_limits<unsigned int>::max();
      bool iscc(false);

      bool foundEvent(false);

      for(unsigned int j = 0; j < n_events; ++j){
     
        t_event->GetEntry(j);

        int event_id = b_event_id->GetLeaf("event_id")->GetValue();
        int time_now = b_time_now->GetLeaf("time_now")->GetValue();
        
        if(event_id != unique_event_list[i].first || time_now != unique_event_list[i].second) continue;

        foundEvent = true;

        temp_r_vertex[0]       = b_r_vertex->GetLeaf("r_vertex")->GetValue(0);
        temp_r_vertex[1]       = b_r_vertex->GetLeaf("r_vertex")->GetValue(1);
        temp_r_vertex[2]       = b_r_vertex->GetLeaf("r_vertex")->GetValue(2);
        temp_t_vertex[0]       = b_t_vertex->GetLeaf("t_vertex")->GetValue(0);
        temp_t_vertex[1]       = b_t_vertex->GetLeaf("t_vertex")->GetValue(1);
        temp_t_vertex[2]       = b_t_vertex->GetLeaf("t_vertex")->GetValue(2);
        nuance                 = b_t_nuance->GetLeaf("t_interaction")->GetValue();
        iscc                   = b_t_iscc->GetLeaf("t_iscc")->GetValue();
      
        break;

      }

      if(!foundEvent) throw 3;
      
      TVector3 r_vertex(temp_r_vertex);
      TVector3 t_vertex(temp_t_vertex);

      EventSelectionTool::GetTrackList(t_track, unique_event_list[i], tracks);
      EventSelectionTool::GetShowerList(t_shower, unique_event_list[i], showers);
      EventSelectionTool::GetMCParticleList(t_particle, unique_event_list[i], mcparticles);
      EventSelectionTool::GetRecoParticleFromTrack(tracks, recoparticles);
      EventSelectionTool::GetRecoParticleFromShower(showers, r_vertex, recoparticles);

      event_list.push_back(Event(mcparticles, recoparticles, nuance, iscc, t_vertex, r_vertex));
   
    }
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetUniqueEventList(TTree *event_tree, UniqueEventIdList &unique_event_list){
  
    TBranch *b_event_id = event_tree->GetBranch("event_id");
    TBranch *b_time_now = event_tree->GetBranch("time_now");
    
    unsigned int n_events = event_tree->GetEntries();

    for(unsigned int i = 0; i < n_events; ++i){
   
      event_tree->GetEntry(i);

      int event_id_a = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now_a = b_time_now->GetLeaf("time_now")->GetValue();
  
      bool shouldAdd(true);
    
      for(unsigned int j = 0; j < n_events; ++j){
        event_tree->GetEntry(j);

        int event_id_b = b_event_id->GetLeaf("event_id")->GetValue();
        int time_now_b = b_time_now->GetLeaf("time_now")->GetValue();
        
        if (event_id_a == event_id_b && time_now_a == time_now_b && i != j){
          shouldAdd = false;
          break;
        }
      }
       
      if (shouldAdd)
        unique_event_list.push_back(std::pair<int, int>(event_id_a, time_now_a));
    }        
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetTrackList(TTree *track_tree, const std::pair<int, int> &unique_event, TrackList &track_list){
   
    TBranch *b_event_id       = track_tree->GetBranch("event_id");
    TBranch *b_time_now       = track_tree->GetBranch("time_now");
    TBranch *b_vertex         = track_tree->GetBranch("tr_vertex");
    TBranch *b_end            = track_tree->GetBranch("tr_end");
    TBranch *b_pida           = track_tree->GetBranch("tr_pida");
    TBranch *b_chi2_mu        = track_tree->GetBranch("tr_chi2_mu");
    TBranch *b_chi2_pi        = track_tree->GetBranch("tr_chi2_pi");
    TBranch *b_chi2_pr        = track_tree->GetBranch("tr_chi2_pr");
    TBranch *b_chi2_ka        = track_tree->GetBranch("tr_chi2_ka");
    TBranch *b_length         = track_tree->GetBranch("tr_length");
    TBranch *b_kinetic_energy = track_tree->GetBranch("tr_kinetic_energy");
    
    unsigned int n_entries = track_tree->GetEntries();

    for(unsigned int i = 0; i < n_entries; ++i){
    
      track_tree->GetEntry(i);

      double temp_vertex[3];
      double temp_end[3];
      
      int event_id         = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now         = b_time_now->GetLeaf("time_now")->GetValue();
      temp_vertex[0]       = b_vertex->GetLeaf("tr_vertex")->GetValue(0);
      temp_vertex[1]       = b_vertex->GetLeaf("tr_vertex")->GetValue(1);
      temp_vertex[2]       = b_vertex->GetLeaf("tr_vertex")->GetValue(2);
      temp_end[0]          = b_end->GetLeaf("tr_end")->GetValue(0);
      temp_end[1]          = b_end->GetLeaf("tr_end")->GetValue(1);
      temp_end[2]          = b_end->GetLeaf("tr_end")->GetValue(2);
      float pida           = b_pida->GetLeaf("tr_pida")->GetValue();
      float chi2_mu        = b_chi2_mu->GetLeaf("tr_chi2_mu")->GetValue();
      float chi2_pi        = b_chi2_pi->GetLeaf("tr_chi2_pi")->GetValue();
      float chi2_pr        = b_chi2_pr->GetLeaf("tr_chi2_pr")->GetValue();
      float chi2_ka        = b_chi2_ka->GetLeaf("tr_chi2_ka")->GetValue();
      float length         = b_length->GetLeaf("tr_length")->GetValue();
      float kinetic_energy = b_kinetic_energy->GetLeaf("tr_kinetic_energy")->GetValue();

      TVector3 vertex(temp_vertex);
      TVector3 end(temp_end);

      if(event_id == unique_event.first && time_now == unique_event.second)
        track_list.push_back(Track(pida, chi2_mu, chi2_pi, chi2_pr, chi2_ka, length, kinetic_energy, vertex, end));
    } 
  
  }
  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetShowerList(TTree *shower_tree, const std::pair<int, int> &unique_event, ShowerList &shower_list){
  
    TBranch *b_event_id   = shower_tree->GetBranch("event_id");
    TBranch *b_time_now   = shower_tree->GetBranch("time_now");
    TBranch *b_vertex     = shower_tree->GetBranch("sh_start");
    TBranch *b_direction  = shower_tree->GetBranch("sh_direction");
    TBranch *b_open_angle = shower_tree->GetBranch("sh_open_angle");
    TBranch *b_length     = shower_tree->GetBranch("sh_length");
    
    unsigned int n_entries = shower_tree->GetEntries();

    for(unsigned int i = 0; i < n_entries; ++i){
    
      shower_tree->GetEntry(i);

      double temp_vertex[3];
      double temp_direction[3];
      
      int event_id      = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now      = b_time_now->GetLeaf("time_now")->GetValue();
      temp_vertex[0]    = b_vertex->GetLeaf("sh_start")->GetValue(0);
      temp_vertex[1]    = b_vertex->GetLeaf("sh_start")->GetValue(1);
      temp_vertex[2]    = b_vertex->GetLeaf("sh_start")->GetValue(2);
      temp_direction[0] = b_direction->GetLeaf("sh_direction")->GetValue(0);
      temp_direction[1] = b_direction->GetLeaf("sh_direction")->GetValue(1);
      temp_direction[2] = b_direction->GetLeaf("sh_direction")->GetValue(2);
      float open_angle  = b_open_angle->GetLeaf("sh_open_angle")->GetValue();
      float length      = b_length->GetLeaf("sh_length")->GetValue();
 
      TVector3 vertex(temp_vertex);
      TVector3 direction(temp_direction);

      if(event_id == unique_event.first && time_now == unique_event.second)
        shower_list.push_back(Shower(vertex, direction, open_angle, length));
    } 
  
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetMCParticleList(TTree *mcparticle_tree, const std::pair<int, int> &unique_event, ParticleList &mcparticle_list){
    
    TBranch *b_event_id = mcparticle_tree->GetBranch("event_id");
    TBranch *b_time_now = mcparticle_tree->GetBranch("time_now");
    TBranch *b_pdgcode  = mcparticle_tree->GetBranch("p_pdgcode");
    TBranch *b_mass     = mcparticle_tree->GetBranch("p_mass");
    TBranch *b_energy   = mcparticle_tree->GetBranch("p_energy");
    TBranch *b_vertex   = mcparticle_tree->GetBranch("p_vertex");
    TBranch *b_end      = mcparticle_tree->GetBranch("p_end");
    TBranch *b_momentum = mcparticle_tree->GetBranch("p_momentum");
    
    unsigned int n_entries = mcparticle_tree->GetEntries();

    for(unsigned int i = 0; i < n_entries; ++i){
    
      mcparticle_tree->GetEntry(i);

      double temp_vertex[3];
      double temp_end[3];
      double temp_momentum[3];
      
      int event_id          = b_event_id->GetLeaf("event_id")->GetValue();
      int time_now          = b_time_now->GetLeaf("time_now")->GetValue();
      int pdgcode           = b_pdgcode->GetLeaf("p_pdgcode")->GetValue();
      float mass            = b_mass->GetLeaf("p_mass")->GetValue();
      float energy          = b_energy->GetLeaf("p_energy")->GetValue();
      temp_vertex[0]        = b_vertex->GetLeaf("p_vertex")->GetValue(0);
      temp_vertex[1]        = b_vertex->GetLeaf("p_vertex")->GetValue(1);
      temp_vertex[2]        = b_vertex->GetLeaf("p_vertex")->GetValue(2);
      temp_end[0]           = b_end->GetLeaf("p_end")->GetValue(0);
      temp_end[1]           = b_end->GetLeaf("p_end")->GetValue(1);
      temp_end[2]           = b_end->GetLeaf("p_end")->GetValue(2);
      temp_momentum[0]      = b_momentum->GetLeaf("p_momentum")->GetValue(0);
      temp_momentum[1]      = b_momentum->GetLeaf("p_momentum")->GetValue(1);
      temp_momentum[2]      = b_momentum->GetLeaf("p_momentum")->GetValue(2);
 
      TVector3 vertex(temp_vertex);
      TVector3 end(temp_end);
      TVector3 momentum(temp_momentum);

      if(event_id == unique_event.first && time_now == unique_event.second)
        mcparticle_list.push_back(Particle(pdgcode, mass, energy, vertex, end, momentum));
      
      /*
      std::cout << "-------------------------------------------------------------------------" << std::endl;
      std::cout << " Event   : " << event_id  << std::endl;
      std::cout << " PdgCode : " << pdgcode   << std::endl;
      std::cout << " Vertex  : " << vertex[0] << std::endl;
      std::cout << " End     : " << end[0]    << std::endl;
      */
    }
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetRecoParticleFromTrack(const TrackList &track_list, ParticleList &recoparticle_list){
 
    // Muon candidates 
    std::vector<unsigned int> mu_candidates;

    // Loop over track list
    for(unsigned int id = 0; id < track_list.size(); ++id){
   
      const Track &track(track_list[id]);

      // Use PIDA values to find pions, protons and kaons
      int pida_pdg = EventSelectionTool::GetPdgByPIDA(track);

      if(pida_pdg == 13) 
        mu_candidates.push_back(id);
      else if(pida_pdg == 211 || pida_pdg == 321 || pida_pdg == 2212) 
        recoparticle_list.push_back(Particle(pida_pdg, track.m_kinetic_energy, track.m_vertex, track.m_end));
      else
        recoparticle_list.push_back(Particle(EventSelectionTool::GetPdgByChi2(track), track.m_kinetic_energy, track.m_vertex, track.m_end)); 
    }

    if(mu_candidates.size() == 0) return;
    if(mu_candidates.size() == 1) {
      const Track &muon(track_list[mu_candidates[0]]);
      recoparticle_list.push_back(Particle(13, muon.m_kinetic_energy, muon.m_vertex, muon.m_end));
      return;
    }
   
    // If more than one muon candidate exists
    bool foundOneMuon(false);
    unsigned int muonID = std::numeric_limits<unsigned int>::max();

    for(unsigned int i = 0; i < mu_candidates.size(); ++i){
    
      unsigned int id = mu_candidates[i];
      const Track &candidate(track_list[id]);

      if(candidate.m_pida*candidate.m_chi2_mu >= 5 && candidate.m_pida*candidate.m_chi2_mu < 9) {
      
        if(!foundOneMuon) {
          muonID = id;
          foundOneMuon = true;
        }
        else{
          foundOneMuon = false;
          break;
        }
      }
    }
    if(!foundOneMuon) {
    
      float min_chi2 = std::numeric_limits<float>::max();
      unsigned int best_id  = std::numeric_limits<unsigned int>::max();

      for(unsigned int i = 0; i < mu_candidates.size(); ++i){

        unsigned int id = mu_candidates[i];
        const Track &candidate(track_list[id]);

        if(candidate.m_chi2_mu < min_chi2) {
          best_id  = id;
          min_chi2 = candidate.m_chi2_mu; 
        }
      }
      muonID = best_id;
    } 
    
    const Track &muon(track_list[muonID]);
    recoparticle_list.push_back(Particle(13, muon.m_kinetic_energy, muon.m_vertex, muon.m_end));

    for(unsigned int id = 0; id < mu_candidates.size(); ++id){
      
      if(id == muonID) continue;
      
      const Track &track(track_list[id]);
      recoparticle_list.push_back(Particle(EventSelectionTool::GetPdgByChi2(track), track.m_kinetic_energy, track.m_vertex, track.m_end)); 
    } 
  }

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetRecoParticleFromShower(const ShowerList &shower_list, const TVector3 &reco_vertex, ParticleList &recoparticle_list){
 
    // To start with, 
    //  If an event has 2 showers
    //  Push back a pi0 at that point
    if(shower_list.size() == 2) recoparticle_list.push_back(Particle(111,reco_vertex,reco_vertex));
 
  }

  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetPdgByChi2(const Track &track){

    // Push the chi2 values onto a vector to find the minimum
    // NOT MUON
    std::map<int, float> chi2_map;

    chi2_map.insert(std::map<int, float>::value_type(211,  track.m_chi2_pi));
    chi2_map.insert(std::map<int, float>::value_type(321,  track.m_chi2_ka));
    chi2_map.insert(std::map<int, float>::value_type(2212, track.m_chi2_pr));

    float min_chi2 = std::numeric_limits<float>::max();
    int best_pdg   = std::numeric_limits<int>::max();

    for(std::map<int, float>::const_iterator it = chi2_map.begin(); it != chi2_map.end(); ++it){
    
      if(it->second < min_chi2){
        min_chi2 = it->second; 
        best_pdg = it->first;
      }
    } 
    return best_pdg;
  }
  
  //------------------------------------------------------------------------------------------ 
  
  int EventSelectionTool::GetPdgByPIDA(const Track &track){

    // Muon
    if(track.m_pida >= 5 && track.m_pida < 9) return 13;

    //Pion
    if(track.m_pida >= 9 && track.m_pida < 13) return 211;
  
    //Kaon
    if(track.m_pida >= 13 && track.m_pida < 13.5) return 321;
  
    //Proton
    if(track.m_pida >= 13.5 && track.m_pida < 21) return 2212;

    return std::numeric_limits<int>::max();
  
  }
  
  //------------------------------------------------------------------------------------------ 
  
  EventSelectionTool::Track::Track(const float &pida, const float &chi2_mu, const float &chi2_pi, const float &chi2_pr, const float &chi2_ka, const float &length, const float &kinetic_energy, const TVector3 &vertex, const TVector3 &end) :
    m_pida(pida),
    m_chi2_mu(chi2_mu),
    m_chi2_pi(chi2_pi),
    m_chi2_pr(chi2_pr),
    m_chi2_ka(chi2_ka),
    m_length(length),
    m_kinetic_energy(kinetic_energy),
    m_vertex(vertex),
    m_end(end) {}

  //------------------------------------------------------------------------------------------ 
  
  EventSelectionTool::Shower::Shower(const TVector3 &vertex, const TVector3 &direction, const float &open_angle, const float &length) :
    m_vertex(vertex),
    m_direction(direction),
    m_open_angle(open_angle),
    m_length(length) {}

} // namespace: selection
