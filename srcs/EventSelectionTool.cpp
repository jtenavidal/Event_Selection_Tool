#include "../include/EventSelectionTool.h"

namespace selection{
/* 
  void EventSelectionTool::LoadEventList(const std::string &file_name, EventList &event_list);

    // Construct event from defined variables above and fill event list
    //Event e();
    //event_list.push_back(e);
  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetUniqueEventList(const TTree &event_tree, UniqueEventIdList &unique_event_list);

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetTrackList(const TTree &track_tree, const UniqueEventIdList &unique_event_list, TrackList &track_list);

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetShowerList(const TTree &shower_tree, const UniqueEventIdList &unique_event_list, ShowerList &shower_list);

  //------------------------------------------------------------------------------------------ 
  
  void EventSelectionTool::GetMCParticleList(const TTree &mcparticle_tree, const UniqueEventIdList &unique_event_list, ParticleList &mcparticle_list);
*/
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
      unsigned int   best_id  = std::numeric_limits<unsigned int>::max();

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
    
    const Track &muon(track_list[mu_candidates[muonID]]);
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
  
  EventSelectionTool::Shower::Shower(const TVector3 &vertex, const TVector3 &end, const TVector3 &direction, const float &open_angle) :
    m_vertex(vertex),
    m_end(end),
    m_direction(direction),
    m_open_angle(open_angle) {}

} // namespace: selection
