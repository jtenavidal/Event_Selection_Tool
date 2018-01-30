// Hel  per class 
//    Static member functions only
//
//    No constructor
//    No member data
//
#include "EventSelectionTool.h"

namespace selection{
 
  static void EventSelectionTool::LoadEventList(const std::string &file_name, EventList &event_list);

    /*
    *
    // Construct event from defined variables above and fill event list
    Event e();
    event_list.push_back(e);
    *
    */
  //------------------------------------------------------------------------------------------ 
  
  static void EventSelectionTool::GetUniqueEventList(const TTree &event_tree, UniqueEventIdList &unique_event_list);

  //------------------------------------------------------------------------------------------ 
  
  static void EventSelectionTool::GetTrackList(const TTree &track_tree, const UniqueEventIdList &unique_event_list, TrackList &track_list);

  //------------------------------------------------------------------------------------------ 
  
  static void EventSelectionTool::GetShowerList(const TTree &shower_tree, const UniqueEventIdList &unique_event_list, ShowerList &shower_list);

  //------------------------------------------------------------------------------------------ 
  
  static void EventSelectionTool::GetMCParticleList(const TTree &mcparticle_tree, const UniqueEventIdList &unique_event_list, ParticleList &mcparticle_list);

  //------------------------------------------------------------------------------------------ 
  
  static void EventSelectionTool::GetRecoParticleFromTrack(const TrackList &track_list, ParticleList &recoparticle_list);

  //------------------------------------------------------------------------------------------ 
  
  static void EventSelectionTool::GetRecoParticleFromShower(const ShowerList &shower_list, ParticleList &recoparticle_list);

  //------------------------------------------------------------------------------------------ 
  
  EventSelectionTool::Track::Track(const float &pida, const float &chi2_mu, const float &chi2_pi, const float &chi2_pr, const float &chi_ka, const float &length, const float &kinetic_energy, const TVector3 &vertex, const TVector3 &end) :
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
