#include "../include/Event.h"

namespace selection{
  
  Event::Event(const ParticleList &mc_particles, const ParticleList &reco_particles, const unsigned int nuance, const bool is_cc, const TVector3 &mc_vertex, const TVector3 &reco_vertex, const float neutrino_energy) :
    m_mc_particles(mc_particles),
    m_reco_particles(reco_particles),
    m_nuance(nuance),
    m_is_cc(is_cc),
    m_mc_vertex(mc_vertex),
    m_reco_vertex(reco_vertex), 
    m_neutrino_energy(neutrino_energy) {}

  //------------------------------------------------------------------------------------------ 
    
  unsigned int Event::CountMCParticlesWithPdg(const int pdg) const{
 
    return this->CountParticlesWithPdg(pdg, m_mc_particles);

  }

  //------------------------------------------------------------------------------------------ 

  unsigned int Event::CountRecoParticlesWithPdg(const int pdg) const{
    
    return this->CountParticlesWithPdg(pdg, m_reco_particles);

  }

  //------------------------------------------------------------------------------------------ 

  bool Event::CheckMCTopology(const TopologyMap &topology) const{
  
    return this->CheckTopology(topology, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------ 

  bool Event::CheckRecoTopology(const TopologyMap &topology) const{
  
    return this->CheckTopology(topology, m_reco_particles);

  }

  //------------------------------------------------------------------------------------------ 

  ParticleList Event::GetMCParticleList() const{
  
    return m_mc_particles;

  }

  //------------------------------------------------------------------------------------------ 

  ParticleList Event::GetRecoParticleList() const{
  
    return m_reco_particles;

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNuanceCode() const{
  
    return m_nuance;
  
  }

  //------------------------------------------------------------------------------------------ 

  bool Event::GetIsCC() const{
  
    return m_is_cc;
  
  }

  //------------------------------------------------------------------------------------------ 

  TVector3 Event::GetMCNuVertex() const{
  
    return m_mc_vertex;
  
  }

  //------------------------------------------------------------------------------------------ 

  TVector3 Event::GetRecoNuVertex() const{
  
    return m_reco_vertex;
  
  }

  //------------------------------------------------------------------------------------------ 
  
  float Event::GetTrueNuEnergy() const{
  
    return m_neutrino_energy;

  }
  
  //------------------------------------------------------------------------------------------ 

  unsigned int Event::CountParticlesWithPdg(const int pdg, const ParticleList &particle_list) const{
  
    unsigned int particle_counter = 0;

    for(unsigned int i = 0; i < particle_list.size(); ++i) if(particle_list[i].GetPdgCode() == pdg) particle_counter++;

    return particle_counter;
  }

  //------------------------------------------------------------------------------------------ 
  
  bool Event::CheckTopology(const TopologyMap &topology, const ParticleList &particle_list) const{
   
    // Loop over the map
    for( TopologyMap::const_iterator it = topology.begin(); it != topology.end(); ++it ){

      // Define temporary variables for the current map element
      std::vector< int > pdg_codes = it->first; 
      int n_total                  = it->second;

      // Count the number of particles in the current event with the same PDG codes 
      // as given by the chosen topology
      int counter = 0;
     
      // Loop over particles in current event
      for(unsigned int i = 0; i < pdg_codes.size(); ++i){

        counter += this->CountParticlesWithPdg(pdg_codes[i], particle_list);
      }
      
      if(counter != n_total) return false;
    }

    return true;

  }

  float Event::GetCC0piRecoNeutrinoEnergy(const Particle &particle) const{
    
    // The variables from the branches and get the leaves
    float m_n   = 0.93828;   // Nucleon mass, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float reco, e, p, cth;   // track variables
    
    // Vector of z direction
    TVector3 z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;

    // Get the values needed
    e    = particle.GetEnergy();
    p    = particle.GetMomentum().Mag();
    cth  = (1/p) * (particle.GetMomentum()).Dot(z);
    
    reco = ( 1 / ( 1 - ( ( 1 / m_n ) * ( e - p*cth ) ) ) ) * ( e - ( 1 / ( 2 * m_n) ) * m_mu * m_mu );

    return reco;
  
  }


} // Selection
