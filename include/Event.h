#ifndef EVENT_H
#define EVENT_H

#include "Particle.h"
#include "TVector3.h"
#include <map>

namespace selection{
  
  // Typedef for the map
  typedef std::map< std::vector< int >, int > TopologyMap;
  typedef std::vector<Particle> ParticleList;
  typedef std::vector< vector<double> > ParticleMatrix ;

  /**
   * @brief  Event class
   */
  class Event{

    public :

      /**
       * @brief  Constructor
       *
       * @param  mc_particles list of the MC particle objects in the event
       * @param  reco_particles list of the reconstructed particle objects in the event
       * @param  nuance the nuance code corresponding to the event
       * @param  is_cc is this a charged or neutral current event
       * @param  mc_vertex Monte Carlo neutrino vertex 
       * @param  reco_vertex reconstructed neutrino vertex
       */
      Event(const ParticleList &mc_particles, const ParticleList &reco_particles, const unsigned int nuance, const bool is_cc, const TVector3 &mc_vertex, const TVector3 &reco_vertex);
        
      /**
       * @brief  CountMCParticlesWithPdg
       *
       * @param  pdg pdg code to count
       *
       * @return number of Monte Carlo particles with the given pdg code
       */
      unsigned int CountMCParticlesWithPdg(const int pdg) const;

      /**
       * @brief  CountRecoParticlesWithPdg
       *
       * @param  pdg pdg code to count
       *
       * @return number of reconstructed partices with the given pdg code
       */
      unsigned int CountRecoParticlesWithPdg(const int pdg) const;

      /**
       * @brief  CheckMCTopology
       *
       * @param  topology map of the chosen topology which holds the number of each chosen type of particle to look for
       *
       * @return boolean as to whether the event has the desired Monte Carlo topology
       */
      bool CheckMCTopology(const TopologyMap &topology) const;

      /**
       * @brief  CheckRecoTopology
       *
       * @param  topology map of the chosen topology which holds the number of each chosen type of particle to look for
       *
       * @return boolean as to whether the event has the desired reconstructed topology
       */
      bool CheckRecoTopology(const TopologyMap &topology) const;

      /**
       * @brief  Get the list of MC particls for this event
       */
      ParticleList GetMCParticleList() const;

      /**
       * @brief  Get the list of reconstructed particls for this event
       */
      ParticleList GetRecoParticleList() const;
      
      /**
       * @brief  Get the nuance code - interaction of the event
       */
      int GetNuanceCode() const;

      /**
       * @brief  Get if the event is CC or NC
       */
      bool GetIsCC() const;

      /**
       * @brief  Get the Monte Carlo neutrino vertex position
       */
      TVector3 GetMCNuVertex() const;

      /**
       * @brief  Get the reconstructed neutrino vertex position
       */
      TVector3 GetRecoNuVertex() const;

      /**                                                              
       * @brief Define topologies (Event objects)                      
       * @param signal_map_NC, signal_map_cc_inclusive_signal_map_cc_0pi, signal_map_cc_1pi, signal_map_cc_pi0;                              
       */

      TopologyMap signal_map_NC;
      TopologyMap signal_map_cc_inclusive;
      TopologyMap signal_map_cc_0pi;
      TopologyMap signal_map_cc_1pi;
      TopologyMap signal_map_cc_pi0;
      void SetTopologies( );

      /**                                                              
       * @brief Gives the number of MC, Reco and Coincidences for a given topology                                                           
       * @param signal_map_topology, Count_MC, Count_TReco, Count_Reco     
       */
      void Count_per_Topology ( const TopologyMap signal_map_topology, double & Count_MC, double & Count_TReco, double & Count_Reco ) ;

      /**                                                              
       * @brief Obtains the Topology matrix for a specific set of even\
ts                                                                     
       * @param Count_MC_Topology, Count_Reco_Topology, Count_Topology\
 matrix                                                                
       */

      ParticleMatrix TopologyMatrix( ParticleMatrix & Count_MC_Topology, ParticleMatrix & Count_TReco_Topology, ParticleMatrix & Count_Reco_Topology );

      /**                                                              
       * @brief Returns the true track length for a given particle     
       * @param pdg                                                    
       */
      float GetMCLengthWithPdg(const int pdg) const;

      /**                                                              
       * @brief Returns the reconstructed track length for a given particle                                                                  
       * @param pdg                                                    
       */
      float GetRecoLengthWithPdg(const int pdg) const;

      /**                                                              
       *  @brief Gives the number of times the Muon has the longest track and the pion or proton have                                         
       *  the second longest trak.                                     
       *  @param signal_map_topology                                   
       */
      ParticleMatrix CountLength_topology( const TopologyMap & signal_map_topology , ParticleMatrix & Count_L );


      /**                                                              
       * @brief Returns the cos(theta) for a MC event                  
       * @param pdg                                                    
       */
      float GetMCCosThetaWithPdg(const int pdg) const ;

      /**                                                              
       * @brief Returns the cos(theta) for a reconstructed event       
       * @param pdg                                                    
       */
      float GetRecoCosThetaWithPdg(const int pdg) const;

      /**                                                              
       * @brief Returns number of reconstructed particles per ID ( MC , TReco and Reco )      
       * @param Count_MC_ID, Count_TReco_ID, Count_Reco_ID, pdg                                           
       */
      ParticleMatrix ParticleReconstruction(  ParticleMatrix & Count_MC_ID,  ParticleMatrix & Count_TReco_ID,  ParticleMatrix & Count_Reco_ID,  const int pdg);

      /**                                                              
       * @brief Returns number of times two particles have been exchanged, i.e. miss identified.       
       * @param Count_ExChange_MC, Count_ExChange_TReco, Count_ExChange_Reco                                          
       */
      ParticleMatrix ParticleExChange( ParticleMatrix & Count_ExChange_MC, ParticleMatrix & Count_ExChange_TReco, ParticleMatrix & Count_ExChange_Reco );

    private : 

      /**
       * @brief  CountParticlesWithPdg
       *
       * @param  pdg pdg code to count
       *
       * @return number of partices with the given pdg code
       */
      unsigned int CountParticlesWithPdg(const int pdg, const ParticleList &particle_list) const;
      
      /**
       * @brief  CheckTopology
       *
       * @param  topology map of the chosen topology which holds the number of each chosen type of particle to look for
       *
       * @return boolean as to whether the event has the desired topology
       */
      bool CheckTopology(const TopologyMap &topology, const ParticleList &particle_list) const;


      /**                                                              
       * @brief Returns the track length for a given particle          
       * @param pdg, particle_list                                     
       */
      float LengthWithPdg(const int pdg, const ParticleList &particle_list) const;

      /**                                                              
       * @brief Returns the cos(theta) of the track direction regarding u_z                                                                  
       * @param pdg, particle_list                                     
       */
      float CosThetaWithPdg(const int pdg, const ParticleList &particle_list) const;



      // Member variables
      ParticleList       m_mc_particles;       ///< vector of Monte Carlo particles
      ParticleList       m_reco_particles;     ///< vector of reconstructed particles
      unsigned int       m_nuance;             ///< Nuance code/interaction of the event
      bool               m_is_cc;              ///< whether the event contains and CC or NC interaction
      TVector3           m_reco_vertex;        ///< reconstructed neutrino vertex
      TVector3           m_mc_vertex;          ///< reconstructed neutrino vertex


  }; // Event
} // Selection

#endif
