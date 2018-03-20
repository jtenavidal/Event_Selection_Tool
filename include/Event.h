#ifndef EVENT_H
#define EVENT_H

#include "Particle.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include "TH1.h"
#include "TF1.h"

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
      Event(const ParticleList &mc_particles, const ParticleList &reco_particles, const unsigned int nuance, const int neutrino_pdg, const unsigned int charged_pi, const unsigned int neutral_pi, const bool is_cc, const TVector3 &mc_vertex, const TVector3 &reco_vertex, const float neutrino_energy);

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
       * @brief  Get the neutrino pdg code in the event
       */
      int GetNeutrinoPdgCode() const;

      /**
       * @brief  Get the number of charged pions
       */
      int GetNChargedPions() const;
      
      /**
       * @brief  Get the number of neutral pions
       */
      int GetNNeutralPions() const;
      
      /**
       * @brief  Get the physical process
       */
      int GetPhysicalProcess() const;
      
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
      ParticleMatrix CountLength_topology( const TopologyMap & signal_map_topology , ParticleMatrix & Count_L ,ParticleMatrix & Count_2L );


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
       * @brief Returns MC Energy of a given particle     
       * @param pdg                        
       */
      float GetMCEnergyWithPdg(const int pdg) const;

      /**                                                              
       * @brief Returns MC Energy of a given particle     
       * @param pdg                        
       */
      float GetRecoEnergyWithPdg(const int pdg) const;

      /**                                                              
       * @brief Returns MC Kinetic Energy of a given particle     
       * @param pdg                        
       */
      float GetMCKineticEnergyWithPdg(const int pdg) const;

      /**                                                              
       * @brief Returns MC Kinetic Energy of a given particle     
       * @param pdg                        
       */
      float GetRecoKineticEnergyWithPdg(const int pdg) const;

      /**                                                              
       * @brief Returns the momentum module of a given particle     
       * @param pdg                        
       */
      float GetMCModuleMomentumWithPdg(const int pdg) const;

      /**                                                              
       * @brief Returns the momentum module of a given particle     
       * @param pdg                        
       */
      float GetRecoModuleMomentumWithPdg(const int pdg) const;
      /**                                                              
       * @brief Returns MC Energy of the particle with longest track for cc1p                         
       */
      float GetMCEnergyLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with longest track for cc1p    
       */
      float GetRecoEnergyLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with the second longest track for cc1p    
       */
      float GetMCEnergySecondLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with the second longest track for cc1p    
       */
      float GetRecoEnergySecondLongest_cc1pi( ) const;

      float GetMCKineticEnergyLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns Reco Energy  of the particle with longest track for cc1p    
       */
      float GetRecoKineticEnergyLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns MC kinetic Energy  of the particle with the second longest track for cc1p    
       */
      float GetMCKineticEnergySecondLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns Reco kinetic Energy  of the particle with the second longest track for cc1p    
       */
      float GetRecoKineticEnergySecondLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns MC momentum module  of the particle with longest track for cc1p    
       */
      float GetMCModuleMomentumLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns Reco momentum module  of the particle with longest track for cc1p    
       */
      float GetRecoModuleMomentumLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns momentum module  of the particle with the second longest track for cc1p    
       */
      float GetMCModuleMomentumSecondLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns momentum module of the particle with the second longest track for cc1p    
       */
      float GetRecoModuleMomentumSecondLongest_cc1pi( ) const;
      /**                                                              
       * @brief Returns the energy of the delta particle produced in a resonance                                                     
       */
      float GetMCDeltaEnergy( ) const;

      /**                                                              
       * @brief Returns the energy of the delta particle produced in a resonance                                                     
       */
      float GetRecoDeltaEnergy( ) const;

      /**
       * @brief  Get the true neutrino energy
       */
      float GetTrueNuEnergy() const;
      /**
       * @brief  Get the true neutrino energy for CC 0pi interactions
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      float GetMCCC0piNeutrinoEnergy() const;

      /**
       * @brief  Get the true neutrino energy for CC 1pi interactions
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      float GetMCCC1piNeutrinoEnergy() const;

      /**
       * @brief  Get the reconstructed neutrino energy for CC 0pi interactions
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      float GetRecoCC0piNeutrinoEnergy() const;

      /**
       * @brief  Get the reconstructed neutrino energy for CC 1pi interactions
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      float GetRecoCC1piNeutrinoEnergy() const;
      /**
       * @brief  Get the reconstructed neutrino energy for CC 1pi interactions (METHOD 2)
       *
       * @param  track muon track to use in the calculation
       *
       * @return float reconstructed neutrino energy 
       */
      float GetRecoCC1piNeutrinoEnergyMethod2(  ) const;
      /**
       * @brief  Get the most energetic reconstructed particle
       *
       * @return Particle most energetic reco
       */
      Particle GetMostEnergeticRecoParticle() const;

      /**
       * @brief  Get the most energetic true particle
       *
       * @return Particle most energetic true
       */
      Particle GetMostEnergeticTrueParticle() const;

      /**
       * @brief Calculates the Efficiency, Purity, Background Rejection
       * Parameters for Efficiency calculation ( MC, TReco and Reco ) for a given topology : 
       * 0-> No muon , 
       * 1 -> CCinclusive,
       * 2-> CC0pi, 
       * 3-> CC1pi+/-,
       * 4-> CC1pi0
       **/
      double Efficiency( const std::vector< double > & CountMC, const std::vector< double > & CountTReco, const std::vector< double > & CountReco, const TopologyMap &topology  ) const;


      /**
       * @brief Save Topology Matrix into a file
       * BackGround Study : topology mis identification table 
       * 0-> No muon , 
       * 1 -> CCinclusive,
       * 2-> CC0pi, 
       * 3-> CC1pi+/-,
       * 4-> CC1pi0
       **/
      void SaveTopologyMatrix( const ParticleMatrix & Count_MC_Topology, const ParticleMatrix & Count_TReco_Topology, const ParticleMatrix & Count_Reco_Topology ) const ;

      /**
       * @brief Saves the event information in a file ( types of particles in the event 
       * and topology for True and Selected
       **/
      void EventInformationParticles( const std::string name, const int event_number ) const;

      /**
       * @brief Saves the event characteristics in a file ( length , angle and kinetic energy ) 
       * for the selected topology
       **/
      void EventProperties( const TopologyMap &topology, std::string event_file, const int event_number ) const;

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
      /**                                                              
       * @brief Returns the energy of a particle given a pdg                        
       * @param pdg, particle_list                                     
       */
      float EnergyWithPdg(const int pdg, const ParticleList &particle_list) const;
      /**                                                              
       * @brief Returns the Kinetic energy of a particle given a pdg                        
       * @param pdg, particle_list                                     
       */
      float KineticEnergyWithPdg(const int pdg, const ParticleList &particle_list) const;

      /**                                                              
       * @brief Returns the energy of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      float GetEnergyLongest_cc1pi( const ParticleList & particle_list) const;
      /**                                                              
       * @brief Returns the energy of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      float GetEnergySecondLongest_cc1pi( const ParticleList & particle_list) const;
      /**                                                              
       * @brief Returns the kinetic energy of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      float GetKineticEnergyLongest_cc1pi( const ParticleList & particle_list) const;
      /**                                                              
       * @brief Returns the kinetic energy of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      float GetKineticEnergySecondLongest_cc1pi( const ParticleList & particle_list) const;
      /**                                                              
       * @brief Returns the momentum module of a particle with longest track lenght                        
       * @param pdg, particle_list                                     
       */
      float GetModuleMomentumLongest_cc1pi( const ParticleList & particle_list) const;
      /**                                                              
       * @brief Returns the momentum module of a particle with the second longest track lenght                 
       * @param pdg, particle_list                                     
       */
      float GetModuleMomentumSecondLongest_cc1pi( const ParticleList & particle_list) const;
      /**                                                              
       * @brief Returns the energy of the delta particle produced in a resonance                 
       * @param particle_list                                     
       */
      float GetDeltaEnergy( const ParticleList & particle_list ) const;
      /**                                                                                                                                                   * @brief Returns the energy of the delta particle produced in a resonance with a proton in the final state
       * @param particle_list                                                                                                                    
       */
      float GetDeltaEnergy_p( const ParticleList & particle_list ) const;

      /**                                                              
       * @brief Returns the momentum module of a particle given a pdg                        
       * @param pdg, particle_list                                     
       */
      float ModuleMomentumWithPdg(const int pdg, const ParticleList &particle_list) const;

      /**
       * @brief  Get the most energetic particle
       *
       * @return Particle most energetic
       */
      Particle GetMostEnergeticParticle(const ParticleList &particle_list) const;
      /**
       * @brief  Get the neutrino energy for CC0pi
       */      
      float GetCC0piNeutrinoEnergy( const ParticleList & particle_list ) const;
      /**
       * @brief  Get the neutrino energy for CC1pi
       */
      float GetCC1piNeutrinoEnergy( const ParticleList & particle_list ) const;
      /**
       * @brief  Get the neutrino energy for CC1pi (METHOD 2)
       */
      float GetCC1piNeutrinoEnergyMethod2( const ParticleList & particle_list ) const;

      // Member variables
      ParticleList       m_mc_particles;       ///< vector of Monte Carlo particles
      ParticleList       m_reco_particles;     ///< vector of reconstructed particles
      unsigned int       m_nuance;             ///< Nuance code/interaction of the event
      int                m_nu_pdg;             ///< Neutrino pdg code of the event
      unsigned int       m_charged_pi;         ///< Number of charged pions in the event
      unsigned int       m_neutral_pi;         ///< Number of neutral pions in the event
      bool               m_is_cc;              ///< whether the event contains and CC or NC interaction
      TVector3           m_reco_vertex;        ///< reconstructed neutrino vertex
      TVector3           m_mc_vertex;          ///< reconstructed neutrino vertex
      float              m_neutrino_energy;    ///< true neutrino energy



  }; // Event
} // selection

#endif
