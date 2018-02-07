#include "../include/Event.h"

namespace selection{
  
  Event::Event(const ParticleList &mc_particles, const ParticleList &reco_particles, const unsigned int nuance, const bool is_cc, const TVector3 &mc_vertex, const TVector3 &reco_vertex) :
    m_mc_particles(mc_particles),
    m_reco_particles(reco_particles),
    m_nuance(nuance),
    m_is_cc(is_cc),
    m_mc_vertex(mc_vertex),
    m_reco_vertex(reco_vertex) {}

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
  //----------------------------------------------------------------------------------------
  void Event::SetTopologies( ) {
    /*                                                               
     * NC topology  ( Topology number 0 )                            
     * TopologyMap signal_map_NC;                                    
     */
    std::vector< int > NC_mu;

    NC_mu.push_back( 13 );
    signal_map_NC.insert( std::make_pair( NC_mu, 0 ) );

    /*                                                                 
     * CC inclusive topology ( Topology number 1 )                     
     * TopologyMap signal_map_cc_inclusive;                            
     */
    std::vector< int > cc_inclusive_mu;

    cc_inclusive_mu.push_back( 13 );
    signal_map_cc_inclusive.insert( std::make_pair( cc_inclusive_mu, 1 ) );

    /*                                                                 
     * CC0pi topology ( Topology number 2 )                            
     * TopologyMap signal_map_cc_0pi;                                  
     */
    std::vector< int > cc_0pi_mu;
    std::vector< int > cc_0pi_pi;

    cc_0pi_mu.push_back( 13 );
    cc_0pi_pi.push_back( 211 );
    cc_0pi_pi.push_back(-211 );
    cc_0pi_pi.push_back( 111 );

    signal_map_cc_0pi.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map_cc_0pi.insert( std::make_pair( cc_0pi_pi, 0 ) );

    /*                                                                 
     * CC1pi+/- topology  ( Topology number 3 )                        
     * TopologyMap signal_map_cc_1pi;                                  
     */
    std::vector< int > cc_1pi_mu;
    std::vector< int > cc_1pi_pi;

    cc_1pi_mu.push_back( 13 );
    cc_1pi_pi.push_back( 211 );
    cc_1pi_pi.push_back(-211 );

    signal_map_cc_1pi.insert( std::make_pair( cc_1pi_mu, 1 ) );
    signal_map_cc_1pi.insert( std::make_pair( cc_1pi_pi, 1 ) );

    /*                                                                 
     * CC1pi0 topology ( Topology number 4 )                           
     * TopologyMap signal_map_cc_pi0;                                  
     */
    std::vector< int > cc_pi0_mu;
    std::vector< int > cc_pi0_pi;

    cc_pi0_mu.push_back( 13 );
    cc_pi0_pi.push_back( 111 );

    signal_map_cc_pi0.insert( std::make_pair( cc_pi0_mu, 1 ) );
    signal_map_cc_pi0.insert( std::make_pair( cc_pi0_pi, 1 ) );


  }
  
  //-------------------------------------------------------------------------------------------------
  void Event::Count_per_Topology( const TopologyMap signal_map_topology, double & Count_MC, double & Count_TReco, double & Count_Reco ){

    if( CheckMCTopology( signal_map_topology ) == 1 ) {
      Count_MC++;
    }
    if( CheckRecoTopology( signal_map_topology ) == 1 ) {
    Count_Reco++;
    }
    if( CheckMCTopology( signal_map_topology ) == 1 && CheckRecoTopology( signal_map_topology ) == 1 ){
      Count_TReco++ ;
    }
  }
  
  //-------------------------------------------------------------------------------------------------
  ParticleMatrix Event::TopologyMatrix( ParticleMatrix & Count_MC_Topology, ParticleMatrix & Count_TReco_Topology, ParticleMatrix & Count_Reco_Topology ) \
  {

    std::vector< TopologyMap > topology_vector (5);
    topology_vector[0] = signal_map_NC;
    topology_vector[1] = signal_map_cc_inclusive;
    topology_vector[2] = signal_map_cc_0pi;
    topology_vector[3] = signal_map_cc_1pi;
    topology_vector[4] = signal_map_cc_pi0;

    for( unsigned int i=0 ; i < 5 ; ++i ){
      for(unsigned int j=0; j < 5 ; ++j ){

	if ( CheckMCTopology( topology_vector[i] ) == 1 ){
	  Count_MC_Topology[i][j]++;
	}
	if ( CheckRecoTopology( topology_vector[i] ) == 1 ){
	  Count_Reco_Topology[i][j]++;
	}
	if ( CheckMCTopology( topology_vector[i] ) == 1 && CheckRecoTopology( topology_vector[j] ) == 1 ){
	  Count_TReco_Topology[i][j]++;
	}

      }

    }
    return Count_TReco_Topology;
  }


  //------------------------------------------------------------------------------------------------


  float Event::GetRecoLengthWithPdg(const int pdg) const{

    return this -> LengthWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCLengthWithPdg(const int pdg) const{

    return this -> LengthWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------

  float Event::LengthWithPdg(const int pdg, const ParticleList &particle_list) const{

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) return particle_list[i].GetLength();
    }
  }

  //-----------------------------------------------------------------------------------------------
  ParticleMatrix Event::CountLength_topology( const TopologyMap & signal_map_topology , ParticleMatrix & Count_L ) {

    if( CheckMCTopology(signal_map_topology) ){

      if( GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 2212 ) && GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 211 ) ){ Count_L[0][0]++ ; }

      if( GetMCLengthWithPdg( 211 ) > GetMCLengthWithPdg( 2212 ) ) {  Count_L[0][1]++; } else  Count_L[0][2]++;

      if( CheckRecoTopology(signal_map_cc_1pi) == 1 ){

	if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 ) ){
	  Count_L[1][0]++ ;}

	if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) ) { Count_L[1][1]++; } else Count_L[1][2]++ ;
      }
    }
    if( CheckRecoTopology( signal_map_cc_1pi ) == 1 ) {

      if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 ) ){Count_L[2][0]++ ;}

      if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) ) { Count_L[2][1]++; }else Count_L[2][2]++;
    }
    return Count_L;
  }

  //------------------------------------------------------------------------------------------------


  float Event::GetRecoCosThetaWithPdg(const int pdg) const{

    return this -> CosThetaWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCCosThetaWithPdg(const int pdg) const{

    return this -> CosThetaWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------

  float Event::CosThetaWithPdg(const int pdg, const ParticleList &particle_list) const{

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) {

	return particle_list[i].GetCosTheta();
      }
    }
  }

    //--------------------------------------------------------------------------------------------                                               
  ParticleMatrix Event::ParticleReconstruction(  ParticleMatrix & Count_MC_ID,  ParticleMatrix & Count_TReco_ID,  ParticleMatrix & Count_Reco_ID,  const int pdg) {

    unsigned int k = 0 ; //default muon                                                                                                        
    if( pdg == 2212 ) { //Proton                                                                                                               
      k=3;
    }else if ( pdg == 211 || pdg == -211 ){ //Pion                                                                                             
      k=1;
    }else if ( pdg == 111 ){// Pi0                                                                                                             
      k=2;
    }
        //PARTICLE RECONSTRUCTION                                                                                                              
    if ( CountMCParticlesWithPdg(pdg) == 1){
      Count_MC_ID[1][k]++;

      if ( CountRecoParticlesWithPdg(pdg) == 1 ){
        Count_TReco_ID[1][k]++;
      }
    } else if  ( CountMCParticlesWithPdg(pdg) == 0 ){
      Count_MC_ID[0][k]++;

      if ( CountRecoParticlesWithPdg(pdg) == 0 ){
        Count_TReco_ID[0][k]++;
      }
    } else if  ( CountMCParticlesWithPdg(pdg) == 2 ){
      Count_MC_ID[2][k]++;
      if ( CountRecoParticlesWithPdg(pdg) == 2 ){
	Count_TReco_ID[2][k]++;
      }
    } else if  ( CountMCParticlesWithPdg(pdg) > 2 ){
      Count_MC_ID[3][k]++;

      if ( CountRecoParticlesWithPdg(pdg) >2 ){
        Count_TReco_ID[3][k]++;
      }
    }
    if ( CountRecoParticlesWithPdg(pdg) == 0 ){
      Count_Reco_ID[0][k]++;

    } else if  ( CountRecoParticlesWithPdg(pdg) == 1 ){
      Count_Reco_ID[1][k]++;

    } else if  ( CountRecoParticlesWithPdg(pdg) == 2 ){
      Count_Reco_ID[2][k]++;

    } else if  ( CountRecoParticlesWithPdg(pdg) > 2 ){
      Count_Reco_ID[3][k]++;
    }

    return Count_MC_ID;
  }

  //-----------------------------------------------------------------------------------------------                                            
                                          

  ParticleMatrix Event::ParticleExChange( ParticleMatrix & Count_ExChange_MC, ParticleMatrix & Count_ExChange_TReco, ParticleMatrix & Count_ExChange_Reco ){

    int pdg[] = {13, 211, 111, 2212} ;

    for( unsigned int i=0 ; i < 4 ; ++i ){
      for(unsigned int j=0; j < 4 ; ++j ){
           if ( i!= j && CountMCParticlesWithPdg( pdg[i] ) == 1 && CountMCParticlesWithPdg( pdg[j] ) == 0 ){
             Count_ExChange_MC[i][j]++;
             if ( CountRecoParticlesWithPdg( pdg[i] ) == 1 && CountRecoParticlesWithPdg( pdg[j] ) == 0 ){
               Count_ExChange_TReco[i][j]++;
             }
             if ( CountRecoParticlesWithPdg( pdg[i] ) == 0 && CountRecoParticlesWithPdg( pdg[j] ) == 1 ){
               Count_ExChange_Reco[i][j]++;
             }
           }
      }
    }

    return Count_ExChange_TReco;
  }




} // Selection
