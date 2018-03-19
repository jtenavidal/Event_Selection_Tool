#include "../include/Event.h"

namespace selection{
  
  Event::Event(const ParticleList &mc_particles, const ParticleList &reco_particles, const unsigned int nuance, const int neutrino_pdg, const unsigned int charged_pi, const unsigned int neutral_pi, const bool is_cc, const TVector3 &mc_vertex, const TVector3 &reco_vertex, const float neutrino_energy, const unsigned int file_number, const unsigned int event_id) :
    m_mc_particles(mc_particles),
    m_reco_particles(reco_particles),
    m_nuance(nuance),
    m_nu_pdg(neutrino_pdg),
    m_charged_pi(charged_pi),
    m_neutral_pi(neutral_pi),
    m_is_cc(is_cc),
    m_mc_vertex(mc_vertex),
    m_reco_vertex(reco_vertex), 
    m_neutrino_energy(neutrino_energy), 
    m_file_number(file_number), 
    m_event_id(event_id) {}

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

  Particle Event::GetMostEnergeticRecoParticle() const{
 
    return this->GetMostEnergeticParticle(m_reco_particles);

  }
  
  //------------------------------------------------------------------------------------------ 

  Particle Event::GetMostEnergeticTrueParticle() const{
 
    return this->GetMostEnergeticParticle(m_mc_particles);

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
  
  int Event::GetNeutrinoPdgCode() const{
  
    return m_nu_pdg;

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNChargedPions() const{
  
    return m_charged_pi;

  }

  //------------------------------------------------------------------------------------------ 
  
  int Event::GetNNeutralPions() const{
  
    return m_neutral_pi;

  }
  //------------------------------------------------------------------------------------------ 
  
  int Event::GetPhysicalProcess() const{

    // QEL
    if(m_nuance == 0 
    || m_nuance == 1001 
    || m_nuance == 1002) return 0;
    // MEC
    else if(m_nuance == 10) return 1;
    // RES
    else if(m_nuance == 1 
         || m_nuance == 1003 
         || m_nuance == 1004
         || m_nuance == 1005
         || m_nuance == 1006
         || m_nuance == 1007
         || m_nuance == 1008
         || m_nuance == 1009
         || m_nuance == 1010) return 2;
    // DIS
    else if(m_nuance == 2
         || m_nuance == 1091) return 3;
    // COH
    else if(m_nuance == 1097) return 4;
    // Non RES 1pi
    else if(m_charged_pi + m_neutral_pi == 1) return 5;
    // Other
    else return 6;
  
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

  unsigned int Event::GetID() const{

    return m_event_id;
  }
  
  //------------------------------------------------------------------------------------------ 

  unsigned int Event::GetFileNumber() const{

    return m_file_number;
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



  float Event::GetCC0piRecoNeutrinoEnergy(const Particle &particle) const{

    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float V     = 0.02950;   // Nucleon removal energy, GeV
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

    reco = (1/(m_n - V - e + p*cth))*((m_n - V)*e - (m_mu*m_mu*0.5) + m_n*V - (V*V*0.5) + (m_p*m_p - m_n*m_n)*0.5);
    
    return reco;
  
  }

  Particle Event::GetMostEnergeticParticle(const ParticleList &particle_list) const{

    float highest_energy   = -std::numeric_limits<float>::max();
    unsigned int energy_id = std::numeric_limits<unsigned int >::max();

    for(unsigned int i = 0; i < particle_list.size(); ++i){
    
      if(!particle_list[i].GetHasCalorimetry()) continue;

      if(particle_list[i].GetEnergy() > highest_energy) energy_id = i;
    }

    return particle_list[energy_id];

  }
} // selection
