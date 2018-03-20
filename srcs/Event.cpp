#include "../include/Event.h"
#include "../include/EventSelectionTool.h"
namespace selection{
  
  Event::Event(const ParticleList &mc_particles, const ParticleList &reco_particles, const unsigned int nuance, const int neutrino_pdg, const unsigned int charged_pi, const unsigned int neutral_pi, const bool is_cc, const TVector3 &mc_vertex, const TVector3 &reco_vertex, const float neutrino_energy) :
    m_mc_particles(mc_particles),
    m_reco_particles(reco_particles),
    m_nuance(nuance),
    m_nu_pdg(neutrino_pdg),
    m_charged_pi(charged_pi),
    m_neutral_pi(neutral_pi),
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

    if( CheckMCTopology( signal_map_topology ) == 1 ) { Count_MC++; }
    if( CheckRecoTopology( signal_map_topology ) == 1  && GetRecoCC0piNeutrinoEnergy()>0 &&  GetRecoCC1piNeutrinoEnergy()>0) { 
       Count_Reco++; } 
    if( CheckMCTopology( signal_map_topology ) == 1 && CheckRecoTopology( signal_map_topology ) == 1 && GetRecoCC0piNeutrinoEnergy()>0 &&  GetRecoCC1piNeutrinoEnergy()>0){
       Count_TReco++ ;  }
  }

  //-------------------------------------------------------------------------------------------------
  ParticleMatrix Event::TopologyMatrix( ParticleMatrix & Count_MC_Topology, ParticleMatrix & Count_TReco_Topology, ParticleMatrix & Count_Reco_Topology ) {

    std::vector< TopologyMap > topology_vector (5);
    topology_vector[0] = signal_map_NC;
    topology_vector[1] = signal_map_cc_inclusive;
    topology_vector[2] = signal_map_cc_0pi;
    topology_vector[3] = signal_map_cc_1pi;
    topology_vector[4] = signal_map_cc_pi0;

    for( unsigned int i=0 ; i < 5 ; ++i ){
      for(unsigned int j=0; j < 5 ; ++j ){

	if ( CheckMCTopology( topology_vector[i] ) == 1 ){ Count_MC_Topology[i][j]++; }
	if ( CheckRecoTopology( topology_vector[i] ) == 1 && GetRecoCC0piNeutrinoEnergy()>0 &&  GetRecoCC1piNeutrinoEnergy()>0){ 
	  Count_Reco_Topology[i][j]++; }
	if ( CheckMCTopology( topology_vector[i] ) == 1 && CheckRecoTopology( topology_vector[j] ) == 1 && GetRecoCC0piNeutrinoEnergy()>0 &&  GetRecoCC1piNeutrinoEnergy()>0 ){ 
	  Count_TReco_Topology[i][j]++; }
      }
    }
    return Count_TReco_Topology;
  }
 //------------------------------------------------------------------------------------------------


  float Event::GetRecoLengthWithPdg(const int pdg) const{
    return this -> LengthWithPdg( pdg, m_reco_particles );
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCLengthWithPdg(const int pdg) const{
    return this -> LengthWithPdg( pdg, m_mc_particles );
  }

  //------------------------------------------------------------------------------------------------

  float Event::LengthWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg ){ return particle_list[i].GetLength();}
    }
    return 0.;
  }

 //-----------------------------------------------------------------------------------------------
  ParticleMatrix Event::CountLength_topology( const TopologyMap & signal_map_topology , ParticleMatrix & Count_L , ParticleMatrix & Count_2L ) {

    if( CheckMCTopology(signal_map_topology)==1 ){
      //LONGEST
      if( GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 2212 ) && GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 211 ) ) {Count_L[0][0]++ ;}
      else if( GetMCLengthWithPdg( 211 ) > GetMCLengthWithPdg( 2212 ) && GetMCLengthWithPdg( 211 ) > GetMCLengthWithPdg( 13 ) ) { Count_L[0][1]++; } 
      else if( GetMCLengthWithPdg( 2212 ) > GetMCLengthWithPdg( 211 ) && GetMCLengthWithPdg( 2212 ) > GetMCLengthWithPdg( 13 ) ) { Count_L[0][2]++;}

      if( CheckRecoTopology(signal_map_cc_1pi) == 1 ){

	if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 ) ){ Count_L[1][0]++ ;}
	else if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 13 ) ) { Count_L[1][1]++; } 
	else if( GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 211 ) && GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 13 ) ) { Count_L[1][2]++;}
      }
      if( CheckRecoTopology( signal_map_cc_1pi ) == 1 ) {
	if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 ) ){ Count_L[2][0]++ ;}
	else if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 13 ) ) { Count_L[2][1]++; } 
	else if( GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 211 ) && GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 13 ) ) { Count_L[2][2]++;}

	}
      //SECOND LONGEST
        if( GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 211 ) && GetMCLengthWithPdg( 211 ) > GetMCLengthWithPdg( 2212 ) ) {Count_2L[0][1]++ ;
      } else if( GetMCLengthWithPdg( 2212 ) > GetMCLengthWithPdg( 211 ) && GetMCLengthWithPdg( 211 ) > GetMCLengthWithPdg( 13 ) ) {Count_2L[0][1]++ ;
      } else if( GetMCLengthWithPdg( 211 ) > GetMCLengthWithPdg( 13 ) && GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 2212 )  ) { Count_2L[0][0]++;
      } else if( GetMCLengthWithPdg( 2212 ) > GetMCLengthWithPdg( 13 ) && GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 211 )  ) { Count_2L[0][0]++; 
	} else if( GetMCLengthWithPdg( 13 ) > GetMCLengthWithPdg( 2212 ) && GetMCLengthWithPdg( 2212 ) > GetMCLengthWithPdg( 211 )  ) { Count_2L[0][2]++; std::cout<<"Hello"<<std::endl;
	} else if( GetMCLengthWithPdg( 211 ) > GetMCLengthWithPdg( 2212 ) && GetMCLengthWithPdg( 2212 ) > GetMCLengthWithPdg( 13 )  ) { Count_2L[0][2]++; std::cout<<"Hello"<<std::endl; }
	
	if( CheckRecoTopology(signal_map_cc_1pi) == 1 ){
	if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 ) && GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) ) {Count_2L[1][1]++ ;
	} else if( GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 211 ) && GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 13 ) ) {Count_2L[1][1]++ ;
	} else if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 13 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 )  ) { Count_2L[1][0]++;
	} else if( GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 13 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 )  ) { Count_2L[1][0]++; 
	} else if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 211 )  ) { Count_2L[1][2]++; 
	} else if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 13 )  ) { Count_2L[1][2]++; }

      }
      if( CheckRecoTopology( signal_map_cc_1pi ) == 1 ) {
	if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 ) && GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) ) {Count_2L[2][1]++ ;
	} else if( GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 211 ) && GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 13 ) ) {Count_2L[2][1]++ ;
	} else if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 13 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 )  ) { Count_2L[2][0]++;
	} else if( GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 13 ) && GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 211 )  ) { Count_2L[2][0]++; 
	} else if( GetRecoLengthWithPdg( 13 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 211 )  ) { Count_2L[2][2]++; 
	} else if( GetRecoLengthWithPdg( 211 ) > GetRecoLengthWithPdg( 2212 ) && GetRecoLengthWithPdg( 2212 ) > GetRecoLengthWithPdg( 13 )  ) { Count_2L[2][2]++; } }
	
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
    return 0.;
  }
 //------------------------------------------------------------------------------------------------


  float Event::GetMCEnergyWithPdg(const int pdg) const{
    return this -> EnergyWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------


  float Event::GetRecoEnergyWithPdg(const int pdg) const{
    return this -> EnergyWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::EnergyWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) return particle_list[i].GetEnergy();
    }
    return 0.;
 }

  //------------------------------------------------------------------------------------------------


  float Event::GetMCKineticEnergyWithPdg(const int pdg) const{
    return this -> KineticEnergyWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------


  float Event::GetRecoKineticEnergyWithPdg(const int pdg) const{
    return this -> KineticEnergyWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::KineticEnergyWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) return particle_list[i].GetKineticEnergy();
    }
    return 0.;
 }
  //------------------------------------------------------------------------------------------------


  float Event::GetMCModuleMomentumWithPdg(const int pdg) const{
    return this ->ModuleMomentumWithPdg( pdg, m_mc_particles);
  }

  //------------------------------------------------------------------------------------------------


  float Event::GetRecoModuleMomentumWithPdg(const int pdg) const{
    return this -> ModuleMomentumWithPdg( pdg, m_reco_particles);
  }

  //-----------------------------------------------------------------------------------------------

  float Event::ModuleMomentumWithPdg(const int pdg, const ParticleList &particle_list) const{
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(particle_list[i].GetPdgCode() == pdg) return particle_list[i].GetModuleMomentum();
    }
    return 0.;
 }


  //-----------------------------------------------------------------------------------------------

  float Event::GetEnergyLongest_cc1pi( const ParticleList & particle_list) const {
 
    float MaxEnergy = 0.;
    float MaxLength = 0.;
    if( CheckMCTopology( signal_map_cc_1pi ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
	if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
	  MaxEnergy = particle_list[i].GetEnergy();
	  MaxLength = particle_list[i].GetLength();
	}
      }
    }
    return MaxEnergy;
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCEnergyLongest_cc1pi( ) const {
    return this -> GetEnergyLongest_cc1pi( m_mc_particles );
  }


  //-----------------------------------------------------------------------------------------------

  float Event::GetRecoEnergyLongest_cc1pi( ) const {
    return this -> GetEnergyLongest_cc1pi( m_reco_particles );
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetEnergySecondLongest_cc1pi( const ParticleList & particle_list ) const {
 
  float SecondMaxEnergy = 0.;
  float SecondMaxLength = 0;
  float MaxLength = 0.;
  if( CheckMCTopology( signal_map_cc_1pi ) == 1 ){
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
	MaxLength = particle_list[i].GetLength();
      }else if ( MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111  ){
	SecondMaxEnergy = particle_list[i].GetEnergy();
	SecondMaxLength = particle_list[i].GetLength();
      }
    }
    return SecondMaxEnergy;
  }
  return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCEnergySecondLongest_cc1pi( ) const {
    return this -> GetEnergySecondLongest_cc1pi( m_mc_particles );
  }


  //-----------------------------------------------------------------------------------------------

  float Event::GetRecoEnergySecondLongest_cc1pi( ) const {
    return this -> GetEnergySecondLongest_cc1pi( m_reco_particles );
  }


  //-----------------------------------------------------------------------------------------------

  float Event::GetKineticEnergyLongest_cc1pi( const ParticleList & particle_list) const {
 
    float MaxKineticEnergy = 0.;
    float MaxLength = 0.;
    if( CheckMCTopology( signal_map_cc_1pi ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
	if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
	  MaxKineticEnergy = particle_list[i].GetKineticEnergy();
	  MaxLength = particle_list[i].GetLength();
	}
      }
    }
    return MaxKineticEnergy;
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCKineticEnergyLongest_cc1pi( ) const {
    return this -> GetKineticEnergyLongest_cc1pi( m_mc_particles );
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetRecoKineticEnergyLongest_cc1pi( ) const {
    return this -> GetKineticEnergyLongest_cc1pi( m_reco_particles );
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetKineticEnergySecondLongest_cc1pi( const ParticleList & particle_list ) const {
 
  float SecondMaxKineticEnergy = 0.;
  float SecondMaxLength = 0;
  float MaxLength = 0.;
  if( CheckMCTopology( signal_map_cc_1pi ) == 1 ){
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
	MaxLength = particle_list[i].GetLength();
      }else if ( MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111  ){
	SecondMaxKineticEnergy = particle_list[i].GetKineticEnergy();
	SecondMaxLength = particle_list[i].GetLength();
      }
    }
    return SecondMaxKineticEnergy;
  }
  return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCKineticEnergySecondLongest_cc1pi( ) const {
    return this -> GetKineticEnergySecondLongest_cc1pi( m_mc_particles );
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetRecoKineticEnergySecondLongest_cc1pi( ) const {
    return this -> GetKineticEnergySecondLongest_cc1pi( m_reco_particles );
  }
  //--------------------------------------------------------
  float Event::GetModuleMomentumLongest_cc1pi( const ParticleList & particle_list) const {
 
    float MaxMomentum = 0.;
    float MaxLength = 0.;
    if( CheckMCTopology( signal_map_cc_1pi ) == 1  ){
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
	if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
	  MaxMomentum = particle_list[i].GetModuleMomentum();
	  MaxLength = particle_list[i].GetLength();
	}
      }
    }
    return MaxMomentum;
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCModuleMomentumLongest_cc1pi( ) const {
    return this -> GetModuleMomentumLongest_cc1pi( m_mc_particles );
  }


  //-----------------------------------------------------------------------------------------------

  float Event::GetRecoModuleMomentumLongest_cc1pi( ) const {
    return this -> GetModuleMomentumLongest_cc1pi( m_reco_particles );
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetModuleMomentumSecondLongest_cc1pi( const ParticleList & particle_list ) const {
 
  float SecondMaxModuleMomentum = 0.;
  float SecondMaxLength = 0;
  float MaxLength = 0.;
  if( CheckMCTopology( signal_map_cc_1pi ) == 1 ){
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if(MaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111 ){
	MaxLength = particle_list[i].GetLength();
      }else if ( MaxLength > particle_list[i].GetLength() && SecondMaxLength < particle_list[i].GetLength() && particle_list[i].GetPdgCode() != 111  ){
	SecondMaxModuleMomentum = particle_list[i].GetModuleMomentum();
	SecondMaxLength = particle_list[i].GetLength();
      }
    }
    return SecondMaxModuleMomentum;
  }
  return 0.;
  }

  //-----------------------------------------------------------------------------------------------

  float Event::GetMCModuleMomentumSecondLongest_cc1pi( ) const {
    return this -> GetModuleMomentumSecondLongest_cc1pi( m_mc_particles );
  }


  //-----------------------------------------------------------------------------------------------

  float Event::GetRecoModuleMomentumSecondLongest_cc1pi( ) const {
    return this -> GetModuleMomentumSecondLongest_cc1pi( m_reco_particles );
  }

  //----------------------------------------------
  //---------------------------------------------------------------------------------------------

  float Event::GetDeltaEnergy( const ParticleList & particle_list ) const{
    TLorentzVector pion, p2;
      for(unsigned int i = 0; i < particle_list.size(); ++i) {
	if( GetNuanceCode() == 1003 || GetNuanceCode() == 1009 ){
	  if( particle_list[i].GetPdgCode() ==  211 ){
	    pion.SetE( particle_list[i].GetEnergy() );
	    pion.SetVect( particle_list[i].GetMomentum() );
	  }else if( particle_list[i].GetPdgCode() == 2212 ) {
	    p2.SetE( particle_list[i].GetEnergy() );
	    p2.SetVect( particle_list[i].GetMomentum() );
	  }
	  
       	} else if( GetNuanceCode() == 1005 ){
	  if( particle_list[i].GetPdgCode() == 211 ){
	    pion.SetE( particle_list[i].GetEnergy() );
	    pion.SetVect( particle_list[i].GetMomentum() );
	  }else if( particle_list[i].GetPdgCode() == 2112 ) {
	    p2.SetE( particle_list[i].GetEnergy() );
	    p2.SetVect( particle_list[i].GetMomentum() );
	  }
	  
      }
      }
      return (pion + p2).M();
  }
    
  //-----------------------------------------------------------------------------------------------
  float Event::GetMCDeltaEnergy( ) const{
     return this -> GetDeltaEnergy( m_mc_particles );
  }
  //---------------------------------------------------------------------------------------------                                                     

  float Event::GetDeltaEnergy_p( const ParticleList & particle_list ) const{
    TLorentzVector pion, p2;                                                                                          
    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      	if( particle_list[i].GetPdgCode() ==  211 ){
	  pion.SetE( particle_list[i].GetEnergy() );
	  pion.SetVect( particle_list[i].GetMomentum() );
	}else if( particle_list[i].GetPdgCode() == 2212 ) {
	  p2.SetE( particle_list[i].GetEnergy() );
	  p2.SetVect( particle_list[i].GetMomentum() );
	}
    }
    return (pion + p2).M();
  }

  
  //-------------------------------------------------------------------------------------------------

  float Event::GetCC0piNeutrinoEnergy( const ParticleList & particle_list ) const{

    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float V     = 0.02950;   // Nucleon removal energy, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float reco=0, e, p, cth;   // track variables
    
    // Vector of z direction
    TVector3 z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if( particle_list[i].GetPdgCode() == 13 ){
	// Get the values needed
	e    = particle_list[i].GetEnergy();
	p    = particle_list[i].GetMomentum().Mag();
	cth  = (1/p) * (particle_list[i].GetMomentum()).Dot(z);
	
	reco = (1/(m_n - V - e + p*cth))*((m_n - V)*e - (m_mu*m_mu*0.5) + m_n*V - (V*V*0.5) + (m_p*m_p - m_n*m_n)*0.5);
	return reco;
      }
    }
    return reco;
  
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetCC1piNeutrinoEnergy( const ParticleList & particle_list ) const{

    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float m_D   =   1.232;   // Delta mass, GeV
    float V     = 0.02950;   // Nucleon removal energy, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float reco=0, e, p, cth;   // track variables
    // Assuming D+ production ( xsec(n) bigger in Ar ). Modified m_n to consider the D++ production.
    m_n=(22*m_n + 18*m_p)/40.;
    // Vector of z direction
    TVector3 z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      if( particle_list[i].GetPdgCode() == 13 ){
	// Get the values needed
	e    = particle_list[i].GetEnergy();
	p    = particle_list[i].GetMomentum().Mag();
	cth  = (1/p) * (particle_list[i].GetMomentum()).Dot(z);
	
	reco = (1/(m_n - V - e + p*cth))*((m_n - V)*e - (m_mu*m_mu*0.5) + m_n*V - (V*V*0.5) + (m_D*m_D - m_n*m_n)*0.5);
	return reco;
      }
    }
    return reco;
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetCC1piNeutrinoEnergyMethod2( const ParticleList & particle_list ) const{

    // JLAB measured V on Argon - 0.0295 GeV
    // The variables from the branches and get the leaves
    float m_n   = 0.93957;   // Neutron mass, GeV
    float m_p   = 0.93828;   // Neutron mass, GeV
    float m_mu  = 0.10566;   // Muon mass, GeV
    float m_pi  = 0.13957;   // Pion mass, GeV
    float reco=0,p_mu, p_pi, e_mu, e_pi;   // track variables

    // Vector of z direction
    TVector3  z;
    z[0] = 0;
    z[1] = 0;
    z[2] = 1;
    double angle_mu, angle_pi; //cos angle_mu_pi;

    for(unsigned int i = 0; i < particle_list.size(); ++i) {
      for(unsigned int j = 0; j < particle_list.size(); ++j) {
	if( particle_list[i].GetPdgCode() == 13 ){
	  e_mu    = particle_list[i].GetEnergy();
	  p_mu    = particle_list[i].GetMomentum().Mag();
	  angle_mu = particle_list[i].GetCosTheta();
	  if( particle_list[j].GetPdgCode() == 211 ){
	    e_pi    = particle_list[j].GetEnergy();
	    p_pi    = particle_list[j].GetMomentum().Mag();
	    angle_pi =  particle_list[j].GetCosTheta();
	    reco = (m_mu*m_mu+m_pi*m_pi-2*m_p*(e_mu+e_pi)+2*p_mu*p_pi)/(2*(e_mu+e_pi-p_mu*angle_mu-p_pi*angle_pi-m_p));
	    return reco;
	    }
	}
      }
    }
    return reco;
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetMCCC0piNeutrinoEnergy(  ) const{
    return this -> GetCC0piNeutrinoEnergy( m_mc_particles ); 
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetMCCC1piNeutrinoEnergy(  ) const{
    return this -> GetCC1piNeutrinoEnergy( m_mc_particles );   
  }
   //-------------------------------------------------------------------------------------------------
 
  float Event::GetRecoCC0piNeutrinoEnergy(  ) const{
    return this -> GetCC0piNeutrinoEnergy( m_reco_particles ); 
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetRecoCC1piNeutrinoEnergy(  ) const{
    return this -> GetCC1piNeutrinoEnergy( m_reco_particles );   
  }

  //-------------------------------------------------------------------------------------------------

  float Event::GetRecoCC1piNeutrinoEnergyMethod2(  ) const{
    return this -> GetCC1piNeutrinoEnergyMethod2( m_reco_particles );   
  }
  //-------------------------------------------------------------------------------------------------
 
  Particle Event::GetMostEnergeticParticle(const ParticleList &particle_list) const{

    float highest_energy   = -std::numeric_limits<float>::max();
    unsigned int energy_id = std::numeric_limits<unsigned int >::max();

    for(unsigned int i = 0; i < particle_list.size(); ++i){
    
      if(!particle_list[i].GetHasCalorimetry()) continue;

      if(particle_list[i].GetEnergy() > highest_energy) energy_id = i;
    }

    return particle_list[energy_id];

  }
  
  double Event::Efficiency( const std::vector< double > & CountMC, const std::vector< double > & CountTReco, const std::vector< double > & CountReco, const TopologyMap &topology  ) const {
   ofstream rfile ;
   rfile.open( "~/Desktop/Output_Selection_Tool/results.txt" ) ;

   for( int i = 0; i<5; ++i ){
    rfile << "__________________________________________________________"                                                     << "\n";
    rfile                                                                                                                     << "\n";
    rfile << "                 TOPOLOGY NUMBER " << i                                                                         << "\n";
    rfile << "__________________________________________________________"                                                     << "\n";
    rfile << "CountMC = "                       << CountMC[i]                                                                 << "\n";
    rfile << "CountReco = "                     << CountReco[i]                                                               << "\n";
    rfile << "SameCount = "                     << CountTReco[i]                                                              << "\n";
    rfile << "Background = "                    << CountReco[i] - CountTReco[i]                                               << "\n";
    rfile << "Correct Reconstructed Events[%]=" << ( CountTReco[i] / CountMC[i] ) * 100                                       << "\n";
    rfile << "Purity[%]="                       << (( CountTReco[i] ) / CountReco[i] ) * 100                                  << "\n";
    rfile << "Background_Rejection[%]="         << (1-( CountReco[i] - CountTReco[i] ) / ( CountMC[0]+CountMC[1]-CountMC[i] ) ) * 100 << "\n";
    rfile << "__________________________________________________________"                                                     << "\n";
 }

   if( topology == signal_map_NC ) return ( CountTReco[0] / CountMC[0] ) * 100 ;
   if( topology == signal_map_cc_inclusive ) return ( CountTReco[1] / CountMC[1] ) * 100 ;
   if( topology == signal_map_cc_0pi ) return ( CountTReco[2] / CountMC[2] ) * 100 ;
   if( topology == signal_map_cc_1pi ) return ( CountTReco[3] / CountMC[3] ) * 100 ;
   if( topology == signal_map_cc_pi0 ) return ( CountTReco[4] / CountMC[4] ) * 100 ;
   return 0;
  }

  void Event::SaveTopologyMatrix( const ParticleMatrix & Count_MC_Topology, const ParticleMatrix & Count_TReco_Topology, const ParticleMatrix & Count_Reco_Topology ) const {
   ofstream TMfile ;
   TMfile.open( "~/Desktop/Output_Selection_Tool/TopologyMatrix.txt" ) ;
    TMfile                                                                                                                       << "\n";
    TMfile << "____________________________________________________________"                                                     << "\n";
    TMfile                                                                                                                       << "\n";
    TMfile << "  TOPOLOGY MATRIX - TRUE RECO  (#_TReco / #_Total_MC) : "                                                         << "\n";
    TMfile << "____________________________________________________________"                                                     << "\n";
    for( unsigned int i = 0 ; i < 5; ++i ){
      TMfile << "(";
      for( unsigned int k = 0 ; k < 5 ; ++k ) {
	if( Count_TReco_Topology[i][k]!=0 ){
	  TMfile << ( Count_TReco_Topology[i][k] / Count_MC_Topology[i][k] ) * 100 << "   ,   ";
	} else TMfile << "   --   ";
    }
      TMfile <<")"                                                                                                               << "\n";
    }
    
  }

  void Event::EventInformationParticles( std::string event_file, const int event_number ) const {
    ofstream efile ;
    efile.open( event_file , std::ofstream::app);

    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << " EVENT NUMBER =                                            " << event_number                 << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "TRUE EVENTS      : "                                                                         << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "   muons         : " << CountMCParticlesWithPdg(13)                                          << "\n";
    efile << "   pi+/-         : " << CountMCParticlesWithPdg(211) + CountMCParticlesWithPdg(-211)         << "\n";
    efile << "   pi0           : " << CountMCParticlesWithPdg(111)                                         << "\n";
    efile << "   protons       : " << CountMCParticlesWithPdg(2212)                                        << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "TRUE TOPOLOGY    : "                                                                         << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    if(CheckMCTopology(signal_map_NC))           efile << "   NC            : " << "TRUE"                  << "\n";
    if(CheckMCTopology(signal_map_cc_inclusive)) efile << "   ccincl.       : " << "TRUE"                  << "\n";
    if(CheckMCTopology(signal_map_cc_0pi))       efile << "   cc0pi         : " << "TRUE"                  << "\n";
    if(CheckMCTopology(signal_map_cc_1pi))       efile << "   cc1pi+/-      : " << "TRUE"                  << "\n";
    if(CheckMCTopology(signal_map_cc_pi0))       efile << "   cc1pi0        : " << "TRUE"                  << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << " SELECTED EVENTS :                                         "                                 << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "   muons         : " << CountRecoParticlesWithPdg(13)                                        << "\n";
    efile << "   pi+/-         : " << CountRecoParticlesWithPdg(211) + CountRecoParticlesWithPdg(-211)     << "\n";
    efile << "   pi0           : " << CountRecoParticlesWithPdg(111)                                       << "\n";
    efile << "   protons       : " << CountRecoParticlesWithPdg(2212)                                      << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    efile << "SELECTED TOPOLOGY: "                                                                         << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";
    if(CheckRecoTopology(signal_map_NC))           efile << "   NC            : " << "TRUE"                << "\n";
    if(CheckRecoTopology(signal_map_cc_inclusive)) efile << "   ccincl.       : " << "TRUE"                << "\n";
    if(CheckRecoTopology(signal_map_cc_0pi))       efile << "   cc0pi         : " << "TRUE"                << "\n";
    if(CheckRecoTopology(signal_map_cc_1pi))       efile << "   cc1pi+/-      : " << "TRUE"                << "\n";
    if(CheckRecoTopology(signal_map_cc_pi0))       efile << "   cc1pi0        : " << "TRUE"                << "\n";
    efile << "-----------------------------------------------------------"                                 << "\n";


  }

  void Event::EventProperties( const TopologyMap &topology, std::string event_file, const int event_number ) const {
    ofstream lfile, afile, Kfile ;
    lfile.open( event_file+="_length.txt" , std::ofstream::app);
    afile.open( event_file+="_angle.txt"  , std::ofstream::app);
    Kfile.open( event_file+="_KineticEnergy.txt" , std::ofstream::app);

    if( CheckMCTopology( topology ) ) { /** Change the topology here **/
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "EVENT NUMBER      : " << event_number                                                        << "\n";
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "LENGTH INFORMATION: "                                                                        << "\n";
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "TRUE EVENTS       : "                                                                        << "\n";
       lfile << "-----------------------------------------------------------"                                 << "\n";
       lfile << "Muon length         : " << GetMCLengthWithPdg(  13  )         << "\n";
       lfile << "pi+/- length        : " << GetMCLengthWithPdg( 211  )  << "\n";
       lfile << "pi0   length        : " << GetMCLengthWithPdg( 111  )  << "\n";
       lfile << "p length            : " << GetMCLengthWithPdg( 2212 ) << "\n";
       
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "EVENT NUMBER       : " << event_number                                                       << "\n";
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "ANGULAR INFORMATION: "                                                                       << "\n";
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "TRUE EVENTS        : "                                                                       << "\n"; 
       afile << "-----------------------------------------------------------"                                 << "\n";
       afile << "Muon angle       : " << GetMCCosThetaWithPdg(  13  )  << "\n";
       afile << "Pion angle       : " << GetMCCosThetaWithPdg( 211  )  << "\n";
       afile << "pi0   angle      : " << GetMCCosThetaWithPdg( 111  )  << "\n";
       afile << "Proton angle     : " << GetMCCosThetaWithPdg( 2212 )  << "\n";

       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "EVENT NUMBER       : " << event_number                                                       << "\n";
       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "KINETIC ENERGY [GeV] : "                                                                     << "\n";
       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "TRUE EVENTS        : "                                                                       << "\n"; 
       Kfile << "-----------------------------------------------------------"                                 << "\n";
       Kfile << "Muon Kenergy       : " << GetMCKineticEnergyWithPdg(  13  )  << "\n";
       Kfile << "Pion Kenergy       : " << GetMCKineticEnergyWithPdg( 211  )  << "\n";
       
       if( CountMCParticlesWithPdg(2212)!=0 ) Kfile << "Proton Kenergy     : " << GetMCKineticEnergyWithPdg( 2212 )  << "\n";

       if( CheckRecoTopology( topology ) ) {

	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 lfile << "SIGNAL EVENTS     : "                                                                        << "\n";
	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckMCTopology(signal_map_NC) )lfile << "Muon length          : " << GetRecoLengthWithPdg(  13  )  << "\n";
	 if( CheckMCTopology(signal_map_cc_1pi) )lfile << "pi+/- length      : "  << GetRecoLengthWithPdg( 211  ) << "\n";
	 if( CheckMCTopology(signal_map_cc_pi0) )lfile << "pi0   length      : "  << GetRecoLengthWithPdg( 111  ) << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 )lfile << "p length            : " << GetRecoLengthWithPdg( 2212 )  << "\n";
	 

	 afile << "-----------------------------------------------------------"                                 << "\n";
	 afile << "SIGNAL EVENTS      : "                                                                       << "\n"; 
	 afile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckMCTopology(signal_map_NC)) afile << "Muon angle              : " << GetRecoCosThetaWithPdg(  13  )  << "\n";
	 if( CheckMCTopology(signal_map_cc_1pi) )afile << "Pion angle         : " << GetRecoCosThetaWithPdg( 211  )  << "\n";
	 if( CheckMCTopology(signal_map_cc_pi0) )afile << "pi0   angle        : " << GetRecoCosThetaWithPdg( 111  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) afile << "Proton angle        : " << GetRecoCosThetaWithPdg( 2212 )  << "\n";
 	 
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "KINETIC ENERGY [GeV] : "                                                                     << "\n";
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "TRUE EVENTS        : "                                                                       << "\n"; 
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckMCTopology(signal_map_NC) ) Kfile << "Muon Kenergy         : " << GetRecoKineticEnergyWithPdg(  13  )  << "\n";
	 if( CheckMCTopology(signal_map_cc_1pi) )Kfile << "Pion Kenergy       : " << GetRecoKineticEnergyWithPdg( 211  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) Kfile << "Proton Kenergy      : " << GetRecoKineticEnergyWithPdg( 2212 )  << "\n";
       }

     }
    if( CheckRecoTopology( topology ) ) { // Change topology here

	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 lfile << "SELECTED EVENTS   : "                                                                        << "\n";
	 lfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckRecoTopology(signal_map_NC) ) lfile << "Muon length       : " << GetRecoLengthWithPdg(  13  )         << "\n";
	 if( CheckRecoTopology(signal_map_cc_1pi) )lfile << "pi+/- length      : " << GetRecoLengthWithPdg( 211  )  << "\n";
	 if( CheckRecoTopology(signal_map_cc_pi0) )lfile << "pi0   length      : " << GetRecoLengthWithPdg( 111  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) lfile << "p length          : " << GetRecoLengthWithPdg( 2212 ) << "\n";
	 
	 afile << "-----------------------------------------------------------"                                 << "\n";
	 afile << "SELECTED EVENTS      : "                                                                     << "\n"; 
	 afile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckRecoTopology(signal_map_NC) ) afile << "Muon angle         : " << GetRecoCosThetaWithPdg(  13  )  << "\n";
	 if( CheckRecoTopology(signal_map_cc_1pi) )afile << "Pion angle         : " << GetRecoCosThetaWithPdg( 211  )  << "\n";
	 if( CheckRecoTopology(signal_map_cc_pi0) )afile << "pi0   angle        : " << GetRecoCosThetaWithPdg( 111  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) afile << "Proton angle       : " << GetRecoCosThetaWithPdg( 2212 )  << "\n";

	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "KINETIC ENERGY [GeV] : "                                                                     << "\n";
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 Kfile << "TRUE EVENTS        : "                                                                       << "\n"; 
	 Kfile << "-----------------------------------------------------------"                                 << "\n";
	 if( !CheckRecoTopology(signal_map_NC) ) Kfile << "Muon Kenergy       : " << GetRecoKineticEnergyWithPdg(  13  )  << "\n";
	 if( CheckRecoTopology(signal_map_cc_1pi) )Kfile << "Pion Kenergy       : " << GetRecoKineticEnergyWithPdg( 211  )  << "\n";
	 if( CountMCParticlesWithPdg(2212)!=0 ) Kfile << "Proton Kenergy     : " << GetRecoKineticEnergyWithPdg( 2212 )  << "\n";
     }

  }

} // selection
