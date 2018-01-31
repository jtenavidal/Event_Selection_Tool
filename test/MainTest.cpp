#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>

using namespace selection;

int MainTest(){

  std::string filename = "/hepstore/rjones/Samples/LArSoft_reconstruction_temp/30_tree_test.root";

  EventSelectionTool::EventList events;
  EventSelectionTool::LoadEventList(filename, events);

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    TopologyMap signal_map_all;
   
    std::vector< int > cc_0pi_mu;
    std::vector< int > cc_0pi_pi;
    
    cc_0pi_mu.push_back( 13 );
    cc_0pi_pi.push_back( 211 );
    cc_0pi_pi.push_back(-211 );
    cc_0pi_pi.push_back( 111 );
    
    signal_map_all.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map_all.insert( std::make_pair( cc_0pi_pi, 0 ) );

    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << " RECO : " << std::endl;
    std::cout << "   muons   : " << e.CountRecoParticlesWithPdg(13) << std::endl;
    std::cout << "   pi+/-   : " << e.CountRecoParticlesWithPdg(211) + e.CountRecoParticlesWithPdg(-211) << std::endl;
    std::cout << "   pi0     : " << e.CountRecoParticlesWithPdg(111) << std::endl;
    std::cout << "   protons : " << e.CountRecoParticlesWithPdg(2212) << std::endl;
    std::cout << "   cc0pi   : " << e.CheckRecoTopology(signal_map_all) << std::endl;

    std::cout << " MC   : " << std::endl;
    std::cout << "   muons   : " << e.CountMCParticlesWithPdg(13) << std::endl;
    std::cout << "   pi+/-   : " << e.CountMCParticlesWithPdg(211) + e.CountMCParticlesWithPdg(-211) << std::endl;
    std::cout << "   pi0     : " << e.CountMCParticlesWithPdg(111) << std::endl;
    std::cout << "   protons : " << e.CountMCParticlesWithPdg(2212) << std::endl;
    std::cout << "   cc0pi   : " << e.CheckMCTopology(signal_map_all) << std::endl;
    
  }
  
  return 0;

} // MainTest()
