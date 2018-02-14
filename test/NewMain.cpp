#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <sstream>
#include <numeric>
#include <time.h>
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TObjArray.h"

using namespace selection;

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start: Local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Running all files " << std::endl;
  
  // Counters
  unsigned int correctly_reconstructed_cc  = 0;
  unsigned int true_topology_cc            = 0;
  unsigned int reco_topology_cc            = 0;

  unsigned int correctly_reconstructed_nc  = 0;
  unsigned int true_topology_nc            = 0;
  unsigned int reco_topology_nc            = 0;

  unsigned int correctly_reconstructed_0pi = 0;
  unsigned int true_topology_0pi           = 0;
  unsigned int reco_topology_0pi           = 0;
  unsigned int true_ccqe                   = 0;
  unsigned int true_ccmec                  = 0;
  unsigned int true_ccdis                  = 0;
  unsigned int true_ccres                  = 0;
  unsigned int true_cccoh                  = 0;
  unsigned int reco_ccqe                   = 0;
  unsigned int reco_ccmec                  = 0;
  unsigned int reco_ccdis                  = 0;
  unsigned int reco_ccres                  = 0;
  unsigned int reco_cccoh                  = 0;
  unsigned int good_ccqe                   = 0;
  unsigned int good_ccmec                  = 0;
  unsigned int good_ccdis                  = 0;
  unsigned int good_ccres                  = 0;
  unsigned int good_cccoh                  = 0;

  unsigned int correctly_reconstructed_1pi = 0;
  unsigned int true_topology_1pi           = 0;
  unsigned int reco_topology_1pi           = 0;

  unsigned int correctly_reconstructed_pi0 = 0;
  unsigned int true_topology_pi0           = 0;
  unsigned int reco_topology_pi0           = 0;

  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap cc_signal_map;
  TopologyMap nc_signal_map;
  TopologyMap cc0pi_signal_map;
  TopologyMap cc1pi_signal_map;
  TopologyMap ccpi0_signal_map;
 
  std::vector< int > mu;
  std::vector< int > pi;
  std::vector< int > pi0;
  
  mu.push_back( 13 );
  pi.push_back( 211 );
  pi.push_back(-211 );
  pi0.push_back( 111 );
  
  cc_signal_map.insert( std::make_pair( mu,  1 ) );
  nc_signal_map.insert( std::make_pair( mu,  0 ) );
  
  cc0pi_signal_map.insert( std::make_pair( mu,  1 ) );
  cc0pi_signal_map.insert( std::make_pair( pi,  0 ) );
  cc0pi_signal_map.insert( std::make_pair( pi0, 0 ) );

  cc1pi_signal_map.insert( std::make_pair( mu,  1 ) );
  cc1pi_signal_map.insert( std::make_pair( pi,  1 ) );
  
  ccpi0_signal_map.insert( std::make_pair( mu,  1 ) );
  ccpi0_signal_map.insert( std::make_pair( pi0, 1 ) );
  
  for( unsigned int i = 0; i < 398; ++i ){
  
    // Get the filename for each 2D histogram
    std::stringstream ss;
    ss.clear();
    
    std::string name;
    name.clear();
    
    char file_name[1024];
    
    ss << "/hepstore/rjones/Samples/FNAL/analysis_trees/all/3486578_" << i <<"/output_file.root";
    name = ss.str();
            
    strcpy( file_name, name.c_str() );
      
    EventSelectionTool::LoadEventList(file_name, events);
  }

  // Neutrino energy histograms
  TH1F *h_reco_energy      = new TH1F("h_reco_energy",      "#nu_{#mu} CC 0#pi neutrino energy",100,-0.5,2);
  TH1F *h_good_reco_energy = new TH1F("h_good_reco_energy", "#nu_{#mu} CC 0#pi neutrino energy",100,-0.5,2);
  TH1F *h_true_energy      = new TH1F("h_true_energy",      "#nu_{#mu} CC 0#pi neutrino energy",100,-0.5,2);
  TH2F *h_energy_cos       = new TH2F("h_energy_cos",       "#nu_{#mu} CC 0#pi neutino energy and opening angle", 50, -0.5, 2, 50, -1, 1);

  // Get CC 0pi true and reconstructed neutrino energies
  std::vector<float> true_neutrino_energy, reco_neutrino_energy, good_reco_neutrino_energy, good_reco_cos_theta;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    if(e.CheckMCTopology(cc_signal_map) && e.CheckRecoTopology(cc_signal_map)) correctly_reconstructed_cc++;
    if(e.CheckRecoTopology(cc_signal_map)) reco_topology_cc++;
    if(e.CheckMCTopology(cc_signal_map))   true_topology_cc++;

    if(e.CheckMCTopology(nc_signal_map) && e.CheckRecoTopology(nc_signal_map)) correctly_reconstructed_nc++;
    if(e.CheckRecoTopology(nc_signal_map)) reco_topology_nc++;
    if(e.CheckMCTopology(nc_signal_map))   true_topology_nc++;
    
    if(e.CheckMCTopology(cc1pi_signal_map) && e.CheckRecoTopology(cc1pi_signal_map)) correctly_reconstructed_1pi++;
    if(e.CheckRecoTopology(cc1pi_signal_map)) reco_topology_1pi++;
    if(e.CheckMCTopology(cc1pi_signal_map))   true_topology_1pi++;
    
    if(e.CheckMCTopology(ccpi0_signal_map) && e.CheckRecoTopology(ccpi0_signal_map)) correctly_reconstructed_pi0++;
    if(e.CheckRecoTopology(ccpi0_signal_map)) reco_topology_pi0++;
    if(e.CheckMCTopology(ccpi0_signal_map))   true_topology_pi0++;
    
    if(e.CheckMCTopology(cc0pi_signal_map) && e.CheckRecoTopology(cc0pi_signal_map)) {
      
      // Counter
      correctly_reconstructed_0pi++;
      
      // Get the well selected reconstructed neutrino energy
      ParticleList parts       = e.GetRecoParticleList();
      unsigned int n_particles = e.GetRecoParticleList().size();
      
      // Get reconstructed energy
      for( unsigned int i = 0; i < n_particles; ++i ) if(parts[i].GetPdgCode() == 13) {
        good_reco_neutrino_energy.push_back(e.GetCC0piRecoNeutrinoEnergy(parts[i]));

        TVector3 z;
        z[0] = 0;
        z[1] = 0;
        z[2] = 1;

        float p    = parts[i].GetMomentum().Mag();
        float cth  = (1/p) * (parts[i].GetMomentum()).Dot(z);
        good_reco_cos_theta.push_back(cth);
      }
      
      // Nuance codes
      if(e.GetIsCC()){
      
        if(e.GetNuanceCode() == 1001) good_ccqe++;
        if(e.GetNuanceCode() == 10)   good_ccmec++;
        if(e.GetNuanceCode() == 2 || e.GetNuanceCode() == 1091) good_ccdis++;
        if(e.GetNuanceCode() == 1003 || e.GetNuanceCode() == 1004 || e.GetNuanceCode() == 1005) good_ccres++;
        if(e.GetNuanceCode() == 1097) good_cccoh++;
          
      }
    }
    if(e.CheckRecoTopology(cc0pi_signal_map)) {
      
      // Counter
      reco_topology_0pi++;

      // Get the well selected reconstructed neutrino energy
      ParticleList parts       = e.GetRecoParticleList();
      unsigned int n_particles = e.GetRecoParticleList().size();
      
      // Get reconstructed energy
      for( unsigned int i = 0; i < n_particles; ++i ) if(parts[i].GetPdgCode() == 13) {
        reco_neutrino_energy.push_back(e.GetCC0piRecoNeutrinoEnergy(parts[i]));
      }
      
      // Nuance codes
      if(e.GetIsCC()){
      
        if(e.GetNuanceCode() == 1001) reco_ccqe++;
        if(e.GetNuanceCode() == 10)   reco_ccmec++;
        if(e.GetNuanceCode() == 2 || e.GetNuanceCode() == 1091) reco_ccdis++;
        if(e.GetNuanceCode() == 1003 || e.GetNuanceCode() == 1004 || e.GetNuanceCode() == 1005) reco_ccres++;
        if(e.GetNuanceCode() == 1097) reco_cccoh++;
          
      }
    }
    if(e.CheckMCTopology(cc0pi_signal_map)) {
     
      // Counter
      true_topology_0pi++;
    
      // Get true neutrino energy
      true_neutrino_energy.push_back(e.GetTrueNuEnergy());
      
      // Nuance codes
      if(e.GetIsCC()){
      
        if(e.GetNuanceCode() == 1001) true_ccqe++;
        if(e.GetNuanceCode() == 10)   true_ccmec++;
        if(e.GetNuanceCode() == 2 || e.GetNuanceCode() == 1091) true_ccdis++;
        if(e.GetNuanceCode() == 1003 || e.GetNuanceCode() == 1004 || e.GetNuanceCode() == 1005) true_ccres++;
        if(e.GetNuanceCode() == 1097) true_cccoh++;
          
      }
    }
  }
  for( unsigned int i = 0; i < reco_neutrino_energy.size(); ++i)      h_reco_energy->Fill(reco_neutrino_energy[i]);
  for( unsigned int i = 0; i < true_neutrino_energy.size(); ++i)      h_true_energy->Fill(true_neutrino_energy[i]);
  for( unsigned int i = 0; i < good_reco_neutrino_energy.size(); ++i) {
    h_good_reco_energy->Fill(good_reco_neutrino_energy[i]);
    h_energy_cos->Fill(good_reco_neutrino_energy[i],good_reco_cos_theta[i]);
  }

  std::cout << "===========================================================" << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed NC events     : " << correctly_reconstructed_nc/double(true_topology_nc) * 100                         << std::endl;
  std::cout << " Impurity of reconstructed NC events                 : " << (reco_topology_nc - correctly_reconstructed_nc)/double(reco_topology_nc) * 100    << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC events     : " << correctly_reconstructed_cc/double(true_topology_cc) * 100                         << std::endl;
  std::cout << " Impurity of reconstructed CC events                 : " << (reco_topology_cc - correctly_reconstructed_cc)/double(reco_topology_cc) * 100    << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC 1Pi events : " << correctly_reconstructed_1pi/double(true_topology_1pi) * 100                       << std::endl;
  std::cout << " Impurity of reconstructed CC 1Pi events             : " << (reco_topology_1pi - correctly_reconstructed_1pi)/double(reco_topology_1pi) * 100 << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC Pi0 events : " << correctly_reconstructed_pi0/double(true_topology_pi0) * 100                       << std::endl;
  std::cout << " Impurity of reconstructed CC Pi0 events             : " << (reco_topology_pi0 - correctly_reconstructed_pi0)/double(reco_topology_pi0) * 100 << std::endl; 
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Percentage of correctly reconstructed CC 0Pi events : " << correctly_reconstructed_0pi/double(true_topology_0pi) * 100                       << std::endl;
  std::cout << " Impurity of reconstructed CC 0Pi events             : " << (reco_topology_0pi - correctly_reconstructed_0pi)/double(reco_topology_0pi) * 100 << std::endl;
  std::cout << " Percentage of true " << std::endl;
  std::cout << "      CC QE                                          : " << true_ccqe  / double(true_topology_0pi) * 100 << std::endl;
  std::cout << "      CC MEC                                         : " << true_ccmec / double(true_topology_0pi) * 100 << std::endl;
  std::cout << "      CC RES                                         : " << true_ccres / double(true_topology_0pi) * 100 << std::endl;
  std::cout << "      CC DIS                                         : " << true_ccdis / double(true_topology_0pi) * 100 << std::endl;
  std::cout << "      CC COH                                         : " << true_cccoh / double(true_topology_0pi) * 100 << std::endl;
  std::cout << " Percentage of reconstructed " << std::endl;
  std::cout << "      CC QE                                          : " << reco_ccqe  / double(reco_topology_0pi) * 100 << std::endl;
  std::cout << "      CC MEC                                         : " << reco_ccmec / double(reco_topology_0pi) * 100 << std::endl;
  std::cout << "      CC RES                                         : " << reco_ccres / double(reco_topology_0pi) * 100 << std::endl;
  std::cout << "      CC DIS                                         : " << reco_ccdis / double(reco_topology_0pi) * 100 << std::endl;
  std::cout << "      CC COH                                         : " << reco_cccoh / double(reco_topology_0pi) * 100 << std::endl;
  std::cout << " Percentage of well reconstructed " << std::endl;
  std::cout << "      CC QE                                          : " << good_ccqe  / double(correctly_reconstructed_0pi) * 100 << std::endl;
  std::cout << "      CC MEC                                         : " << good_ccmec / double(correctly_reconstructed_0pi) * 100 << std::endl;
  std::cout << "      CC RES                                         : " << good_ccres / double(correctly_reconstructed_0pi) * 100 << std::endl;
  std::cout << "      CC DIS                                         : " << good_ccdis / double(correctly_reconstructed_0pi) * 100 << std::endl;
  std::cout << "      CC COH                                         : " << good_cccoh / double(correctly_reconstructed_0pi) * 100 << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "===========================================================" << std::endl;

  TCanvas *c = new TCanvas();
  TLegend *l = new TLegend( 0.58, 0.68, 0.88, 0.88 );

  l->AddEntry( h_true_energy,      " True ",               "l" );
  l->AddEntry( h_reco_energy,      " Reconstructed ",      "l" );
  l->AddEntry( h_good_reco_energy, " Well reconstructed ", "l" );
    
  h_true_energy->SetLineColor(2);
  h_true_energy->SetStats(kFALSE);
  h_reco_energy->SetLineColor(4);
  h_reco_energy->SetStats(kFALSE);
  h_reco_energy->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  h_good_reco_energy->SetLineColor(6);
  h_good_reco_energy->SetStats(kFALSE);

  h_reco_energy->Draw();
  h_good_reco_energy->Draw("same");
  h_true_energy->Draw("same");
  l->Draw();

  c->SaveAs("plots/cc0pi_nu_energy.root");
  c->Clear();

  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);

  h_energy_cos->SetStats(kFALSE);
  h_energy_cos->GetXaxis()->SetTitle("Neutrino energy [GeV]");
  h_energy_cos->GetYaxis()->SetTitle("cos#theta_{#mu}");
  h_energy_cos->Draw("colz");
 
  c->SaveAs("plots/cc0pi_energy_cos.root");
  c->Clear();
  
  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
