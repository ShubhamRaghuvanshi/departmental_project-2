
//This is updated version of main130.cc event code 
// It uses TLorentz vector to store 4 momenta and pushes it into root tree branch
// on event by event basis in order to avoid making large array of quantities as was 
// done in main130.cc 
#include <iostream>
#include <vector>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TLorentzVector.h" 
#include "fstream" 

using namespace Pythia8;
using namespace std;


const int top_id=6;
const int Zp_id = 32;
const int e_id = 11;
const int mu_id = 13;
const int tau_id = 15;

const int Wplus_id = 24;
const int b_id = 5;
const int m_w = 80;
const int m_top =172;




int main(){

  Pythia pythia;
  pythia.readFile("newgaugeboson.cmnd");
  int nEvent = pythia.mode("Main:numberOfEvents");
  pythia.init();

  // these are 4 momenta of particles in our hard scatter of interest
  // hadron are 4 momenta of all the hadrons produced in the hard scatter of interest
  // hadrons20 are stored only when there is one lepton with pT>20Gev in an event, it is used for jet reconstruction.  

  TLorentzVector P_z,P_top,P_topbar,P_wplus,P_wminus,P_b,P_bbar,P_lepton,P_neutrino,P_quark,P_antiquark;
  TLorentzVector P_hadron;
  vector<double> Px_hadrons, Py_hadrons, Pz_hadrons, E_hadrons;
  vector<double> Px_hadrons20, Py_hadrons20, Pz_hadrons20, E_hadrons20;
  int id_wh;

//  ofstream mfile;  
//  mfile.open ("example.txt");

  TTree *tree = new TTree("z2ttbar_1l_4j","z2ttbar_1l_4j");  

  TBranch *v_z  = tree->Branch("Z", &P_z);
  TBranch *v_top = tree->Branch("top", &P_top);
  TBranch *v_topbar = tree->Branch("topbar", &P_topbar);
  TBranch *v_wplus = tree->Branch("wplus", &P_wplus);
  TBranch *v_wminus = tree->Branch("wminus",&P_wminus);
  TBranch *v_b = tree->Branch("b", &P_b);  
  TBranch *v_bbar = tree->Branch("bbar",&P_bbar);
  TBranch *v_lepton = tree->Branch("lepton",&P_lepton);
  TBranch *v_neutrino = tree->Branch("neutrino",&P_neutrino);
  TBranch *v_quark = tree->Branch("quark",&P_quark);
  TBranch *v_antiquark = tree->Branch("antiquark",&P_antiquark);

  TBranch *v_hadrons_px = tree->Branch("hadrons_px", &Px_hadrons);
  TBranch *v_hadrons_py = tree->Branch("hadrons_py", &Py_hadrons);
  TBranch *v_hadrons_pz = tree->Branch("hadrons_pz", &Pz_hadrons);
  TBranch *v_hadrons_e  = tree->Branch("hadrons_e", &E_hadrons);
  
  TBranch *id_Wh  = tree->Branch("id_wh", &id_wh);

  // the decaying particles  
  Particle particle;
  Particle top,topbar;
  Particle wplus,wminus;
  // final state hadron and lepton  
  Particle hadron,lepton;  
  
  int z2ttbar_events=0;
  int lepton_events=0, not_lepton_events=0, lepton20_events=0, n_hadrons=0;
  int failed_events=0;
  bool event_lepton,event_lepton20;

  

  int w0=0,w1=0,w2=0;

  cout<<"##############  "<<nEvent<<endl;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    id_wh=0;
    if (!pythia.next()) {failed_events++; continue;}
     
 
    for (int i = 0; i < pythia.event.size(); ++i){ 

      particle = pythia.event[i];
      // identify a Z' which decays to ttbar
  
      //check the if
      if(particle.id()==Zp_id && pythia.event[particle.daughter1()].id() == top_id && pythia.event[particle.daughter2()].id() == -top_id )  { 
        
        z2ttbar_events++;

        top     = pythia.event[ particle.daughter1() ];
        topbar  = pythia.event[ particle.daughter2() ];
        
        // find top at it's final decay stage            
        int m=15;
        for(int i=0;i<m;i++){  
          if(top.id() == pythia.event[top.daughter1()].id() )             
            top = pythia.event[top.daughter1()];          
          else break;
        }
        // print if it didnt decay to W+ and b 

        //check this if
        if(pythia.event[top.daughter1()].id() != Wplus_id || pythia.event[top.daughter2()].id() != b_id  ){
          cout<<"t didnt decay to W+ and b : "<<iEvent<<endl;
          cout<<"t decays to : "<<pythia.event[top.daughter1()].name()<<setw(20)<<pythia.event[top.daughter2()].name()<<endl;
          cout<<"ignoring event : "<<iEvent<<endl;
          failed_events++;
          continue;
        }
        
        // find top at it's final decay stage               
        for(int i=0;i<m;i++){  
          if(topbar.id() == pythia.event[topbar.daughter1()].id() )             
            topbar = pythia.event[topbar.daughter1()];          
          else break;
          }
        // print if it didnt decay to W- and bbar 
        if(pythia.event[topbar.daughter1()].id() != -Wplus_id || pythia.event[topbar.daughter2()].id() != -b_id  ){

          cout<<"tb didnt decay to W+ and b : "<<iEvent<<endl;
          cout<<"tb decays to : "<<pythia.event[topbar.daughter1()].name()<<setw(20)<<pythia.event[topbar.daughter2()].name()<<endl;
          cout<<"ignoring event : "<<iEvent<<endl;
          failed_events++;
          continue;
        }

        wplus =   pythia.event[top.daughter1()];
        wminus = pythia.event[topbar.daughter1()];  
        

        // find W+ at it's final decay stage
        for(int i=0;i<m;i++){                       
          if(wplus.id() == pythia.event[wplus.daughter1()].id() ) 
            wplus = pythia.event[wplus.daughter1()];
          else break; 
        }

        // find W- at it's final decay stage
        for(int i=0;i<m;i++){                       
          if(wminus.id() == pythia.event[wminus.daughter1()].id() ) 
            wminus = pythia.event[wminus.daughter1()];
          else break; 
        }
         

        bool quark1  = pythia.event[wplus.daughter1()].isQuark() && pythia.event[wplus.daughter2()].isQuark(); 
        bool lepton1=  pythia.event[wplus.daughter1()].isLepton() && pythia.event[wplus.daughter2()].isLepton();

        bool quark2  = pythia.event[wminus.daughter1()].isQuark() && pythia.event[wminus.daughter2()].isQuark() ;
        bool lepton2=  pythia.event[wminus.daughter1()].isLepton() && pythia.event[wminus.daughter2()].isLepton() ;    
 
        event_lepton = (quark1 && lepton2) || (quark2 && lepton1) ;

        // FIll data for the event
        if(event_lepton){

          
    
          if(quark2 && lepton1 && pythia.event[wplus.daughter1()].p().pT()>20 && abs( pythia.event[wplus.daughter1()].p().eta())<2.4 ) {

            id_wh = -1;  w2++; 

//cout<<iEvent<<setw(20)<<topbar.p().px()<<setw(20)<<topbar.p().py()<<setw(20)<<topbar.p().pz()<<setw(20)<<topbar.p().e()<<setw(20)<<id_wh<<endl; 

            P_z.SetPxPyPzE      ( particle.p().px() ,particle.p().py(), particle.p().pz(), particle.p().e()  ) ;
            P_top.SetPxPyPzE    ( top.p().px() ,top.p().py(), top.p().pz(), top.p().e()  ) ;
            P_topbar.SetPxPyPzE ( topbar.p().px() ,topbar.p().py(), topbar.p().pz(), topbar.p().e()  ) ;


            P_wplus.SetPxPyPzE  (wplus.p().px(), wplus.p().py(), wplus.p().pz(), wplus.p().e()  ) ;      
            P_wminus.SetPxPyPzE (wminus.p().px(), wminus.p().py(), wminus.p().pz(), wminus.p().e()  ) ;      
            P_b.SetPxPyPzE(pythia.event[top.daughter2()].p().px(),pythia.event[top.daughter2()].p().py(),pythia.event[top.daughter2()].p().pz(),pythia.event[top.daughter2()].p().e()) ; 
        P_bbar.SetPxPyPzE(pythia.event[topbar.daughter2()].p().px(),pythia.event[topbar.daughter2()].p().py(),pythia.event[topbar.daughter2()].p().pz(),pythia.event[topbar.daughter2()].p().e()) ; 


            P_lepton.SetPxPyPzE(pythia.event[wplus.daughter1()].p().px(), pythia.event[wplus.daughter1()].p().py(), pythia.event[wplus.daughter1()].p().pz(), pythia.event[wplus.daughter1()].p().e()); 
            P_neutrino.SetPxPyPzE(pythia.event[wplus.daughter2()].p().px(), pythia.event[wplus.daughter2()].p().py(), pythia.event[wplus.daughter2()].p().pz(), pythia.event[wplus.daughter2()].p().e()) ; 

            P_quark.SetPxPyPzE(pythia.event[wminus.daughter1()].p().px(), pythia.event[wminus.daughter1()].p().py(), pythia.event[wminus.daughter1()].p().pz(), pythia.event[wminus.daughter1()].p().e()); 
            P_antiquark.SetPxPyPzE(pythia.event[wminus.daughter2()].p().px() ,pythia.event[wminus.daughter2()].p().py(), pythia.event[wminus.daughter2()].p().pz(), pythia.event[wminus.daughter2()].p().e()); 

 
          lepton_events++;    
          v_z->Fill();
          v_top->Fill();  
          v_topbar->Fill();
          v_wplus->Fill();  
          v_wminus->Fill();
          v_b->Fill();
          v_bbar->Fill();
          v_lepton->Fill();
          v_neutrino->Fill();
          v_quark->Fill();
          v_antiquark->Fill();
          id_Wh->Fill();

        //Fill hadron data 
          for (int j = 0; j < pythia.event.size(); ++j){ 
            hadron = pythia.event[j]; 


           if(  (hadron.isHadron() || (hadron.id()==22 && pythia.event[hadron.mother1()].id()==111  ) )  && hadron.isFinal()  ){

              Px_hadrons.push_back( hadron.p().px() );  
              Py_hadrons.push_back( hadron.p().py() );
              Pz_hadrons.push_back( hadron.p().pz() );
              E_hadrons.push_back( hadron.p().e() );
            }
          }
          v_hadrons_px->Fill();
          v_hadrons_py->Fill();
          v_hadrons_pz->Fill();
          v_hadrons_e->Fill();
          
          Px_hadrons.clear();
          Py_hadrons.clear();
          Pz_hadrons.clear();
          E_hadrons.clear();
  

          }
             
           
          if(quark1 && lepton2 && pythia.event[wminus.daughter1()].p().pT()>20 && abs( pythia.event[wminus.daughter1()].p().eta())<2.4  ){


            
            id_wh = 1;  w1++; 
      

            P_z.SetPxPyPzE      ( particle.p().px() ,particle.p().py(), particle.p().pz(), particle.p().e()  ) ;
            P_top.SetPxPyPzE    ( top.p().px() ,top.p().py(), top.p().pz(), top.p().e()  ) ;
            P_topbar.SetPxPyPzE ( topbar.p().px() ,topbar.p().py(), topbar.p().pz(), topbar.p().e()  ) ;


            P_wplus.SetPxPyPzE  (wplus.p().px(), wplus.p().py(), wplus.p().pz(), wplus.p().e()  ) ;      
            P_wminus.SetPxPyPzE (wminus.p().px(), wminus.p().py(), wminus.p().pz(), wminus.p().e()  ) ;      
            P_b.SetPxPyPzE(pythia.event[top.daughter2()].p().px(),pythia.event[top.daughter2()].p().py(),pythia.event[top.daughter2()].p().pz(),pythia.event[top.daughter2()].p().e()) ; 
        P_bbar.SetPxPyPzE(pythia.event[topbar.daughter2()].p().px(),pythia.event[topbar.daughter2()].p().py(),pythia.event[topbar.daughter2()].p().pz(),pythia.event[topbar.daughter2()].p().e()) ; 




            P_lepton.SetPxPyPzE(pythia.event[wminus.daughter1()].p().px(), pythia.event[wminus.daughter1()].p().py(), pythia.event[wminus.daughter1()].p().pz(), pythia.event[wminus.daughter1()].p().e()); 
            P_neutrino.SetPxPyPzE(pythia.event[wminus.daughter2()].p().px(), pythia.event[wminus.daughter2()].p().py(), pythia.event[wminus.daughter2()].p().pz(), pythia.event[wminus.daughter2()].p().e()) ; 

            P_quark.SetPxPyPzE(pythia.event[wplus.daughter1()].p().px(), pythia.event[wplus.daughter1()].p().py(), pythia.event[wplus.daughter1()].p().pz(), pythia.event[wplus.daughter1()].p().e()); 
            P_antiquark.SetPxPyPzE(pythia.event[wplus.daughter2()].p().px() ,pythia.event[wplus.daughter2()].p().py(), pythia.event[wplus.daughter2()].p().pz(), pythia.event[wplus.daughter2()].p().e()); 

          lepton_events++;    
          v_z->Fill();
          v_top->Fill();  
          v_topbar->Fill();
          v_wplus->Fill();  
          v_wminus->Fill();
          v_b->Fill();
          v_bbar->Fill();
          v_lepton->Fill();
          v_neutrino->Fill();
          v_quark->Fill();
          v_antiquark->Fill();
          id_Wh->Fill();

         //Fill hadron data 
          for (int j = 0; j < pythia.event.size(); ++j){ 
            hadron = pythia.event[j]; 
              if(  (hadron.isHadron() || (hadron.id()==22 && pythia.event[hadron.mother1()].id()==111  ) )  && hadron.isFinal()  ){
              Px_hadrons.push_back( hadron.p().px() );  
              Py_hadrons.push_back( hadron.p().py() );
              Pz_hadrons.push_back( hadron.p().pz() );
              E_hadrons.push_back( hadron.p().e() );
            }
          }
          v_hadrons_px->Fill();
          v_hadrons_py->Fill();
          v_hadrons_pz->Fill();
          v_hadrons_e->Fill();
          
          Px_hadrons.clear();
          Py_hadrons.clear();
          Pz_hadrons.clear();
          E_hadrons.clear();
 
          }
          }
              
        
        }  //z2tt case
  } //eventsize        break;

    if(id_wh ==0) w0++;
      }  //loop over events
  
 
  ofstream myfile;
  myfile.open ("Zp2ttbar_1l_Event2000.txt");
 
  cout<<setw(40)<<"Total events : "<<setw(30)<<nEvent<<endl;
  cout<<setw(40)<<"ww -> l1q4 events : "<<setw(30)<<lepton_events<<endl;
  cout<<setw(40)<<"expected ww->1l3h events : "<<setw(30)<<2.0*nEvent*(3.0*0.108660)*0.675581<<endl;
  cout<<setw(40)<<"other events : "<<setw(30)<<not_lepton_events<<endl;
  cout<<setw(40)<<"ww -> l1q3 events + other events : "<<setw(30)<<lepton_events + not_lepton_events<<endl;
  cout<<setw(40)<<"Z' -> ttbar events : "<<setw(30)<<z2ttbar_events<<endl; 
  cout<<setw(40)<<"l20 events : "<<setw(30)<<lepton20_events<<endl;
  cout<<setw(40)<<"Hadron fill events : "<<setw(30)<<lepton_events<<endl;
  cout<<setw(40)<<"Total hadrons : "<<setw(30)<<n_hadrons<<endl;

  myfile<<setw(40)<<"Total events : "<<setw(30)<<nEvent<<endl;
  myfile<<setw(40)<<"ww -> l1q4 events : "<<setw(30)<<lepton_events<<endl;
  myfile<<setw(40)<<"expected ww->1l3h events : "<<setw(30)<<2.0*nEvent*(3.0*0.108660)*0.675581<<endl;
  myfile<<setw(40)<<"other events : "<<setw(30)<<not_lepton_events<<endl;
  myfile<<setw(40)<<"ww -> l1q3 events + other events : "<<setw(30)<<lepton_events + not_lepton_events<<endl;
  myfile<<setw(40)<<"Z' -> ttbar events : "<<setw(30)<<z2ttbar_events<<endl; 
  myfile<<setw(40)<<"l20 events : "<<setw(30)<<lepton20_events<<endl;
  myfile<<setw(40)<<"Hadron fill events : "<<setw(30)<<lepton_events<<endl;
  myfile<<setw(40)<<"Total hadrons : "<<setw(30)<<n_hadrons<<endl;

  myfile.close();

  cout<<setw(20)<<"w :"<<w0<<setw(20)<<w1<<setw(20)<<w2<<setw(20)<<w0+w1+w2<<endl;

  //write tree  to file

  TFile f("Zp2ttbar_1l_Event2000.root","RECREATE");  
  tree->Fill();
  tree->Write();
  f.Close();

  pythia.stat();
  return 0; 
}































