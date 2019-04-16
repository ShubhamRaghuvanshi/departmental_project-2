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


int main(){

  Pythia pythia;
  pythia.readFile("qcdmultijet.cmnd");
  int nEvent = pythia.mode("Main:numberOfEvents");
  pythia.init();


  TTree *tree = new TTree("qcdmultijet","qcdmultijet");

  Particle hadron;
  Particle photon;
  vector<double> Px_hadrons, Py_hadrons, Pz_hadrons, E_hadrons;

  TBranch *v_hadrons_px = tree->Branch("hadrons_px", &Px_hadrons);
  TBranch *v_hadrons_py = tree->Branch("hadrons_py", &Py_hadrons);
  TBranch *v_hadrons_pz = tree->Branch("hadrons_pz", &Pz_hadrons);
  TBranch *v_hadrons_e  = tree->Branch("hadrons_e", &E_hadrons);


 
  int n_fill=0;  


  cout<<"##############  "<<nEvent<<endl;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {continue;}
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
    n_fill++; 
    Px_hadrons.clear();
    Py_hadrons.clear();
    Pz_hadrons.clear();
    E_hadrons.clear();
  }   // loop over events

  cout<<"filled events : "<<n_fill<<endl;
  TFile f("qcd_rawhadrons.root","RECREATE");  
  tree->Fill();
  tree->Write();
  f.Close();

  pythia.stat();


return 0;
}


/*
    for (int i = 0; i < pythia.event.size(); ++i){
      
      if(pythia.event[i].isFinal()){
        if(  abs(pythia.event[i].id())==11 || abs(pythia.event[i].id())==13  ){
          if(pythia.event[i].p().pT() > 20.0 && pythia.event[i].p().eta()<2.4 ){
}}}




   /*           if(hadron.isFinal() && hadron.id()==22 && pythia.event[hadron.mother1()].id()==hadroHn.id() ) {
              cout<<setw(20)<<"radating  "<<setw(20)<<hadron.name()<<setw(20)<<pythia.event[hadron.mother1()].name()<<setw(20)<<endl; 
              photon = pythia.event[j];
                while(pythia.event[photon.mother1()].id()==photon.id())                
                  photon = pythia.event[photon.mother1()];  
              cout<<setw(20)<<"radating  "<<setw(20)<<photon.name()<<setw(20)<<pythia.event[photon.mother1()].name()<<setw(20)<<endl;   
            }
    */



//    cout<<setw(20)<<i<<setw(20)<<pythia.event[i].name()<<setw(20)<<pythia.event[i].id()<<setw(20)<<pythia.event[i].p().pT()<<setw(20)<<pythia.event[i].p().eta()<<endl;

















