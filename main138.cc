
//QCD background

#include <iostream>
#include <vector>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TLorentzVector.h" 
#include "fstream" 
#include "Pythia8/Pythia.h"

#include "pythiaFunc.cc"
using namespace Pythia8;
using namespace std;

int main() {

  char filename[200];
  vector<double> Px_hadrons, Py_hadrons, Pz_hadrons, E_hadrons;
   
  TTree *tree1 = new TTree("hadrons","hadrons");  
  
  TBranch *b_hadrons_px = tree1->Branch("px", &Px_hadrons);
  TBranch *b_hadrons_py = tree1->Branch("py", &Py_hadrons);
  TBranch *b_hadrons_pz = tree1->Branch("pz", &Pz_hadrons);
  TBranch *b_hadrons_e  = tree1->Branch("e", &E_hadrons);
    
  Pythia pythia;    
  
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 200.");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Init:showChangedParticleData=0");      

  pythia.init();

  int nEvent=10;
 
  nEvent=10000;
 
  for (int iEvent = 0;iEvent<nEvent; ++iEvent) {
  
    if (!pythia.next())  { cout<<"Go away, no time to say hi."<<endl; continue;}   

    FillHadrons(&pythia, &Px_hadrons, &Py_hadrons, &Pz_hadrons, &E_hadrons);

    b_hadrons_px->Fill();
    b_hadrons_py->Fill();
    b_hadrons_pz->Fill();
    b_hadrons_e->Fill();

    Px_hadrons.clear();
    Py_hadrons.clear();
    Pz_hadrons.clear();
    E_hadrons.clear();

    cout<<"Events Loading......"<<100.0*float(iEvent+1)/float(nEvent)<<" % "<<"\r";
 
   }  //event loop 
   cout<<endl;
    
    
  sprintf(filename, "./HEPTopTagger2/QCD20.root");

  cout<<"number of events : "<<nEvent<<endl;  

  cout<<"creating : " <<filename<<endl;
  cout<<"filling : " <<b_hadrons_px->GetEntries()<<" events"<<endl<<endl;
  TFile f(filename,"RECREATE");  
  tree1->Fill();
  tree1->Write();    
  f.Close();
    
 	return 0;
}



































