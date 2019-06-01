
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
  Particle z, t, tb, wp, b, wm, bm, q, qb;
  TLorentzVector Pz, Pt, Ptb, Pwp, Pb, Pwm, Pbm, Pq, Pqb, Pl;
  vector<double> Px_hadrons, Py_hadrons, Pz_hadrons, E_hadrons;
 
 
  int zmi=2;
 // for(zmi=0; zmi<3; zmi++){ 
 
   TTree *tree = new TTree("partons","partons");  
  TTree *tree1 = new TTree("hadrons","hadrons");  
 
  TBranch *b_z   = tree->Branch("parton[0]", &Pz);
  TBranch *b_t   = tree->Branch("parton[1]", &Pt); 
  TBranch *b_tb  = tree->Branch("parton[2]", &Ptb);
  TBranch *b_wp  = tree->Branch("parton[3]", &Pwp);
  TBranch *b_b   = tree->Branch("parton[4]", &Pb);
  TBranch *b_wm  = tree->Branch("parton[5]", &Pwm);
  TBranch *b_bm  = tree->Branch("parton[6]", &Pbm);
  TBranch *b_q   = tree->Branch("parton[7]", &Pq);
  TBranch *b_qb  = tree->Branch("parton[8]", &Pqb);
  TBranch *b_l  = tree->Branch("parton[9]", &Pl);
  
  
  TBranch *b_hadrons_px = tree1->Branch("px", &Px_hadrons);
  TBranch *b_hadrons_py = tree1->Branch("py", &Py_hadrons);
  TBranch *b_hadrons_pz = tree1->Branch("pz", &Pz_hadrons);
  TBranch *b_hadrons_e  = tree1->Branch("e", &E_hadrons);

    
    Pythia pythia;    
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Init:showChangedParticleData=0");      

    sprintf(filename,"Beams:LHEF=/home/ehep/Downloads/products/MG5_aMC_v2_6_5/pp2zp2tt2qqblv_%d/Events/run_01/unweighted_events.lhe", z_mass[zmi]);  
    pythia.readString("Beams:frameType = 4");
    pythia.readString(filename);
    pythia.init();

    int nEvent=10, nFill=0, nToFill=0; 
 
    nEvent=1000;
 
    bool lepton20;  
    for (int iEvent = 0;iEvent<nEvent; ++iEvent) {  
      lepton20 = false;
      
      if (!pythia.next())  { cout<<"Go away, no time to say hi."<<endl; continue;}    

      for(int i=0; i<pythia.event.size(); i++){
        if(pythia.event[i].isFinal() && pythia.event[i].isLepton()){
          Pl.SetPxPyPzE(pythia.event[i].p().px(), pythia.event[i].p().py(), pythia.event[i].p().pz(), pythia.event[i].p().e() );
          b_l->Fill();
          if(Pl.Pt() > 20 ) { lepton20 = true; }
        }
      }

      if(lepton20){
        nToFill++;
        for(int i=0; i< pythia.event.size(); i++){

        if( pythia.event[i].id() == Zp_id ){

          z  = pythia.event[i];
          t  = FinalDaughter(&pythia, z, 0);
          tb = FinalDaughter(&pythia, z, 1);
          wp = FinalDaughter(&pythia, t, 0);
          b  = FinalDaughter(&pythia, t, 1);
          wm = FinalDaughter(&pythia, tb, 0);
          bm = FinalDaughter(&pythia, tb, 1);
          q  = FinalDaughter(&pythia, wp, 0);
          qb = FinalDaughter(&pythia, wp, 1);

          
          if( t.id() == top_id && tb.id() == -top_id && wp.id() == w_id && wm.id() == -w_id && b.id() == b_id && bm.id() == -b_id){
            nFill++;
            Pz.SetPxPyPzE( z.p().px(), z.p().py(), z.p().pz(), z.p().e()   );
            Pt.SetPxPyPzE( t.p().px(), t.p().py(), t.p().pz(), t.p().e()   );
            Ptb.SetPxPyPzE( tb.p().px(), tb.p().py(), tb.p().pz(), tb.p().e()   );
            Pwp.SetPxPyPzE( wp.p().px(), wp.p().py(), wp.p().pz(), wp.p().e()   );
            Pb.SetPxPyPzE( b.p().px(), b.p().py(), b.p().pz(), b.p().e()   );
            Pwm.SetPxPyPzE( wm.p().px(), wm.p().py(), wm.p().pz(), wm.p().e()   );
            Pbm.SetPxPyPzE( bm.p().px(), bm.p().py(), bm.p().pz(), bm.p().e()   );
            Pq.SetPxPyPzE( q.p().px(), q.p().py(), q.p().pz(), q.p().e()   );
            Pqb.SetPxPyPzE( qb.p().px(), qb.p().py(), qb.p().pz(), qb.p().e()   );
          
            FillHadrons(&pythia, &Px_hadrons, &Py_hadrons, &Pz_hadrons, &E_hadrons);    
    
            b_z->Fill();
            b_t->Fill();
            b_tb->Fill();
            b_wp->Fill();
            b_b->Fill();
            b_wm->Fill();
            b_bm->Fill();
            b_q->Fill();
            b_qb->Fill();
            
            b_hadrons_px->Fill();
            b_hadrons_py->Fill();
            b_hadrons_pz->Fill();
            b_hadrons_e->Fill();

            Px_hadrons.clear();
            Py_hadrons.clear();
            Pz_hadrons.clear();
            E_hadrons.clear();
          } //all particle checked
            break;
          } //Z dikha    
        }  //eventsize
      }  //lepton20  
      cout<<"Events Loading......"<<100.0*float(iEvent+1)/float(nEvent)<<" % "<<"\r";
 
    }  //event loop 
    cout<<endl;
    
    
    sprintf(filename, "./pp2zp2tt2qqblv_%d.root", z_mass[zmi]);

    cout<<"number of events : "<<nEvent<<endl;  
    cout<<"number of event with atleast one lepton with pT>20 GeV : "<<nToFill<<endl;
    cout<<"number of events to write on disk : "<<nFill<<endl;
    cout<<"creating : " <<filename<<endl;
    cout<<"filling : " <<nFill<<" events"<<endl<<endl;
    TFile f(filename,"RECREATE");  
    tree->Fill();
    tree1->Fill();
    tree->Write();
    tree1->Write();    
    f.Close();
      
   	return 0;
}




/*
     if( pythia.event[i].id() ==21 && !pythia.event[i].isQuark() && pythia.event[i].daughterList().size() >3 ){        

        if(match_found(gluon_decayvector, pythia.event[i].daughterList()) == 0 ){ 
          nGluon++; 
          for(int k=pythia.event[i].daughter1()  ; k  <=pythia.event[i].daughter2() ; k++)
           gluon_decayvector.push_back(k);
                        
           cout<<nGluon<<"  gluon jetsize : "<<setw(20)<<pythia.event[pythia.event[i].mother1()].name()<<setw(20)<<pythia.event[i].mother2()<<endl;            
        }
      }
      
      if( pythia.event[i].isQuark() && pythia.event[i].daughterList().size() > 3 ){        

        if(  match_found(quark_decayvector, pythia.event[i].daughterList())==0  ){ 
          nQuark++; 
          //cout<<"new quark "<<pythia.event[i].name()<<setw(20)<<pythia.event[pythia.event[pythia.event[i].daughter1()].mother2()].name()<<endl;
          for(int k=pythia.event[i].daughter1()  ; k  <=pythia.event[i].daughter2() ; k++)
           quark_decayvector.push_back(k);             
          cout<<nQuark<<"  quark jetsize : "<<setw(20)<<pythia.event[pythia.event[i].daughter1()].name()<<setw(20)<<pythia.event[i].daughter2()<<endl; 
        }
      }
*/




































