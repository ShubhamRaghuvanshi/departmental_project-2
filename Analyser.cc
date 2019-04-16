#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>  
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TLegendEntry.h"
#include "TLegend.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "particleproperty.cc"
#include "./HEPTopTagger.hh"
#include "Func.cc"

using namespace std;

using namespace fastjet;

float top_mass = 172.0;
float w_mass = 80;

class Analyser{

  private:
    //outputs
     //inputs
    string filename,treename;
    bool all;     

  public:
  
    int n1,n2;  
 
    ParticleProperty Z,Top, Topbar, B, Bbar, Wplus, Wminus, Quark, AntiQuark;
    ParticleProperty Hadronic_Top, Hadronic_B, Hadronic_W;   
  //  ParticleProperty Final_State_Hardons;
    ParticleProperty Top_Chi, Top_Kin, Top_Hep; 
    ParticleProperty B_Chi, B_Kin, B_Hep;
    ParticleProperty W_Chi, W_Kin, W_Hep; 
    
    ParticleProperty W_Chireco, W_Kinreco ,W_Hepreco;     
    ParticleProperty Top_Chireco, Top_Kinreco, Top_Hepreco;
    ParticleProperty B_Chireco, B_Kinreco, B_Hepreco;       

    ParticleProperty W_Chireco_tagged, W_Kinreco_tagged, W_Hepreco_tagged;     
    ParticleProperty Top_Chireco_tagged, Top_Kinreco_tagged, Top_Hepreco_tagged;
    ParticleProperty B_Chireco_tagged, B_Kinreco_tagged, B_Hepreco_tagged;       

    int n_tagged_kin, n_tagged_chi, n_tagged_hep;
  
    Analyser(){}

    Analyser(string file_name, string tree_name,int iEvent, int fEvent, bool All ){
         filename = file_name;
         treename  = tree_name;     
         n1 = iEvent;   
         n2 = fEvent;
         all = All;   
         n_tagged_chi=0;
         n_tagged_kin=0;
         n_tagged_hep=0;  
      } 

// -1 for error in reading
    int RunAnalyser(float R, int mass_index){
      TTree* tree;
      TFile f(filename.c_str(), "READ");    
      f.GetObject(treename.c_str(), tree);        

      TBranch  *bvpx=0,*bvpy=0,*bvpz=0,*bve=0;
      TBranch *br_z=0, *br_top=0, *br_topbar=0, *br_b=0, *br_bbar=0, *br_idwh=0, *br_wplus=0, *br_wminus=0, *br_quark=0, *br_antiquark=0;

      vector<double>  *hardon_px=0,*hardon_py=0,*hardon_pz=0,*hardon_e=0;
      TLorentzVector *P_z=0, *P_top=0 ,*P_topbar=0, *P_b=0, *P_bbar=0, *P_wplus=0, *P_wminus=0, *P_quark=0, *P_antiquark=0;
      TLorentzVector hadronic_top, hadronic_b, hadronic_w;   


      int id_wh ;
      int n_top_to_be_reconstructed=0;

      tree->SetBranchAddress("hadrons_px",&hardon_px,&bvpx);
      tree->SetBranchAddress("hadrons_py",&hardon_py,&bvpy);
      tree->SetBranchAddress("hadrons_pz",&hardon_pz,&bvpz);
      tree->SetBranchAddress("hadrons_e",&hardon_e ,&bve);   
    
      tree->SetBranchAddress("Z",&P_z,&br_z);
      tree->SetBranchAddress("top",&P_top,&br_top);
      tree->SetBranchAddress("topbar",&P_topbar,&br_topbar);
      tree->SetBranchAddress("wplus",&P_wplus,&br_wplus);
      tree->SetBranchAddress("wminus",&P_wminus,&br_wminus);
      tree->SetBranchAddress("b",&P_b,&br_b);
      tree->SetBranchAddress("bbar",&P_bbar,&br_bbar);
      tree->SetBranchAddress("quark",&P_quark,&br_quark);
      tree->SetBranchAddress("antiquark",&P_antiquark,&br_antiquark);
      tree->SetBranchAddress("id_wh",&id_wh,&br_idwh);

      Long64_t tentry  = tree->LoadTree(0);
      bvpx->GetEntry(tentry);       bvpy->GetEntry(tentry);       bvpz->GetEntry(tentry);   bve->GetEntry(tentry);
      br_top->GetEntry(tentry);     br_topbar->GetEntry(tentry);     
      br_wplus->GetEntry(tentry);   br_wminus->GetEntry(tentry);  br_quark->GetEntry(tentry);   br_antiquark->GetEntry(tentry);
      br_b->GetEntry(tentry);       br_bbar->GetEntry(tentry);
      br_idwh->GetEntry(tentry);    br_z->GetEntry(tentry);  

 


      cout<<bvpx->GetEntries()-1<<setw(20)<<bvpy->GetEntries()-1<<setw(20)<<bvpz->GetEntries()-1<<setw(20)<<bve->GetEntries()-1<<endl<<endl;  

      vector<PseudoJet> particles; 
      vector<PseudoJet> jets;
      vector<PseudoJet> hep_jets; 
      TLorentzVector top_chi, w_chi, b_chi;
      TLorentzVector top_kin, w_kin, b_kin;         

     // float R = 1.0;  
      int z_mass[] = {500,1000,2000};
      char pngname[200];
      sprintf(pngname, "njets(p_T^{jet}>20), mz=%d",z_mass[mass_index]);  
      TH1F *njets20 = new TH1F(pngname, pngname, 10 , 0,10);
      sprintf(pngname, "njets(p_T^{jet}>150), mz=%d",z_mass[mass_index]);  
      TH1F *njets150 = new TH1F(pngname, pngname, 10 , 0,10);

      if(all==1) {n1=0;    n2=br_top->GetEntries()-2;  }  
      for(int iEvent=n1;iEvent<=n2;iEvent++){

        bvpx->GetEntry(iEvent);     bvpy->GetEntry(iEvent);        bvpz->GetEntry(iEvent);   bve->GetEntry(iEvent);
        br_top->GetEntry(iEvent);   br_topbar->GetEntry(iEvent); 
        br_wplus->GetEntry(iEvent); br_wminus->GetEntry(iEvent);  br_quark->GetEntry(iEvent);   br_antiquark->GetEntry(iEvent);
        br_b->GetEntry(iEvent);     br_bbar->GetEntry(iEvent);  
        br_idwh->GetEntry(iEvent);  br_z->GetEntry(iEvent);


        Z.push_back_momenta( *P_z);
        Top.push_back_momenta(*P_top);
        Topbar.push_back_momenta(*P_topbar);
        B.push_back_momenta(*P_b);
        Bbar.push_back_momenta(*P_bbar);
        Wplus.push_back_momenta(*P_wplus);
        Wminus.push_back_momenta(*P_wminus);
        Quark.push_back_momenta(*P_quark);
        AntiQuark.push_back_momenta(*P_antiquark);
        
        //get 4 momenta of w,b,t for hadronically decay

        if(id_wh == 1)       { hadronic_top = *P_top;    hadronic_b = *P_b;    hadronic_w = *P_wplus;   n_top_to_be_reconstructed++; }
        else if(id_wh == -1) { hadronic_top = *P_topbar; hadronic_b = *P_bbar; hadronic_w = *P_wminus;  n_top_to_be_reconstructed++; }  
        else return -1;

        Hadronic_B.push_back_momenta(hadronic_b);
        Hadronic_W.push_back_momenta(hadronic_w);
        Hadronic_Top.push_back_momenta(hadronic_top);


        // ######################################## Jet reconstruction ################################### 
        for(int i = 0;i<hardon_px->size();i++){  
          particles.push_back( PseudoJet(hardon_px->at(i),hardon_py->at(i),hardon_pz->at(i),hardon_e->at(i)) );
       //   Final_State_Hardons.push_back_momenta(hardon_px->at(i),hardon_py->at(i),hardon_pz->at(i),hardon_e->at(i));
        }

        // place R loop here

//        JetDefinition jet_def(cambridge_algorithm, R);
         JetDefinition jet_def(antikt_algorithm, R);
       ClusterSequence cs(particles, jet_def);
        jets = sorted_by_pt(cs.inclusive_jets(20));           
        hep_jets = sorted_by_pt(cs.inclusive_jets(150));
  
        njets20->Fill(jets.size() );
        njets150->Fill(hep_jets.size());

//cout<<setw(20)<<jets.size()<<setw(20)<<hep_jets.size()<<endl;
       if(jets.size()>=3){
       // cout<<setw(20)<<topjetreco_chi(jets, top_chi, w_chi, b_chi) <<setw(20)<<topjetreco_kin(jets, top_kin, w_kin, b_kin)<<endl;
  
        if( topjetreco_chi(jets, top_chi, w_chi, b_chi)==1 ) {
          B_Chi.push_back_momenta( b_chi );
          W_Chi.push_back_momenta( w_chi );
          Top_Chi.push_back_momenta( top_chi );

          if( abs( top_chi.M()- top_mass)<20.0 && abs( w_chi.M()- w_mass)<20.0 ){

            n_tagged_chi++;      

            Top_Chireco.push_back_momenta(hadronic_top);
            W_Chireco.push_back_momenta(hadronic_w);
            B_Chireco.push_back_momenta(hadronic_b);        
 
            Top_Chireco_tagged.push_back_momenta( top_chi );
            W_Chireco_tagged.push_back_momenta(w_chi);
            B_Chireco_tagged.push_back_momenta(b_chi);        
         }
       }//tagged

        if( topjetreco_kin(jets, top_kin, w_kin, b_kin)==1 ){
          B_Kin.push_back_momenta( b_kin );
          W_Kin.push_back_momenta( w_kin );
          Top_Kin.push_back_momenta( top_kin );

          if( abs( top_kin.M()- top_mass)<20.0 && abs( w_kin.M()- w_mass)<20.0 ){

            n_tagged_kin++;

            Top_Kinreco.push_back_momenta(hadronic_top);
            W_Kinreco.push_back_momenta(hadronic_w);
            B_Kinreco.push_back_momenta(hadronic_b);  

            Top_Kinreco_tagged.push_back_momenta( top_kin );
            W_Kinreco_tagged.push_back_momenta(w_kin);
            B_Kinreco_tagged.push_back_momenta(b_kin);      
          }
        }//tagged
     }  //njets.>3
     

     for(unsigned ijet=0; ijet<hep_jets.size(); ijet++){      

          HEPTopTagger::HEPTopTagger tagger(hep_jets[ijet]);

          // Unclustering, Filtering & Subjet Settings
          tagger.set_max_subjet_mass(30.);
          tagger.set_mass_drop_threshold(0.8);
          tagger.set_filtering_R(0.3);
          tagger.set_filtering_n(5);
          tagger.set_filtering_minpt_subjet(30.); 

          // How to select among candidates
          tagger.set_mode(HEPTopTagger::TWO_STEP_FILTER);
  
          // Requirements to accept a candidate
          tagger.set_top_minpt(150); 
          tagger.set_top_mass_range(120, 240.); 
          tagger.set_fw(0.15); 

          // Run the tagger
          tagger.run();      
          // Look at output ie have a tag:
          if (tagger.is_tagged()) {
            Top_Hep.push_back_momenta(tagger.t().px(), tagger.t().py(), tagger.t().pz(), tagger.t().e());
            W_Hep.push_back_momenta(tagger.W().px(), tagger.W().py(), tagger.W().pz(), tagger.W().e());
            B_Hep.push_back_momenta(tagger.b().px(), tagger.b().py(), tagger.b().pz(), tagger.b().e());
        
            if( abs(tagger.t().m() - top_mass)<20.0 && abs(tagger.W().m() - w_mass)<20.0 ){
        
              n_tagged_hep++  ;

              Top_Hepreco.push_back_momenta(hadronic_top.Px(), hadronic_top.Py(), hadronic_top.Pz(), hadronic_top.E());
              W_Hepreco.push_back_momenta(hadronic_w.Px(), hadronic_w.Py(), hadronic_w.Pz(), hadronic_w.E());
              B_Hepreco.push_back_momenta(hadronic_b.Px(), hadronic_w.Py(), hadronic_b.Pz(), hadronic_b.E());

              Top_Hepreco_tagged.push_back_momenta(tagger.t().px(), tagger.t().py(), tagger.t().pz(), tagger.t().e());
              W_Hepreco_tagged.push_back_momenta(tagger.W().px(), tagger.W().py(), tagger.W().pz(), tagger.W().e());
              B_Hepreco_tagged.push_back_momenta(tagger.b().px(), tagger.b().py(), tagger.b().pz(), tagger.b().e());

            }
         }//tagged

      }

      particles.clear();  
      jets.clear();
      hep_jets.clear(); 
      cout<<"Processing data......"<<100.0*float(iEvent+1)/float(n2-n1+1)<<" % "<<"\r";


    }  //loop over all events
    cout<<endl;

    Z.setnames("Z","");                 Z.make_properties();                      
    Top.setnames("top","");             Top.make_properties();                      
    Topbar.setnames("topbar","");       Topbar.make_properties();                 
    B.setnames("b","");                 B.make_properties();                      
    Bbar.setnames("bbar","");           Bbar.make_properties();                   
    Wplus.setnames("Wplus","");         Wplus.make_properties();                  
    Wminus.setnames("minus","");        Wminus.make_properties();                 
    Quark.setnames("q1","");             Quark.make_properties();                  
    AntiQuark.setnames("q2","");     AntiQuark.make_properties();              

    Top_Chireco.setnames("RecoTopChi","");    Top_Chireco.make_properties();            
    W_Chireco.setnames("RecoWChi","");        W_Chireco.make_properties();               
    B_Chireco.setnames("RecoBChi","");        B_Chireco.make_properties();              

    Top_Chireco_tagged.setnames("TopChi_tagged","");      Top_Chireco_tagged.make_properties();            
    W_Chireco_tagged.setnames("WChi_tagged","");          W_Chireco_tagged.make_properties();               
    B_Chireco_tagged.setnames("BChi_tagged","");          B_Chireco_tagged.make_properties();              


    Top_Kinreco.setnames("RecoTopKin","");    Top_Kinreco.make_properties();            
    W_Kinreco.setnames("RecoWKin","");        W_Kinreco.make_properties();               
    B_Kinreco.setnames("RecoBKin","");        B_Kinreco.make_properties();              

    Top_Kinreco_tagged.setnames("TopKin_tagged","");      Top_Kinreco_tagged.make_properties();            
    W_Kinreco_tagged.setnames("WKin_tagged","");          W_Kinreco_tagged.make_properties();               
    B_Kinreco_tagged.setnames("BKin_tagged","");          B_Kinreco_tagged.make_properties();              
//
    Top_Hepreco.setnames("RecoTopHep","");                Top_Hepreco.make_properties();            
    W_Hepreco.setnames("RecoWHep","");                    W_Hepreco.make_properties();               
    B_Hepreco.setnames("RecoBHep","");                    B_Hepreco.make_properties();              

    Top_Hepreco_tagged.setnames("TopHep_tagged","");      Top_Hepreco_tagged.make_properties();            
    W_Hepreco_tagged.setnames("WHep_tagged","");          W_Hepreco_tagged.make_properties();               
    B_Hepreco_tagged.setnames("Hep_tagged","");           B_Hepreco_tagged.make_properties();              


//


    Hadronic_Top.setnames("Hadronic_Top",""); Hadronic_Top.make_properties();           
    Hadronic_B.setnames("Hadronic_b","");     Hadronic_B.make_properties();             
    Hadronic_W.setnames("Hadronic_w","");     Hadronic_W.make_properties();             

 //   Final_State_Hardons.setnames("Hadrons",""); Final_State_Hardons.make_properties();    
  
    B_Chi.setnames("b_Chi","");                 B_Chi.make_properties();                    
    W_Chi.setnames("w_Chi","");                 W_Chi.make_properties();                    
    Top_Chi.setnames("top_Chi","");             Top_Chi.make_properties();                  

    B_Kin.setnames("b_Kin","");               B_Kin.make_properties();                  
    W_Kin.setnames("w_Kin","");               W_Kin.make_properties();                   
    Top_Kin.setnames("top_Kin","");           Top_Kin.make_properties();                

    B_Hep.setnames("b_Hep","");               B_Hep.make_properties();                   
    W_Hep.setnames("w_Hep","");               W_Hep.make_properties();                   
    Top_Hep.setnames("top_Hep","");           Top_Hep.make_properties();                

//cout<<setw(20)<<Top.GetPartName()<<setw(20)<<endl;
//cout<<setw(20)<<Top.GetPropName(0)<<setw(20)<<endl;

    sprintf(pngname, "njets%d", mass_index );
    TCanvas *c_njets = new TCanvas(pngname,pngname,1800,1000);
    c_njets->Divide(2,1);   
    c_njets->cd(1);

    njets20->SetLineWidth(2);    
    njets20->Draw("hist");
    c_njets->cd(2);  
    njets150->SetLineWidth(2);
    njets150->Draw("hist");
 
    int event_directory = mkdir("./Z2ttbar_EventAnalysis", S_IRUSR | S_IWUSR | S_IXUSR);
    sprintf(pngname, "./Z2ttbar_EventAnalysis/AR_%f",R);
    int cone_directory = mkdir(pngname, S_IRUSR | S_IWUSR | S_IXUSR);

    sprintf(pngname , "./Z2ttbar_EventAnalysis/AR_%f/njets%d.png",R, mass_index);
    c_njets->SaveAs(pngname);  

return 0;

    } //run analyser 

  //jet reconstruction only
   int RunJetReco(float R){

    TTree* tree;
    TFile f(filename.c_str(), "READ");    
    f.GetObject(treename.c_str(), tree);        

    TBranch  *bvpx=0,*bvpy=0,*bvpz=0,*bve=0;
    vector<double>  *hardon_px=0,*hardon_py=0,*hardon_pz=0,*hardon_e=0;
 
    tree->SetBranchAddress("hadrons_px",&hardon_px,&bvpx);
    tree->SetBranchAddress("hadrons_py",&hardon_py,&bvpy);
    tree->SetBranchAddress("hadrons_pz",&hardon_pz,&bvpz);
    tree->SetBranchAddress("hadrons_e",&hardon_e ,&bve);   
    
    Long64_t tentry  = tree->LoadTree(0);
    bvpx->GetEntry(tentry);       bvpy->GetEntry(tentry);       bvpz->GetEntry(tentry);   bve->GetEntry(tentry);

    cout<<bvpx->GetEntries()-1<<setw(20)<<bvpy->GetEntries()-1<<setw(20)<<bvpz->GetEntries()-1<<setw(20)<<bve->GetEntries()-1<<endl<<endl;  

    vector<PseudoJet> particles; 
    vector<PseudoJet> jets;
    vector<PseudoJet> hep_jets; 
    TLorentzVector top_chi, w_chi, b_chi;
    TLorentzVector top_kin, w_kin, b_kin;         

    char pngname[200];
 
    sprintf(pngname, "njets(p_T^{jet}>20)");  
    TH1F *njets20 = new TH1F(pngname, pngname, 10 , 0,10);
    sprintf(pngname, "njets(p_T^{jet}>150)");  
    TH1F *njets150 = new TH1F(pngname, pngname, 10 , 0,10);

    if(all==1) {n1=0;    n2=bvpx->GetEntries()-2;  }  
    for(int iEvent=n1;iEvent<=n2;iEvent++){

      bvpx->GetEntry(iEvent);     bvpy->GetEntry(iEvent);        bvpz->GetEntry(iEvent);   bve->GetEntry(iEvent);
 
      for(int i = 0;i<hardon_px->size();i++){  
        particles.push_back( PseudoJet(hardon_px->at(i),hardon_py->at(i),hardon_pz->at(i),hardon_e->at(i)) );
      }

//    JetDefinition jet_def(cambridge_algorithm, R);
      JetDefinition jet_def(antikt_algorithm, R);
      ClusterSequence cs(particles, jet_def);
      jets = sorted_by_pt(cs.inclusive_jets(20));           
      hep_jets = sorted_by_pt(cs.inclusive_jets(150));
  
      njets20->Fill(jets.size() );
      njets150->Fill(hep_jets.size());
 
      if(jets.size()>=3){
  
      if( topjetreco_chi(jets, top_chi, w_chi, b_chi)==1 ) {
          B_Chi.push_back_momenta( b_chi );
          W_Chi.push_back_momenta( w_chi );
          Top_Chi.push_back_momenta( top_chi );

          if( abs( top_chi.M()- top_mass)<20.0 && abs( w_chi.M()- w_mass)<20.0 ){

            n_tagged_chi++;      

            Top_Chireco_tagged.push_back_momenta( top_chi );
            W_Chireco_tagged.push_back_momenta(w_chi);
            B_Chireco_tagged.push_back_momenta(b_chi);        
         }
       }//tagged

        if( topjetreco_kin(jets, top_kin, w_kin, b_kin)==1 ){
          B_Kin.push_back_momenta( b_kin );
          W_Kin.push_back_momenta( w_kin );
          Top_Kin.push_back_momenta( top_kin );

          if( abs( top_kin.M()- top_mass)<20.0 && abs( w_kin.M()- w_mass)<20.0 ){

            n_tagged_kin++;

            Top_Kinreco_tagged.push_back_momenta( top_kin );
            W_Kinreco_tagged.push_back_momenta(w_kin);
            B_Kinreco_tagged.push_back_momenta(b_kin);      
          }
        }//tagged
     }  //njets.>3

     for(unsigned ijet=0; ijet<hep_jets.size(); ijet++){      
       HEPTopTagger::HEPTopTagger tagger(hep_jets[ijet]);

       // Unclustering, Filtering & Subjet Settings
       tagger.set_max_subjet_mass(30.);
       tagger.set_mass_drop_threshold(0.8);
       tagger.set_filtering_R(0.3);
       tagger.set_filtering_n(5);
       tagger.set_filtering_minpt_subjet(30.); 

       // How to select among candidates
       tagger.set_mode(HEPTopTagger::TWO_STEP_FILTER);
  
       // Requirements to accept a candidate
       tagger.set_top_minpt(150); 
       tagger.set_top_mass_range(120, 240.); 
       tagger.set_fw(0.15); 

       // Run the tagger
       tagger.run();      

       // Look at output ie have a tag:
       if (tagger.is_tagged()) {
         Top_Hep.push_back_momenta(tagger.t().px(), tagger.t().py(), tagger.t().pz(), tagger.t().e());
         W_Hep.push_back_momenta(tagger.W().px(), tagger.W().py(), tagger.W().pz(), tagger.W().e());
         B_Hep.push_back_momenta(tagger.b().px(), tagger.b().py(), tagger.b().pz(), tagger.b().e());
         
         if( abs(tagger.t().m() - top_mass)<20.0 && abs(tagger.W().m() - w_mass)<20.0 ){
           n_tagged_hep++  ;

           Top_Hepreco_tagged.push_back_momenta(tagger.t().px(), tagger.t().py(), tagger.t().pz(), tagger.t().e());
           W_Hepreco_tagged.push_back_momenta(tagger.W().px(), tagger.W().py(), tagger.W().pz(), tagger.W().e());
           B_Hepreco_tagged.push_back_momenta(tagger.b().px(), tagger.b().py(), tagger.b().pz(), tagger.b().e());
         }
       }//tagged
     } //HepTT reconstruction

     particles.clear();  
     jets.clear();
     hep_jets.clear(); 
     cout<<"Processing data......"<<100.0*float(iEvent+1)/float(n2-n1+1)<<" % "<<"\r";

    } //loop over all events
    cout<<endl;

    B_Chi.setnames("b_Chi","");                 B_Chi.make_properties();                    
    W_Chi.setnames("w_Chi","");                 W_Chi.make_properties();                    
    Top_Chi.setnames("top_Chi","");             Top_Chi.make_properties();                  

    B_Kin.setnames("b_Kin","");               B_Kin.make_properties();                  
    W_Kin.setnames("w_Kin","");               W_Kin.make_properties();                   
    Top_Kin.setnames("top_Kin","");           Top_Kin.make_properties();                

    B_Hep.setnames("b_Hep","");               B_Hep.make_properties();                   
    W_Hep.setnames("w_Hep","");               W_Hep.make_properties();                   
    Top_Hep.setnames("top_Hep","");           Top_Hep.make_properties();                


    Top_Kinreco_tagged.setnames("TopKin_tagged","");      Top_Kinreco_tagged.make_properties();            
    W_Kinreco_tagged.setnames("WKin_tagged","");          W_Kinreco_tagged.make_properties();               
    B_Kinreco_tagged.setnames("BKin_tagged","");          B_Kinreco_tagged.make_properties();              
//
    Top_Hepreco.setnames("RecoTopHep","");                Top_Hepreco.make_properties();            
    W_Hepreco.setnames("RecoWHep","");                    W_Hepreco.make_properties();               
    B_Hepreco.setnames("RecoBHep","");                    B_Hepreco.make_properties();              

    Top_Hepreco_tagged.setnames("TopHep_tagged","");      Top_Hepreco_tagged.make_properties();            
    W_Hepreco_tagged.setnames("WHep_tagged","");          W_Hepreco_tagged.make_properties();               
    B_Hepreco_tagged.setnames("Hep_tagged","");           B_Hepreco_tagged.make_properties();              

    sprintf(pngname, "njets");
    TCanvas *c_njets = new TCanvas(pngname,pngname,1800,1000);
    c_njets->Divide(2,1);   
    c_njets->cd(1);
    njets20->SetLineWidth(2);    
    njets20->Draw("hist");
    c_njets->cd(2);  
    njets150->SetLineWidth(2);
    njets150->Draw("hist");
 
    int event_directory = mkdir("./qcd_jetanalysis", S_IRUSR | S_IWUSR | S_IXUSR);
    sprintf(pngname, "./qcd_jetanalysis/AR_%f",R);
    int cone_directory = mkdir(pngname, S_IRUSR | S_IWUSR | S_IXUSR);

    sprintf(pngname , "./qcd_jetanalysis/AR_%f/njets.png",R);
    c_njets->SaveAs(pngname);  

  } //run jet reco


       
     void EventSettings(string file_name, string tree_name, int iEvent, int fEvent, bool All ){
         filename = file_name;
         treename  = tree_name;     
         n1 = iEvent;   
         n2 = fEvent;
         all = All;      
     } 

     ~Analyser(){}
  
};



     void HistAngularDistance(ParticleProperty part1, ParticleProperty part2, TCanvas *canvas_dist, TH1F *hist_dist[], THStack *stack[], TLegend *leg[], int mass_index, float coneradius ){

      char name[50], outfilename[200], outfoldername[200];
      int z_mass[3] = {500,1000,2000};
     
      sprintf(name,"#Delta R_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[0]  = new TH1F(name,name, 50, 0, 10 );

      sprintf(name,"#Delta #Phi_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[1]  = new TH1F(name,name, 100, -6, 6 );

      sprintf(name,"#Delta y_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[2]  = new TH1F(name,name, 200, -10, 10 );

      sprintf(name,"#Delta #eta_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[3]  = new TH1F(name,name, 200, -10, 10 );

      float delR, delPhi, delRap, delEta;  

      for(int i=0;i<part1.prop[7].size();i++){          

        delRap = part1.prop[7][i] - part2.prop[7][i];
        delEta = part1.prop[6][i] - part2.prop[6][i];
        delPhi = part1.prop[9][i] - part2.prop[9][i];

        while(delPhi > M_PI) delPhi -= 2.0*M_PI;
        while(delPhi <= -M_PI) delPhi += 2.0*M_PI; 

        delR = sqrt( delRap*delRap + delPhi*delPhi );

        hist_dist[0]->Fill( delR ) ;
        hist_dist[1]->Fill( delPhi ) ;   
        hist_dist[2]->Fill( delRap ) ;
        hist_dist[3]->Fill( delEta ) ;

      }

      for(int iprop=0; iprop<4; iprop++){


        sprintf(name, "mean [mZ = %d] = %0.2f", z_mass[mass_index], hist_dist[iprop]->GetMean());
        leg[iprop]->SetTextSize(0.045);
        leg[iprop]->AddEntry(hist_dist[iprop],name,"l" );
        sprintf(name, "Entries [mZ = %d] = %0.0f", z_mass[mass_index], hist_dist[iprop]->GetEntries()-1);
        leg[iprop]->AddEntry(hist_dist[iprop],name,"l" );

      }

      for(int iprop=0; iprop<4; iprop++){
        FormatHist(hist_dist[iprop], mass_index+1, 2, false);
        stack[iprop]->Add(hist_dist[iprop]);    
      }

      if(mass_index==2){
        TLegendEntry *header[n_prop];
        canvas_dist->Divide(2,2);
  
  
        //for(int iprop=0; iprop<4; iprop++)
        for(int iprop=0; iprop<3; iprop++){
        canvas_dist->cd(iprop+1);
        leg[iprop]->SetHeader(hist_dist[iprop]->GetName());
        header[iprop] = (TLegendEntry*)leg[iprop]->GetListOfPrimitives()->First(); 
        header[iprop]->SetTextSize(.05);  
        header[iprop]->SetTextAlign(22);         
        stack[iprop]->Draw("hist nostack");
        leg[iprop]->Draw();
       }
                 
       sprintf(outfoldername, "./Z2ttbar_EventAnalysis");   
       sprintf(outfoldername, "%s/AR_%f",outfoldername,coneradius);
       sprintf(outfilename, "%s/%s.png",outfoldername,canvas_dist->GetName()); 
       canvas_dist->SaveAs(outfilename);  
     }

     for(int iprop=0;iprop<4;iprop++){
      hist_dist[iprop]->SetName("hist");      
     }
    } 












































