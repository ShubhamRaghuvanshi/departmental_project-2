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
#include "TLorentzVector.h"

using namespace std;

   const int n_prop=10;

class ParticleProperty{

  private:
    string particle_name, particle_surname; 
    string property[n_prop] = { "px",   "py",   "pz",  "e", "m",  "pt", "eta", "rap", "phi", "theta"};
  
  public:
    
    vector<double> prop[n_prop];
    
    ParticleProperty(){} 

    ~ParticleProperty(){}

    ParticleProperty(vector<double> &px, vector<double> &py, vector<double> &pz, vector<double> &e, string name, string append){

      particle_name = name;
      particle_surname = append;

      prop[0] = px; 
      prop[1] = py;
      prop[2] = pz;
      prop[3]  = e; 
 
      TLorentzVector v;

      for(int i=0;i<prop[0].size();i++){
        v.SetPxPyPzE( px[i],py[i],pz[i],e[i] );
        prop[4].push_back(v.M());            
        prop[5].push_back(v.Perp());
        prop[6].push_back(v.Eta());
        prop[7].push_back(v.Rapidity());    
        prop[8].push_back(v.Phi());
        prop[9].push_back(v.Theta());
      } 
    }

    void push_back_momenta( TLorentzVector v){
        prop[0].push_back(v.Px());            
        prop[1].push_back(v.Py());            
        prop[2].push_back(v.Pz());            
        prop[3].push_back(v.E());            
    }

    // constructs all properties given 4 momenta vectors
    void make_properties(){

      vector<double> px,py,pz,e;
      px = prop[0];
      py = prop[1];
      pz = prop[2];
      e  = prop[3];    

      TLorentzVector v;

      for(int i=0;i<prop[0].size();i++){
        v.SetPxPyPzE( px[i],py[i],pz[i],e[i] );
        prop[4].push_back(v.M());            
        prop[5].push_back(v.Perp());
        prop[6].push_back(v.Eta());
        prop[7].push_back(v.Rapidity());    
        prop[8].push_back(v.Phi());
        prop[9].push_back(v.Theta());
      } 
//cout<<"properies made for : "<<particle_name<<endl;
    }


    void push_back_momenta(double px, double py, double pz, double e){
        prop[0].push_back(px);            
        prop[1].push_back(py);            
        prop[2].push_back(pz);            
        prop[3].push_back(e);            
    }


    void setnames( string name, string append ){
      particle_name = name;
      particle_surname = append;
    }


    void SetFromFile(string filename, string treename, int treeindex){
  
      TFile f(filename.c_str(),"READ");
      TTree *tree; 
      f.GetObject(treename.c_str(),tree);
      
      vector<double> *p[4]={0,0,0,0};
      TBranch *b_p[4] = {0,0,0,0};

      char name[50];
      for(int iprop=0; iprop<4;iprop++){
        sprintf(name,"%s_%s%s",particle_name.c_str(), property[iprop].c_str(), particle_surname.c_str());  
        tree->SetBranchAddress(name,&p[iprop], &b_p[iprop]) ;
      }
        
      Long64_t tentry  = tree->LoadTree(treeindex);

      for(int iprop=0; iprop<4;iprop++)
        b_p[iprop]->GetEntry(tentry);

      for(int iprop=0; iprop<4;iprop++)
        b_p[iprop]->GetEntry(0);
  
      prop[0] = *p[0];
      prop[1] = *p[1];
      prop[2] = *p[2];
      prop[3] = *p[3];

      TLorentzVector v;

      for(int i=0;i<prop[0].size();i++){
        v.SetPxPyPzE( prop[0][i], prop[1][i],  prop[2][i], prop[3][i] );
        prop[4].push_back(v.M());            
        prop[5].push_back(v.Perp());
        prop[6].push_back(v.Eta());
        prop[7].push_back(v.Rapidity());    
        prop[8].push_back(v.Phi());
        prop[9].push_back(v.Theta());
      }
      f.Close();
    }

            
    void PrintVectorSizes(){
      cout<<setw(20)<<particle_name<<setw(20)<<particle_surname;
      for(int i=0;i<n_prop;i++)
        cout<<setw(8)<<prop[i].size();
      cout<<endl;
    }

    void PrintEntry(int n){
      for(int i=0;i<n_prop;i++)
        cout<<setw(15)<<prop[i][n]  ;
      cout<<endl;
    } 
     
    //iprop is property index: iprop<nprop
    TH1F* HistProp(int iprop, int particle_id){
    
     float  prop_low[n_prop]  =   {-1500.0,  -1500.0,  -2500.0,  0.0,    0.0,    0.0 ,  -10.0, -10.0, -8.0,   0 };
     float  prop_high[n_prop]  =  { 1500.0,   1500.0,   2500.0,  2500.0, 100.0,  1500.0,  10.0 , 10.0,  8.0, 6.0 };
      
    //z
    if(particle_id == 0 ){
          prop_low[4] = 450.0;   prop_high[4] = 2500.0;
          prop_low[3] = 300.0;   prop_high[3] = 4500.0;
          prop_low[5] = 0.0;     prop_high[5] = 600.0;
          prop_low[2] = -4500.0; prop_high[2] = 4500.0;          
    }
    //top
   else if(particle_id == 1 ){
        prop_low[4] = 150.0;   prop_high[4] = 200.0;
    }
    //w
    else if(particle_id == 2 ){
       prop_low[4] = 20.0;   prop_high[4] = 140.0; 
    }
     //b 
    else if(particle_id == 3 ){
      prop_low[4] = 0.0;   prop_high[4] = 150.0; 
    }

    else if(particle_id == 4 ){
      prop_low[4] = 100.0;   prop_high[4] = 250.0; 
    }
    // hadron
    else if(particle_id == 5 ){
      prop_low[0] = -100.0;   prop_high[0] = 100.0; 
      prop_low[1] = -100.0;   prop_high[1] = 100.0;
      prop_low[2] = -500.0;   prop_high[2] = 500.0;
      prop_low[3] = 0.0;      prop_high[3] = 300.0;
      prop_low[4] = 0.0;      prop_high[4] = 20.0;
      prop_low[5] = 0.0;      prop_high[5] = 50.0;
    }
    else if(particle_id == 6 ){
      prop_low[0] = -2500.0;   prop_high[0] = 2500.0; 
      prop_low[1] = -2500.0;   prop_high[1] = 2500.0;
      prop_low[2] = 0.0;   prop_high[2] = 2500.0;
      prop_low[3] = 0.0;      prop_high[3] = 2500.0;
      prop_low[4] = 0.0;      prop_high[4] = 20.0;
      prop_low[5] = 0.0;      prop_high[5] = 1500.0;
    }

    else{cout<<"new particle "<<endl;}

      int prop_bin[n_prop]  =  {100,  100,  100, 100, 100, 100 , 100, 100, 50, 50 };

      char name[50];
      sprintf(name,"%s_%s%s",particle_name.c_str(), property[iprop].c_str(), particle_surname.c_str());  
      TH1F *hist  = new TH1F(name,name, prop_bin[iprop], prop_low[iprop], prop_high[iprop] );
  
      for(int i=0;i<prop[iprop].size();i++){
        hist->Fill( prop[iprop][i] ) ;
      }

//cout<<hist->GetName()<<setw(20)<<particle_name[iprop]<<setw(20)<<property[iprop]<<setw(20)<<name<<endl;
      return hist;
    }


      
    string GetPropName(int iprop){ return property[iprop]; }
    string GetPartName(){ return particle_name; }
    string GetPartSurname(){ return particle_surname; }
};  //class definition



