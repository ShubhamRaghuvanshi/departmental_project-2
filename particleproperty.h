

#ifndef PARTICLEPROPERTY_H
#define PARTICLEPROPERTY_H

#include <iostream>
#include <fstream>
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


using namespace std;
  

const int n_prop =11;
class ParticleProperty{

  private:
    string particle_name, particle_surname; 
    
    string property[n_prop] = { "px",   "py",   "pz",  "E", "M",  "pT", "eta", "y", "phi", "theta", "jetsize"};
    
  public:
    
    string propXaxis[n_prop] = { "P_{x}", "P_{y}", "P_{z}", "E", "M_{0}",  "p_{T}", "#eta", "y", "#phi", "#theta", "jetsize" };
     
    vector<double> prop[n_prop];
    
    ParticleProperty(){}; 

    ~ParticleProperty(){};

    ParticleProperty(vector<double> &px, vector<double> &py, vector<double> &pz, vector<double> &e, string name, string append);

    void push_back_momenta( TLorentzVector v);
	
		void push_back_momenta(double px, double py, double pz, double e);
    void push_back_momenta( vector<TLorentzVector> v);

    TLorentzVector GetLorentzVector(int i);
    vector<TLorentzVector> GetLorentzVector();

    int push_back_momenta(vector<double> px, vector<double> py, vector<double> pz, vector<double> e);
        
    void add_momenta(ParticleProperty p1, ParticleProperty p2);    
    void make_properties();
    void clear_properties();
    void push_property(int iprop, float val);

    void setnames( string name, string append );

    void SetFromFile(string filename, string treename, int treeindex); 
            
    void PrintVectorSizes();

    void PrintEntry(int n);
     
    //iprop is property index: iprop<nprop
    TH1F* HistProp(int iprop, int def);
    
    string GetPropName(int iprop); 
    string GetPartName(); 
    string GetPartSurname();
};  //class definition

#endif















