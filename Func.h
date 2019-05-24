#ifndef FUNC_H
#define FUNC_H

  #include <iostream>
  #include <fstream>
  #include <sys/stat.h>
  #include <sys/types.h>  
  #include "TGraph.h"
  #include "TTree.h"
  #include "TFile.h"
  #include "TH1.h"
  #include "TH2.h"
  #include "TCanvas.h"
  #include "TStyle.h"
  #include "TGaxis.h"
  #include "TLorentzVector.h"
  #include "fastjet/ClusterSequence.hh"
  #include "fastjet/PseudoJet.hh"
  #include "HEPTopTagger.hh"
  #include "particleproperty.h"
  

  using namespace std;
  using namespace fastjet;

  const double m_Zp = 1000;
  const double fatjet_pt = 200;

  const double m_w = 80.0;
  const double m_top =172.0;
  //only the rlative weights matter, keep them to be the decay widths. 
  const double sigma_w = 2.0;
  const double sigma_top = 1.5;


  double delR( TLorentzVector v1, TLorentzVector v2 );
  
  int topjetreco_kin( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B,ParticleProperty *TaggedTop );
  int topjetreco_chi( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B,ParticleProperty *TaggedTop );										
  int topjetreco_hep( vector<PseudoJet> jets, ParticleProperty *Tjet, ParticleProperty *Wjet, ParticleProperty *Bjet, ParticleProperty *TaggedTop);
  
  int topjetreco_hepTT( vector<PseudoJet> jets, ParticleProperty *Tjet, ParticleProperty *Wjet, ParticleProperty *Bjet, ParticleProperty *Zpjet, int &t1,int &t2);
   
  void DrawdelR( ParticleProperty p1, ParticleProperty p2, string foldername );
  void FormatHist(TH1F *hist, int linecolor, int linewidth, bool stat);
void hist_properties(TCanvas *c[], TH1F *h1[], TH1F *h2[], THStack *stack1[], THStack *stack2[], TLegend *leg1[], TLegend *leg2[], int mass_index, float coneradius );
void hist_properties(TCanvas *c[], TH1F *h1[], THStack *stack1[], TLegend *leg1[], int mass_index, float coneradius,int n_hadronictop,int n_tagged  );
  void HistAngularDistance(ParticleProperty part1, ParticleProperty part2, TCanvas *canvas_dist, TH1F *hist_dist[], THStack *stack[], TLegend *leg[], int mass_index, float coneradius );
	void Draw1Histograms(ParticleProperty particle,  string foldername);
  int DrawHistOnTop(vector<ParticleProperty> p, vector<string> histlabel,  string foldername  );	
  void Draw4Histograms(vector<ParticleProperty> particle,  string foldername);
  void DrawDividedHistograms(ParticleProperty particle1, ParticleProperty particle2,  string foldername);
  void DrawdelRvspT(ParticleProperty top, ParticleProperty w, ParticleProperty b, string foldername);
  
#endif



























