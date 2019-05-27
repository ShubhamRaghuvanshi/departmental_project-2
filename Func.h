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

  const double delta_mw = 15.0;
  const double delta_mtop = 10.0;

  double delR( TLorentzVector v1, TLorentzVector v2 );
  
  int topjetreco_kin( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B, ParticleProperty *TaggedTop );
  int topjetreco_chi( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B,ParticleProperty *TaggedTop );										
  int topjetreco_hep( vector<PseudoJet> jets, ParticleProperty *Tjet, ParticleProperty *Wjet, ParticleProperty *Bjet, ParticleProperty *TaggedTop);
  
     
  void DrawdelR( ParticleProperty p1, ParticleProperty p2, string foldername );
  void FormatHist(TH1F *hist, int linecolor, int linewidth, bool stat);
  int DrawHistograms(vector<ParticleProperty> particle, int drawoption, int mass_index,  string foldername);
   

//  void DrawdelRvspT(ParticleProperty top, ParticleProperty w, ParticleProperty b, string foldername);
  
#endif



























