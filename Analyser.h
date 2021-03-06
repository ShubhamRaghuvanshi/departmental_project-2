#ifndef ANALYSER_H
#define ANALYSER_H

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
#include "particleproperty.h"
#include "./HEPTopTagger.hh"
#include "Func.h"

using namespace std;
  
using namespace fastjet;


class Analyser{

  private:
    string filename;
  public:
 
    ParticleProperty Tophep, Whep, Bhep, Zphep;
    ParticleProperty Topchi, Wchi, Bchi;
    ParticleProperty Topkin, Wkin, Bkin;

    ParticleProperty TophepTagged, TopchiTagged, TopkinTagged;
    ParticleProperty TophepMatched,TopchiMatched, TopkinMatched; 
    ParticleProperty TophepMatch, TopchiMatch,TopkinMatch;
    
    ParticleProperty LeadingJet, top1, top2, top3, top4;
        
    ParticleProperty Top, FatTop;

        
    Analyser(){
      
      Zphep.setnames("Zphep","");
      Top.setnames("top","");
      FatTop.setnames("Fattop","");
      
      top1.setnames("top1_{hepTT}", "");
      top2.setnames("top2_{hepTT}", "");
      top3.setnames("top1_{#chi^{2}}", "");
      top4.setnames("top2_{#chi^{2} }", "");
      
      Topkin.setnames("Topkin","");
      Wkin.setnames("Wkin","");
      Bkin.setnames("Bkin","");
      
      Topchi.setnames("Topchi","");
      Wchi.setnames("Wchi","");
      Bchi.setnames("Bchi","");

      Tophep.setnames("Tophep","");
      Whep.setnames("Whep","");
      Bhep.setnames("Bhep","");

      TopkinTagged.setnames("top_{kin}","");  
      TopchiTagged.setnames("top_{chi}","");    
      TophepTagged.setnames("top_{hepTT}","");  

      TopkinMatched.setnames("top(kin)","");  
      TopchiMatched.setnames("top(#chi^{2})","");    
      TophepMatched.setnames("top(hepTT)","");  

      TopkinMatch.setnames("top_{parton}","");  
      TopchiMatch.setnames("top_{parton}","");    
      TophepMatch.setnames("top_{parton}","");  
        
      LeadingJet.setnames("LeadingJet","");  
        
    };
    ~Analyser(){};

    Analyser(string file_name);

    void SetFile(string file_name);

    int ReadPartons( vector<ParticleProperty> *pp);
   // int RecoJetsztt( );
    int ReadHadrons(ParticleProperty *p );
    int RecoJets(float R, float fatR, bool match);
    int Reco4Tops(float R, float fatR);
    int AKJets(float R, vector<ParticleProperty> *QCDjets);
    
    string FolderName();
};  //analyser class


#endif












































