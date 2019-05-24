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
//#include "./HEPTopTagger.hh"
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
        
    ParticleProperty Top;
        
    Analyser(){
      
      Zphep.setnames("Zphep","");
      Top.setnames("Top","");;
      
      Topkin.setnames("Topkin","");
      Wkin.setnames("Wkin","");
      Bkin.setnames("Bkin","");
      
      Topchi.setnames("Topchi","");
      Wchi.setnames("Wchi","");
      Bchi.setnames("Bchi","");

      Tophep.setnames("Tophep","");
      Whep.setnames("Whep","");
      Bhep.setnames("Bhep","");

      TopkinTagged.setnames("TopkinTagged","");  
      TopchiTagged.setnames("TopchiTagged","");    
      TophepTagged.setnames("TophepTagged","");  

      TopkinMatched.setnames("MatchedTopJet_kin","");  
      TopchiMatched.setnames("MatchedTopJet_chi","");    
      TophepMatched.setnames("MatchedTopJet_hep","");  

      TopkinMatch.setnames("MatchedTop_kin","");  
      TopchiMatch.setnames("MatchedTop_chi","");    
      TophepMatch.setnames("MatchedTop_hep","");  
        
    };
    ~Analyser(){};

    Analyser(string file_name);

    void SetFile(string file_name);

    int ReadPartons( vector<ParticleProperty> *pp);
   // int RecoJetsztt( );
    int RecoJets(float R, float fatR, bool match);
    string FolderName();
};  //analyser class


#endif












































