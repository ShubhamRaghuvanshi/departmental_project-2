
  #include "particleproperty.h"
  #include "Analyser.h"
  #include "Func.h" 
  #include <time.h>

  using namespace std;
  using namespace fastjet;

  int main(){
    clock_t t;
    //######################################## Initialize Analyser object #####################

    Analyser analyser;
    vector<ParticleProperty> leadingjet;
    char filename[100],foldername[200];



    sprintf(filename, "QCD20.root" );
    analyser.SetFile(filename);    
    analyser.RecoJets(0.5, 1.0, false);
	  sprintf(foldername, "./%s/jets", analyser.FolderName().c_str());
  	  
  	leadingjet.push_back(analyser.LeadingJet);  
    DrawHistograms(leadingjet, 0 , -1, foldername);	  
	  leadingjet.clear();
	  
	//##########################################  Jet reconstruction ########################3

// pT and mass of the leading jet > 200
// reco mass HepTT , chi
// size ?
// 

    t = clock() - t;
    cout<<"Time of execution (seconds ): "<<t/float(CLOCKS_PER_SEC)<<endl;
    return 0;
  }
  
  
  // /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/SoftDrop.cc 
  // /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/RecursiveSymmetryCutBase.cc
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
