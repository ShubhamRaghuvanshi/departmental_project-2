
  #include "particleproperty.h"
  #include "Analyser.h"
  #include "Func.h" 
  #include <time.h>

  using namespace std;
  using namespace fastjet;

  int main(){
    clock_t t;
    //######################################## Initialize Analyser object #####################

    vector<ParticleProperty> parton, tops, top1;
      
    string partonnames[9] = {" Z'", "Top", "Topbar", "W+", "b", "W-", "bm", "q", "qb"};
    int z_mass[3] = {500, 1000, 2000};
    Analyser analyser[3];
    char filename[100],foldername[200];
    
        
    for(int iFile=0; iFile<3; iFile++){
    
      sprintf(filename, "pp2zp2tt2qqblv_%d.root", z_mass[iFile] );
      
      analyser[iFile].SetFile(filename);
      sprintf(foldername, "./%s/partons", analyser[iFile].FolderName().c_str());
      analyser[iFile].ReadPartons(&parton);
 
      for(int i=0; i<parton.size(); i++){
        parton[i].setnames(partonnames[i].c_str(), "");
      }
    
      tops.push_back(parton[1]);
      top1.push_back(parton[1]);
      
      cout<<"total partons read :"<<parton.size()<<endl<<endl;
     
      DrawHistograms(top1, 0, iFile , foldername);    top1.clear(); 
//      DrawdelRvspT( parton[1], parton[3], parton[4] ,foldername);
//      DrawdelRvspT( parton[3], parton[7], parton[8] ,foldername);
//      DrawdelRvspT( parton[1], parton[7], parton[8] ,foldername);

      for(int i=0; i<parton.size(); i++){
      	parton[i].clear_properties();
      }
      parton.clear();    
    }//files
    
    DrawHistograms(tops, 3, -1, foldername);    tops.clear() ;


//    DrawHistOnTop(tops, histlabel, foldername );
	
	//##########################################  Jet reconstruction ########################3

/*	
	  analyser.RecoJets(0.5, 1.0, true);
    sprintf(foldername, "./%s/jets", analyser.FolderName().c_str()); 

 
    Draw1Histograms(analyser.Top,  foldername);

    Draw1Histograms(analyser.Topkin,  foldername);
    Draw1Histograms(analyser.TopkinMatched,  foldername);
    Draw1Histograms(analyser.TopkinMatch,  foldername);

    Draw1Histograms(analyser.Topchi,  foldername);
    Draw1Histograms(analyser.TopchiMatched,  foldername);
    Draw1Histograms(analyser.TopchiMatch,  foldername);

    Draw1Histograms(analyser.Tophep,  foldername);
    Draw1Histograms(analyser.TophepMatched,  foldername);
    Draw1Histograms(analyser.TophepMatch,  foldername);

    DrawDividedHistograms(analyser.TopkinMatch, analyser.Top , foldername);
    DrawDividedHistograms(analyser.TopchiMatch, analyser.Top , foldername);
    DrawDividedHistograms(analyser.TophepMatch, analyser.Top , foldername);
    
    
    cout<<"top size : "<<analyser.Top.prop[0].size()<<endl;
*/
    t = clock() - t;
    cout<<"Time of execution (seconds ): "<<t/float(CLOCKS_PER_SEC)<<endl;
    return 0;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
