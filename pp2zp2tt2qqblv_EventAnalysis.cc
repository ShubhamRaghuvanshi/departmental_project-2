
  #include "particleproperty.h"
  #include "Analyser.h"
  #include "Func.h" 

  using namespace std;
  using namespace fastjet;

  int main(){
  
    //######################################## Initialize Analyser object #####################

    vector<ParticleProperty> parton, tops;
       
    string partonnames[9] = {"z", "top", "topbar", "wp", "b", "wm", "bm", "q", "qb"};
    int z_mass[3] = {500, 1000, 2000};
    Analyser analyser[3];
    char filename[100],foldername[200], histname[100];
    
    vector<string> histlabel;
        
    for(int iFile=0; iFile<3; iFile++){
    
      sprintf(filename, "pp2zp2tt2qqblv_%d.root", z_mass[iFile] );
      sprintf(histname, " M_{Z'} = %d GeV", z_mass[iFile] );
      histlabel.push_back(histname);
      
      analyser[iFile].SetFile(filename);
      sprintf(foldername, "./%s/partons", analyser[iFile].FolderName().c_str());
      analyser[iFile].ReadPartons(&parton);
 
      for(int i=0; i<parton.size(); i++){
       parton[i].setnames(partonnames[i].c_str(), "");
      }
    
     
      tops.push_back(parton[1]);
      cout<<"total partons read :"<<parton.size()<<endl;
      
//      Draw1Histograms(parton[0],  foldername); 
//      DrawdelRvspT( parton[1], parton[3], parton[4] ,foldername);
//      DrawdelRvspT( parton[3], parton[7], parton[8] ,foldername);
//      DrawdelRvspT( parton[1], parton[7], parton[8] ,foldername);

      for(int i=0; i<parton.size(); i++){
      	parton[i].clear_properties();
      }
      parton.clear();    
    }
    DrawHistOnTop(tops, histlabel, foldername );
	
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

    return 0;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
