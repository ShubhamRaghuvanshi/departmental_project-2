
  #include "particleproperty.h"
  #include "Analyser.h"
  #include "Func.h" 
  #include <time.h>

  using namespace std;
  using namespace fastjet;

  int main(){
    clock_t t;
    //######################################## Initialize Analyser object #####################

    vector<ParticleProperty> parton, tops, zps, top1;
      
    string partonnames[9] = {" Z'", "Top", "Topbar", "W+", "b", "W-", "bm", "q", "qb"};
    int z_mass[3] = {500, 1000, 2000};
    Analyser analyser[3];
    char filename[100],foldername[200];
    
 /*       
    for(int iFile=0; iFile<3; iFile++){
    
      sprintf(filename, "pp2zp2tt2qqblv_%d.root", z_mass[iFile] );
      
      analyser[iFile].SetFile(filename);
      sprintf(foldername, "./%s/partons", analyser[iFile].FolderName().c_str());
      analyser[iFile].ReadPartons(&parton);
 
      for(int i=0; i<parton.size(); i++){
        parton[i].setnames(partonnames[i].c_str(), "");
      }

      zps.push_back(parton[0]);    
      tops.push_back(parton[1]);      
      top1.push_back(parton[1]);
      
      cout<<"total partons read :"<<parton.size()<<endl<<endl;
     
      DrawHistograms(top1, 0, iFile , foldername);    top1.clear(); 

      for(int i=0; i<parton.size(); i++){
      	parton[i].clear_properties();
      }
      parton.clear();    
    }//files
    
    DrawHistograms(tops, 3, -1, foldername);    tops.clear() ;
    DrawHistograms(zps, 3, -1, foldername);     zps.clear() ;
*/
	
	//##########################################  Jet reconstruction ########################3

	
	   int mass_index=2; 
     sprintf(filename, "pp2zp2tt2qqblv_%d.root", z_mass[mass_index] );
     analyser[mass_index].SetFile(filename);

	   analyser[mass_index].RecoJets(0.5, 1.0, true);
     sprintf(foldername, "./%s/jets", analyser[mass_index].FolderName().c_str()); 

     cout<<"top size : "<<setw(20)<<analyser[mass_index].Top.prop[0].size()<<setw(20)<<analyser[mass_index].Top.prop[1].size()<<setw(20)<<analyser[mass_index].Top.prop[3].size()<<endl;

/*

     vector<ParticleProperty> recoeff, recoeff_kin, recoeff_chi, recoeff_hep ;   
     vector<ParticleProperty> rawtop;
     
     rawtop.push_back(analyser[mass_index].Topkin);
     rawtop.push_back(analyser[mass_index].Topchi);
     rawtop.push_back(analyser[mass_index].Tophep); 

     recoeff.push_back(analyser[mass_index].TopkinMatched); 
     recoeff.push_back(analyser[mass_index].TopchiMatched); 
     recoeff.push_back(analyser[mass_index].TophepMatched); 
     recoeff.push_back(analyser[mass_index].Top); 

     recoeff_kin.push_back(analyser[mass_index].TopkinMatched); 
     recoeff_chi.push_back(analyser[mass_index].TopchiMatched); 
     recoeff_hep.push_back(analyser[mass_index].TophepMatched); 
     recoeff_kin.push_back(analyser[mass_index].Top);
     recoeff_chi.push_back(analyser[mass_index].Top);
     recoeff_hep.push_back(analyser[mass_index].Top); 
     
     DrawHistograms(rawtop, 0 , mass_index, foldername);   rawtop.clear();   
     DrawHistograms(recoeff_kin, 1 , mass_index, foldername);   recoeff_kin.clear(); 
     DrawHistograms(recoeff_chi, 1 , mass_index, foldername);   recoeff_chi.clear();
     DrawHistograms(recoeff_hep, 1 , mass_index, foldername);   recoeff_hep.clear(); 
     DrawHistograms(recoeff, 2 , mass_index, foldername);    recoeff.clear() ;
 */
      
     for(int i=0; i< analyser[mass_index].Top.prop[5].size(); i++)
      cout<<setw(20)<<i<<setw(20)<<analyser[mass_index].Top.prop[5][i]<<endl;
     
     
    t = clock() - t;
    cout<<"Time of execution (seconds ): "<<t/float(CLOCKS_PER_SEC)<<endl;
    return 0;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
