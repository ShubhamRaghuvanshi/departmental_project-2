
  #include "particleproperty.h"
  #include "Analyser.h"
  #include "Func.h" 
  #include <time.h>

  using namespace std;
  using namespace fastjet;

  int main(){
    clock_t t;
    //######################################## Initialize Analyser object #####################

    Analyser analyser[3];
    int z_mass[3] = {500, 1000, 2000};
    char filename[100],foldername[200];




/*
    vector<vector<TLorentzVector>> hadrons;
    vector<ParticleProperty> parton, tops, zps, top1;
    TLorentzVector v_temp;  
    string partonnames[9] = {" Z'", "Top", "Topbar", "w", "b", "W-", "bm", "q", "qb"};

    


    
//    for(int iFile=0; iFile<3; iFile++){
      int iFile = 1;  
      sprintf(filename, "pp2zp2tt2qqblv_%d.root", z_mass[iFile] );
      
      analyser[iFile].SetFile(filename);
      
      sprintf(foldername, "./%s/partons", analyser[iFile].FolderName().c_str());
            
      analyser[iFile].ReadPartons(&parton);      
      for(int i=0; i<parton.size(); i++){
        parton[i].setnames(partonnames[i].c_str(), "");
      }

//
      analyser[iFile].ReadHadrons(&parton[1]);      
      cout<<"size : "<<parton[1].prop[0].size()<<setw(20)<<hadrons.size()<<setw(20)<<parton[1].prop[10].size()<<endl;

      Drawone(parton[1], 10, foldername);
      Draw1v2(parton[1], parton[1], 5, 10, foldername);
//      
      
      
      DrawdelRvspT( parton[1], parton[3], parton[4], foldername );
      DrawdelRvspT( parton[3], parton[7], parton[8], foldername );

      

  //    zps.push_back(parton[0]);    
  //    tops.push_back(parton[1]);      
  //    top1.push_back(parton[1]);
      
      
      
      cout<<"total partons read :"<<parton.size()<<endl<<endl;
     
  //    DrawHistograms(top1, 0, iFile , foldername);    top1.clear(); 

      for(int i=0; i<parton.size(); i++){
      	parton[i].clear_properties();
      }
      parton.clear();    
 //   }//files
    
    //  
    
   // DrawHistograms(top1, 0, -1, foldername);    tops.clear() ;
//    DrawHistograms(zps, 3, -1, foldername);     zps.clear() ;
*/	
	
	
	//##########################################  Jet reconstruction ########################3

	
	   int mass_index; 
     vector<ParticleProperty> toprecokin, toprecochi, toprecohep, top	; 
     vector<ParticleProperty> BigTop;
     for(mass_index =0; mass_index<3; mass_index++){
	   	   	   	   
       sprintf(filename, "pp2zp2tt2qqblv_%d.root", z_mass[mass_index] );
       analyser[mass_index].SetFile(filename);

	     analyser[mass_index].RecoJets(0.5, 1.0, true);
        
       toprecokin.push_back(analyser[mass_index].TopkinMatched);
       
       toprecochi.push_back(analyser[mass_index].TopchiMatched);
       toprecohep.push_back(analyser[mass_index].TophepMatched);         
       top.push_back(analyser[mass_index].Top);       
 
    } 
     
    for(mass_index=0; mass_index<3; mass_index++) 
      BigTop.push_back(toprecokin[mass_index]); 
    for(mass_index=0; mass_index<3; mass_index++) 
      BigTop.push_back(toprecochi[mass_index]); 
    for(mass_index=0; mass_index<3; mass_index++) 
      BigTop.push_back(toprecohep[mass_index]); 
    for(mass_index=0; mass_index<3; mass_index++) 
      BigTop.push_back(top[mass_index]); 
    
    sprintf(foldername, "./%s/jets", analyser[2].FolderName().c_str());     
    DrawHistograms(BigTop, 1, -1, foldername); 
     
    
    
    toprecokin.clear();
    toprecochi.clear();
    toprecohep.clear();
    top.clear();
    BigTop.clear();

    t = clock() - t;
    cout<<"Time of execution (seconds ): "<<t/float(CLOCKS_PER_SEC)<<endl;
    return 0;
  }
  
  
  // /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/SoftDrop.cc 
  // /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/RecursiveSymmetryCutBase.cc
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
