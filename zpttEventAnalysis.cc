

  #include "particleproperty.h"
  #include "Analyser.h"
  #include "Func.h" 

  using namespace std;
  using namespace fastjet;

  int main(){
  
    //######################################## Initialize Analyser object #####################
  
  	//string filename = "zptt.root";
  	//string filename = "pp2tttt.root";
  	
  	
  	string filename = "pp2zp2tt2qqblv_500.root";
		char foldername[50];


    Analyser analyser(filename);
    cout<<"write folder : "<<analyser.FolderName()<<endl;
		sprintf(foldername, "./%s/partons", analyser.FolderName().c_str()); 	

    string partonnames[9] = {"z, top, topbar, wp, b, wm, bm, q, qb  "};     
//    string partonnames[4] = {"TopFromZ", "TopbarFromZ", "t3", "t4"};     
   // string partonnames[4] = {"Top[0]", "Topbar[0]", "Top[1]", "Topbar[1]"};
    
    
    //########################################### Parton Properties ############################   
    vector<ParticleProperty> parton;
    ParticleProperty zp;
   	analyser.ReadPartons(&parton);
    cout<<"total partons read :"<<parton.size()<<endl;
    
    //if(parton.size() != 4) {cout<<"parton numbers dont match, aborting code "<<endl; return 1;}
    if(parton.size() != 9) {cout<<"parton numbers dont match, aborting code "<<endl; return 1;}
    
    for(int i=0; i<parton.size(); i++){
    	parton[i].setnames(partonnames[i].c_str(), "");
    }
    zp.add_momenta(parton[0], parton[1]);
    
		zp.setnames("zp","");
		
		for(int i=0; i<parton.size(); i++){
			 Draw1Histograms(parton[i], 1, foldername);	
		}

		Draw4Histograms(parton, 1, foldername); 
		Draw1Histograms(zp, 0, foldername); 
	
	//##########################################  Jet reconstruction ########################3
	  analyser.RecoJets();
	  sprintf(foldername, "./%s/jets", analyser.FolderName().c_str()); 		
	  Draw1Histograms(analyser.Tophep, 4, foldername);
	  Draw1Histograms(analyser.Topchi, 4, foldername);
	  Draw1Histograms(analyser.Zphep, 0, foldername);	

    return 0;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
