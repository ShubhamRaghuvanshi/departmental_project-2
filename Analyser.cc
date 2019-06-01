
#include "Analyser.h"
#include "Func.h"

using namespace std;
using namespace fastjet;

float top_mass = 172.0;
float w_mass = 80;


    Analyser::Analyser(string file_name){
    
      filename = file_name;
      char name[50];
      file_name.erase( remove(file_name.begin(), file_name.end(), '.'), file_name.end());
      
      sprintf(name,"./%s", file_name.c_str());      
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);
 
      cout<<"output folder : "<<file_name<<endl;      
 
      sprintf(name,"%s/partons", file_name.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);


      sprintf(name,"%s/jets", file_name.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);  

      
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
 
    } 

     void Analyser::SetFile(string file_name){

      filename = file_name;
      file_name.erase( remove(file_name.begin(), file_name.end(), '.'), file_name.end());
 
      char name[30];

      sprintf(name,"./%s", file_name.c_str());      
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);
      cout<<"output folder : "<<file_name<<endl;      

      sprintf(name,"%s/partons", file_name.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);
      sprintf(name,"%s/jets", file_name.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);  
     } 

    string Analyser::FolderName(){
      string file_name = filename;
      file_name.erase( remove(file_name.begin(), file_name.end(), '.'), file_name.end());
      return file_name;
    }
    
    
    int Analyser::givehadrons(int iEvent, vector<TLorentzVector> *particles){
      TFile *f = new TFile(filename.c_str(), "READ");    
      TTree *tree;
      
      f->GetObject("hadrons",tree);          
                   
      int nPartons = tree->GetNbranches();
      
      vector<double>  *px=0, *py=0, *pz=0, *e=0;
      TBranch  *bpx=0, *bpy=0, *bpz=0, *be=0;
    
      tree->SetBranchAddress("px",&px,&bpx);
      tree->SetBranchAddress("py",&py,&bpy);
      tree->SetBranchAddress("pz",&pz,&bpz);
      tree->SetBranchAddress("e",&e ,&be);   

      Long64_t tentry  = tree->LoadTree(0);
      bpx->GetEntry(tentry);       bpy->GetEntry(tentry);       bpz->GetEntry(tentry);   be->GetEntry(tentry);
    
      if(px->size() != py->size() || px->size() != pz->size() || px->size() != e->size()  ) return -666;
         
        bpx->GetEntry(iEvent);       bpy->GetEntry(iEvent);       bpz->GetEntry(iEvent);   be->GetEntry(iEvent);  
        
        for(int i=0; i< px->size() ; i++)
		    	particles->push_back(  TLorentzVector(px->at(i), py->at(i), pz->at(i), e->at(i))  );
        
        return 0;
      }
 	    
    

    int Analyser::ReadHadrons( ParticleProperty *p){
cout<<"reading hadrons"<<endl;
      TFile *f = new TFile(filename.c_str(), "READ");    
      TTree *tree;
      
      f->GetObject("hadrons",tree);          
                   
      int nPartons = tree->GetNbranches();
      
      vector<double>  *px=0, *py=0, *pz=0, *e=0;
      TBranch  *bpx=0, *bpy=0, *bpz=0, *be=0;
    
      tree->SetBranchAddress("px",&px,&bpx);
      tree->SetBranchAddress("py",&py,&bpy);
      tree->SetBranchAddress("pz",&pz,&bpz);
      tree->SetBranchAddress("e",&e ,&be);   

      Long64_t tentry  = tree->LoadTree(0);
      bpx->GetEntry(tentry);       bpy->GetEntry(tentry);       bpz->GetEntry(tentry);   be->GetEntry(tentry);
    
      if(px->size() != py->size() || px->size() != pz->size() || px->size() != e->size()  ) return -666;
      
      vector<TLorentzVector> hadrons;
      
     int nEvent = bpx->GetEntries()-1;
   //  nEvent = 50000;
      if(p->prop[0].size() != nEvent) {cout<<"unequal sizes"<<endl; return -888;}
      for(int iEvent =0; iEvent < nEvent ; iEvent++){
     // cout<<setw(20)<<iEvent<<setw(20)<<px->size()<<endl;
      
        bpx->GetEntry(iEvent);       bpy->GetEntry(iEvent);       bpz->GetEntry(iEvent);   be->GetEntry(iEvent);  

        cout<<"Processing hadron data......"<<100.0*float(iEvent+1)/float(nEvent+1)<<" % "<<"\r"; 
     
        p->prop[10].push_back( sizeofjet( p->GetLorentzVector(iEvent) , px, py, pz, e ));

       hadrons.clear();
       
      } cout<<endl;
 		 
    }

    int Analyser::ReadPartons( vector<ParticleProperty> *pp){

      TFile *f = new TFile(filename.c_str(), "READ");    
      TTree *tree, *tree2;
      f->GetObject("partons",tree);
                
      TBranch *PartonBranch = 0;
      TLorentzVector *PartonVector=0;    

      //this is needed to be set in the root file in the form "parton[integer]";
      char PartonBranchName[20];
      
      int nPartons = tree->GetNbranches();
                    
      cout<<"expected number of partons : "<<nPartons<<endl;    
          
      if(nPartons < 1 || nPartons > 100)  return -1;
     
      for(int part=0; part< nPartons; part++){
  cout<<"Reading parton : "<<part<<endl;    
        ParticleProperty pp_temp;      
        sprintf(PartonBranchName, "parton[%d]", part);
        tree->SetBranchAddress(PartonBranchName, &PartonVector, &PartonBranch);
        
        Long64_t tentry  = tree->LoadTree(0);
        PartonBranch->GetEntry(tentry); 
        int nEvent = PartonBranch->GetEntries()-1;
        //nEvent = 50000;
  
        if(part ==0 )
          cout<<"Number of Events in the file : "<<nEvent<<endl;
  
        for(int iEvent=0; iEvent<nEvent; iEvent++ ){
          PartonBranch->GetEntry(iEvent);
          pp_temp.push_back_momenta(*PartonVector); 
        //  cout<<"loading events......"<< 100.0*float(iEvent + part*nEvent +1)/float( nEvent*nPartons +1 ) <<" % "<<"\r"; 
        }
      //  cout<<"out"<<endl;
        pp_temp.make_properties(); 
        pp->push_back(pp_temp);
      } // for all partons
     // cout<<endl;
     
      return 0;  
    } //read partons


    int Analyser::RecoJets(float R, float fatR, bool match){
      
      float match_dist = 0.3;
    
      TFile f(filename.c_str(), "READ");    
      TTree *tree, *tree2;
      f.GetObject("hadrons", tree);
      f.GetObject("partons", tree2);
            
      TBranch  *bpx=0, *bpy=0, *bpz=0, *be=0, *btop=0;
      vector<double>  *px=0, *py=0, *pz=0, *e=0;
      TLorentzVector *top=0;
    
      tree->SetBranchAddress("px",&px,&bpx);
      tree->SetBranchAddress("py",&py,&bpy);
      tree->SetBranchAddress("pz",&pz,&bpz);
      tree->SetBranchAddress("e",&e ,&be);   
      tree2->SetBranchAddress("parton[1]", &top, &btop);
            
      Long64_t tentry  = tree->LoadTree(0);
      bpx->GetEntry(tentry);       bpy->GetEntry(tentry);       bpz->GetEntry(tentry);   be->GetEntry(tentry);

      Long64_t tentry2  = tree2->LoadTree(0);
      btop->GetEntry(tentry2);

      char name[100];
      sprintf(name, "nJets20");  
      TH1F *njets20 = new TH1F(name, name, 20 , 0, 20);
      njets20->GetXaxis()->SetTitle("N_{jets}(p_{T} > 20 GeV)");
      FormatHist(njets20, 1, 2, false );

      sprintf(name, " nJets200");  
      TH1F *njets200 = new TH1F(name, name, 7 , 0, 7);
      njets200->GetXaxis()->SetTitle("N_{jets}(p_{T} > 200 GeV)");
      FormatHist(njets200, 1, 2, false );

      sprintf(name, "nTaggedkin");      
      TH1F *histTagged_kin = new TH1F(name, name, 7 , 0, 7);
      histTagged_kin->GetXaxis()->SetTitle("nTagged_kin");
      FormatHist(histTagged_kin, 1, 2, false );

      sprintf(name, "nTaggedchi");      
      TH1F *histTagged_chi = new TH1F(name, name, 7 , 0, 7);
      histTagged_chi->GetXaxis()->SetTitle("nTagged_{#chi^{2}}");
      FormatHist(histTagged_chi, 1, 2, false );

      sprintf(name, "nTaggedhep");      
      TH1F *histTagged_hep = new TH1F(name, name, 7 , 0, 7);
      histTagged_hep->GetXaxis()->SetTitle("nTagged_HepTT");
      FormatHist(histTagged_hep, 1, 2, false );


      vector<PseudoJet> particles, jets, fatjets; 
      int nEvent = bpx->GetEntries()-1;
      nEvent = 20000;
          
      cout<<"TOTAL ENTRIES : "<<bpx->GetEntries()-1<<endl;    
      if(nEvent > bpx->GetEntries()){
        cout<<"More than you can handle. "<<endl;
        return -777;
      }
                
              
      int nJet=0, nFatjet=0, nThreeJetEvent=0;    
      for(int iEvent =0; iEvent < nEvent ; iEvent++){

        int iTagged_kin=0, iTagged_chi=0, iTagged_hep=0; 

        bpx->GetEntry(iEvent);       bpy->GetEntry(iEvent);       bpz->GetEntry(iEvent);   be->GetEntry(iEvent);  btop->GetEntry(iEvent);
        
       // cout<<"hadorns "<<px->size()<<endl;
        
        
        Top.push_back_momenta(*top);
        
        if(px->size() != py->size() || px->size() != pz->size() || px->size() != e->size()  ) return -666;
      
        for(int i=0; i< px->size() ; i++)
		    	particles.push_back(  PseudoJet(px->at(i), py->at(i), pz->at(i), e->at(i))  );
		  
		    JetDefinition jet_def(antikt_algorithm, R);  
		    JetDefinition hep_def(cambridge_algorithm, fatR);
		                 
        ClusterSequence cs_jet(particles, jet_def);
        ClusterSequence cs_hep(particles, hep_def);
        
        jets    = sorted_by_pt(cs_jet.inclusive_jets(20));           
        fatjets = sorted_by_pt(cs_hep.inclusive_jets(200));
          
       // jets = UseSoftDrop(jets, R);
  
        if( jets.size() !=0)
          njets20->Fill(jets.size() );
        if( fatjets.size() !=0)
          njets200->Fill(fatjets.size());

        if(jets.size() > 2){  
          nJet+=jets.size();
          nThreeJetEvent++;
        }

        if(fatjets.size() > 0)  
         nFatjet+=fatjets.size();
         
        int lastsize_kin = TopkinTagged.prop[0].size();
        int lastsize_chi = TopchiTagged.prop[0].size();
        int lastsize_hep = TophepTagged.prop[0].size();
                        
        iTagged_kin =topjetreco_kin(jets, &Topkin, &Wkin, &Bkin, &TopkinTagged);
        iTagged_chi =topjetreco_chi(jets, &Topchi, &Wchi, &Bchi, &TopchiTagged);
        iTagged_hep =topjetreco_hep(fatjets, &Tophep, &Whep, &Bhep, &TophepTagged);
          
        if(match){      
          iTagged_kin=0; iTagged_chi=0; iTagged_hep=0;  

          for(int i = lastsize_kin ; i< TopkinTagged.prop[0].size(); i++){          
            TopkinMatch.push_back_momenta(*top);
            
            if( delR(*top, TopkinTagged.GetLorentzVector(i) ) < match_dist ){
              TopkinMatched.push_back_momenta( TopkinTagged.GetLorentzVector(i) );         
              iTagged_kin++;     
            }                   


          }
          
          for(int i = lastsize_chi ; i< TopchiTagged.prop[0].size(); i++){          
            TopchiMatch.push_back_momenta(*top);
            if( delR(*top, TopchiTagged.GetLorentzVector(i) ) < match_dist ){
              TopchiMatched.push_back_momenta( TopchiTagged.GetLorentzVector(i) );          
              iTagged_chi++;    
            }                   
          }
          
          for(int i = lastsize_hep ; i< TophepTagged.prop[0].size(); i++){          
            TophepMatch.push_back_momenta(*top);
            if( delR(*top, TophepTagged.GetLorentzVector(i) ) < match_dist ){
              TophepMatched.push_back_momenta( TophepTagged.GetLorentzVector(i) );          
              iTagged_hep++;     
            }                   
          }
        
        } //match  
        
        if(iTagged_kin != 0)    
         histTagged_kin->Fill(iTagged_kin);
        if(iTagged_chi != 0)    
         histTagged_chi->Fill(iTagged_chi);
        if(iTagged_hep != 0)    
         histTagged_hep->Fill(iTagged_hep);
                       
        particles.clear();
        jets.clear();
        fatjets.clear(); 
        cout<<"Processing data......"<<100.0*float(iEvent+1)/float(nEvent+1)<<" % "<<"\r"; 	
      } //loop over all events 
      cout<<endl;

      Top.make_properties();
      
  
      Topkin.make_properties();
      Topchi.make_properties();
      Tophep.make_properties();

      TopkinTagged.make_properties();
      TopchiTagged.make_properties();
      TophepTagged.make_properties();
      
      Wkin.make_properties();
      Wchi.make_properties();
      Whep.make_properties();
            
      Bkin.make_properties();      
      Bchi.make_properties(); 
      Bhep.make_properties();
          
      if(match){
      
        TopkinMatched.make_properties();
        TopchiMatched.make_properties();
        TophepMatched.make_properties();

        TopkinMatch.make_properties();
        TopchiMatch.make_properties();
        TophepMatch.make_properties();
      
      }    
      
      TCanvas *canvas, *canvas2;
	    sprintf(name,"nJetsCanvas");
	    canvas=  new TCanvas(name,name,1800,1000);
	    sprintf(name,"nTaggedJetsCanvas");
	    canvas2= new TCanvas(name,name,1800,1000);

	    canvas->Divide(2,1);
	    canvas->cd(1);
	    njets20->Draw("hist");
	    canvas->cd(2);
	    njets200->Draw("hist");

      canvas2->cd();
      canvas2->Divide(3,1);
      canvas2->cd(1);
      histTagged_kin->Draw("hist");
      canvas2->cd(2);
      histTagged_chi->Draw("hist");
      canvas2->cd(3);
      histTagged_hep->Draw("hist");

	    sprintf(name, "./%s/jets/nJets.png", FolderName().c_str());
	    canvas->SaveAs(name);
      sprintf(name, "./%s/jets/nTaggedJets.png", FolderName().c_str());
	    canvas2->SaveAs(name);
	    
      cout<<"############################ Analysis Summary #############################"<<endl;

      cout<<"Number of events analysed : "<<nEvent<<endl;
      cout<<"Number of three jet events : "<<nThreeJetEvent<<endl;
      cout<<"Number of Input Jets : "<<nJet<<endl;
      cout<<"Number of Input fATJets : "<<nFatjet<<endl;	    
      cout<<"Number of Tagged top-jets :"<<endl;
      cout<<"    kin : "<<TopkinTagged.prop[0].size()<<"   ,  recoeff_kin : "<<3.0*float(TopkinTagged.prop[0].size())/nJet<<endl;	    
	    cout<<"    chi : "<<TopchiTagged.prop[0].size()<<"   ,  recoeff_chi : "<<3.0*float(TopchiTagged.prop[0].size())/nJet<<endl;
	    cout<<"    hep : "<<TophepTagged.prop[0].size()<<"   ,  recoeff_hep : "<<float(TophepTagged.prop[0].size())/nFatjet<<endl;
	    
	    if(match){
        cout<<"Number of Matched top-jets :"<<endl;
        cout<<setw(20)<<"    kin : "<<setw(20)<<TopkinMatched.prop[0].size()<<setw(20)<<", recoeff_kin : "<<3.0*float(TopkinMatched.prop[0].size())/nJet<<endl;	    
	      cout<<setw(20)<<"    chi : "<<setw(20)<<TopchiMatched.prop[0].size()<<setw(20)<<", recoeff_chi : "<<3.0*float(TopchiMatched.prop[0].size())/nJet<<endl;
	      cout<<setw(20)<<"    hep : "<<setw(20)<<TophepMatched.prop[0].size()<<setw(20)<<", recoeff_hep : "<<float(TophepMatched.prop[0].size())/nFatjet<<endl;
	    }
	
	    delete canvas;
	    delete canvas2;
	    delete njets20;
	    delete njets200;
	    delete histTagged_kin;
	    delete histTagged_chi;
	    delete histTagged_hep;
	    
	    
	    return 0;    
    } //jetreco  
    
  














































