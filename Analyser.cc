
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
    
    
    int Analyser::AKJets(float R, vector<ParticleProperty> *QCDjets){
    
      TFile *f = new TFile(filename.c_str(), "READ");    
      TTree *tree;
            
      f->GetObject("hadrons",tree);          
                        
      vector<double>  *px=0, *py=0, *pz=0, *e=0;
      TBranch  *bpx=0, *bpy=0, *bpz=0, *be=0;
    
      tree->SetBranchAddress("px",&px,&bpx);
      tree->SetBranchAddress("py",&py,&bpy);
      tree->SetBranchAddress("pz",&pz,&bpz);
      tree->SetBranchAddress("e",&e ,&be);   

      Long64_t tentry  = tree->LoadTree(20);
      bpx->GetEntry(tentry);       bpy->GetEntry(tentry);       bpz->GetEntry(tentry);   be->GetEntry(tentry);
           
      char name[100];
      sprintf(name, "nJets200");  
      TH1F *njets200 = new TH1F(name, name, 8 , 0, 8);
      njets200->GetXaxis()->SetTitle("N_{jets}(p_{T} > 200 GeV)");
      njets200->SetTitle("");
      FormatHist(njets200, 1, 2, false );
     
           
      vector<PseudoJet> particles, jets;    
      ParticleProperty pp_jet, pp_leading_jet;     
      int nEvent  = bpx->GetEntries()-1;
   //  nEvent= 50000;

      cout<<"Entries in file : "<<bpx->GetEntries()<<endl;
      cout<<"Events to read : "<<nEvent<<endl;
        
      for(int iEvent =0; iEvent < nEvent ; iEvent++){
        bpx->GetEntry(iEvent);       bpy->GetEntry(iEvent);       bpz->GetEntry(iEvent);   be->GetEntry(iEvent);
               
        if(px->size() != py->size() || px->size() != pz->size() || px->size() != e->size()  ) {cout<<"ERR 666"<<endl;return -666;}
      
        for(int i=0; i< px->size() ; i++)
		    	particles.push_back(  PseudoJet(px->at(i), py->at(i), pz->at(i), e->at(i))  );
		  
		    JetDefinition jet_def(cambridge_algorithm, R);  		                 
        ClusterSequence cs_jet(particles, jet_def);
                          
        jets = sorted_by_pt(cs_jet.inclusive_jets(200)); 
        
        if(jets.size() !=0 ){
        njets200->Fill(jets.size());
          pp_leading_jet.push_back_momenta(jets[0].px(), jets[0].py(), jets[0].pz(), jets[0].e() );
          for(int ijet=0; ijet< jets.size(); ijet++ )
            pp_jet.push_back_momenta(jets[ijet].px(), jets[ijet].py(), jets[ijet].pz(), jets[ijet].e() );
          
          jets.clear();  
          
        } 
         
        particles.clear();
        cout<<"Processing data......"<<100.0*float(iEvent+1)/float(nEvent+1)<<" % "<<"\r";
        } // events
        cout<<endl;
        pp_jet.make_properties();
        pp_leading_jet.make_properties();

        QCDjets->push_back(pp_jet);
        QCDjets->push_back(pp_leading_jet);
 
// cout<<"200 jet n : "<<njets200->GetEntries()<<endl;
     
        TCanvas *canvas;
	      sprintf(name,"nJetsCanvas");
	      canvas=  new TCanvas(name,name,1800,1000);

	      canvas->cd();
	      njets200->Draw("hist");
       
        sprintf(name, "./%s/jets/nJets.png", FolderName().c_str());
	      canvas->SaveAs(name);
     
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
 cout<<"match = "<<match<<endl;   
    
      TFile f(filename.c_str(), "READ");    
      TTree *tree, *tree2;
      f.GetObject("hadrons", tree);
      
            
      TBranch  *bpx=0, *bpy=0, *bpz=0, *be=0, *btop=0;
      vector<double>  *px=0, *py=0, *pz=0, *e=0;
      TLorentzVector *top=0;
     
      tree->SetBranchAddress("px",&px,&bpx);
      tree->SetBranchAddress("py",&py,&bpy);
      tree->SetBranchAddress("pz",&pz,&bpz);
      tree->SetBranchAddress("e",&e ,&be);   

      if(match){
        
        f.GetObject("partons", tree2);    
        tree2->SetBranchAddress("parton[1]", &top, &btop);
        Long64_t tentry2  = tree2->LoadTree(0);
        btop->GetEntry(tentry2);
      }
      
      Long64_t tentry  = tree->LoadTree(0);
      bpx->GetEntry(tentry);       bpy->GetEntry(tentry);       bpz->GetEntry(tentry);   be->GetEntry(tentry);

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
      histTagged_kin->GetXaxis()->SetTitle("nTagged(kin)");
      FormatHist(histTagged_kin, 1, 2, false );

      sprintf(name, "nTaggedchi");      
      TH1F *histTagged_chi = new TH1F(name, name, 7 , 0, 7);
      histTagged_chi->GetXaxis()->SetTitle("nTagged(#chi^{2})");
      FormatHist(histTagged_chi, 1, 2, false );

      sprintf(name, "nTaggedhep");      
      TH1F *histTagged_hep = new TH1F(name, name, 7 , 0, 7);
      histTagged_hep->GetXaxis()->SetTitle("nTagged(HepTT)");
      FormatHist(histTagged_hep, 1, 2, false );


      vector<PseudoJet> particles, jets, fatjets, djets; 
      int nEvent = bpx->GetEntries()-1;
   // nEvent = 10000;
          
      cout<<"TOTAL EVENTS IN FILE : "<<bpx->GetEntries()-1<<endl;
      cout<<"TOTAL EVENTS TO READ : "<<nEvent<<endl;
          
      if(nEvent > bpx->GetEntries()){
        cout<<"More than you can handle. "<<endl;
        return -777;
      }
                
              
      int nJet=0, nFatjet=0, nThreeJetEvent=0, nFatJetEvent=0; 
       int iTagged_kin=0, iTagged_chi=0, iTagged_hep=0;
      for(int iEvent =0; iEvent < nEvent ; iEvent++){

        

        bpx->GetEntry(iEvent);       bpy->GetEntry(iEvent);       bpz->GetEntry(iEvent);   be->GetEntry(iEvent); 
        
        if(match){
          btop->GetEntry(iEvent);        
          
        }
        
        if(px->size() != py->size() || px->size() != pz->size() || px->size() != e->size()  ) return -666;
      
        for(int i=0; i< px->size() ; i++)
		    	particles.push_back(  PseudoJet(px->at(i), py->at(i), pz->at(i), e->at(i))  );
		  
		    JetDefinition jet_def(antikt_algorithm, R);  
		    JetDefinition hep_def(cambridge_algorithm, fatR);
		                 
        ClusterSequence cs_jet(particles, jet_def);
        ClusterSequence cs_hep(particles, hep_def);
        particles.clear();
                  
       // jets = UseSoftDrop(jets, R);
                  
        jets    = sorted_by_pt(cs_jet.inclusive_jets(20)); 
        fatjets = sorted_by_pt(cs_hep.inclusive_jets(200));

        if(jets.size() > 2){
          LeadingJet.push_back_momenta( TLorentzVector(jets[0].px(), jets[0].py(), jets[0].pz(),jets[0].e()  )  );
          
          for(int i=0; i< jets.size(); i++){        
            if(abs ( jets[i].eta()  ) < 2.5 )
              djets.push_back(jets[i]);
          }
       }
        
       if(djets.size() > 2){
                
        nThreeJetEvent++;
        nJet+=djets.size();
        njets20->Fill(djets.size() );
                          
        int lastsize_kin = Topkin.prop[0].size();
        int lastsize_chi = Topchi.prop[0].size();
                 
        iTagged_kin=0; iTagged_chi=0;  
                       
        iTagged_kin =topjetreco_kin(djets,   &Topkin, &Wkin, &Bkin, &Bkin) - lastsize_kin;
        iTagged_chi =topjetreco_chi(djets,   &Topchi, &Wchi, &Bchi, &Bchi) - lastsize_chi;
           
           
                  
        if(match ){      
             
          Top.push_back_momenta(*top);

          if(iTagged_kin >0) {         
            if( delR(*top, Topkin.GetLorentzVector( Topkin.prop[0].size() -1 ) ) < match_dist ){
              TopkinMatch.push_back_momenta(*top);  
              TopkinMatched.push_back_momenta( Topkin.GetLorentzVector(Topkin.prop[0].size() -1) ) ;
            }
          }

          if(iTagged_chi >0) {         
            if( delR(*top, Topchi.GetLorentzVector( Topchi.prop[0].size() -1 ) ) < match_dist ){
              TopchiMatch.push_back_momenta(*top);  
              TopchiMatched.push_back_momenta( Topchi.GetLorentzVector(Topchi.prop[0].size() -1 ));
            }
          }

       }   
          
        if(iTagged_kin > 0)    
         histTagged_kin->Fill(iTagged_kin);
        if(iTagged_chi > 0)    
         histTagged_chi->Fill(iTagged_chi);
                              
        jets.clear();
        djets.clear(); 
       } //jets 
               
       if(fatjets.size() > 0){
          njets200->Fill(fatjets.size());  
          nFatjet+=fatjets.size();
          nFatJetEvent++;          
          iTagged_hep=0;
          int lastsize_hep = Tophep.prop[0].size();
          
          iTagged_hep =topjetreco_hep(&fatjets, &Tophep, &Whep, &Bhep, &Bhep) - lastsize_hep;

          if(!match) {
            for(int i=0; i< fatjets.size(); i++){
              FatTop.push_back_momenta(fatjets[i].px(), fatjets[i].py(), fatjets[i].pz(), fatjets[i].e());
              //cout<<fatjets[i].pt()<<endl;
            }  
          } 

          if(match){

            FatTop.push_back_momenta(*top);

            if(iTagged_hep>0){
            
              if( delR(*top, Tophep.GetLorentzVector(Tophep.prop[0].size()-1) ) < match_dist ){
                TophepMatch.push_back_momenta(*top);
                TophepMatched.push_back_momenta( Tophep.GetLorentzVector(Tophep.prop[0].size()-1) );
                
              }
            
            }
            
          } //match

          if(iTagged_hep > 0) 
            histTagged_hep->Fill(iTagged_hep);

          fatjets.clear(); 
          
        }//fatjets

//cout<<"end of event : "<<iEvent<<setw(20)<<iTagged_kin<<setw(20)<<iTagged_chi<<endl;
        cout<<"Processing data......"<<100.0*float(iEvent+1)/float(nEvent+1)<<" % "<<"\r"; 	
      } //loop over all events 
      cout<<endl;
        
      LeadingJet.make_properties(); 
  
      Topkin.make_properties();
      Topchi.make_properties();
      Tophep.make_properties();
      
      Wkin.make_properties();
      Wchi.make_properties();
      Whep.make_properties();
            
      Bkin.make_properties();      
      Bchi.make_properties(); 
      Bhep.make_properties();
        
      FatTop.make_properties();  
      if(match){

        Top.make_properties();            
          
            
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

      if(!match){
        cout<<"Number of Input Jets : "<<nJet<<setw(20)<<nJet/3<<endl;
        cout<<"Number of Input fATJets : "<<nFatjet<<setw(20)<<FatTop.prop[0].size()<<endl;	
      
        cout<<"    kin : "<<Topkin.prop[0].size()<<"   , mis recoeff_kin : "<<3.0*float(Topkin.prop[0].size())/float(nJet)<<endl;	    
  	    cout<<"    chi : "<<Topchi.prop[0].size()<<"   , mis recoeff_chi : "<<3.0*float(Topchi.prop[0].size())/float(nJet)<<endl;
  	    cout<<"    hep : "<<Tophep.prop[0].size()<<"   , mis recoeff_hep : "<<float(Tophep.prop[0].size())/float(nFatjet)<<endl;
  
      }
    
	    if(match){
	    
	      cout<<"Number of input tops : "<<Top.prop[0].size()<<endl;
        cout<<"Number of input fat tops : "<<FatTop.prop[0].size()<<endl;
      
        cout<<"Number of Tagged top-jets :"<<endl;
        cout<<"    kin : "<<Topkin.prop[0].size()<<"   ,  tageff_kin : "<<float(Topkin.prop[0].size())/float(Top.prop[0].size())<<endl;	    
	      cout<<"    chi : "<<Topchi.prop[0].size()<<"   ,  tageff_chi : "<<float(Topchi.prop[0].size())/float(Top.prop[0].size())<<endl;
        cout<<"    hep : "<<Tophep.prop[0].size()<<"   ,  tageff_hep : "<<float(Tophep.prop[0].size())/float(FatTop.prop[0].size())<<endl;

        cout<<"Number of Matched top-jets :"<<endl;
    
        cout<<"Number of fattops : "<<FatTop.prop[0].size()<<endl;    
cout<<setw(20)<<"    kin : "<<setw(20)<<TopkinMatched.prop[0].size()<<setw(20)<<", recoeff_kin : "<<float(TopkinMatched.prop[0].size())/float(Top.prop[0].size())<<endl;	    
	  cout<<setw(20)<<"    chi : "<<setw(20)<<TopchiMatched.prop[0].size()<<setw(20)<<", recoeff_chi : "<<float(TopchiMatched.prop[0].size())/float(Top.prop[0].size())<<endl;
 cout<<setw(20)<<"    hep : "<<setw(20)<<TophepMatched.prop[0].size()<<setw(20)<<", recoeff_hep : "<<float(TophepMatched.prop[0].size())/float(FatTop.prop[0].size())<<endl;
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
    
  int Analyser::Reco4Tops(float R, float fatR){
    
      TFile f(filename.c_str(), "READ");    
      TTree *tree, *tree2;
      f.GetObject("hadrons", tree);
                 
      TBranch  *bpx=0, *bpy=0, *bpz=0, *be=0; 
      vector<double>  *px=0, *py=0, *pz=0, *e=0;
     
      tree->SetBranchAddress("px",&px,&bpx);
      tree->SetBranchAddress("py",&py,&bpy);
      tree->SetBranchAddress("pz",&pz,&bpz);
      tree->SetBranchAddress("e",&e ,&be);   
      
      Long64_t tentry  = tree->LoadTree(0);
      bpx->GetEntry(tentry);       bpy->GetEntry(tentry);       bpz->GetEntry(tentry);   be->GetEntry(tentry);

      char name[100];
      sprintf(name, "nJets20");  
      TH1F *njets20 = new TH1F(name, name, 20 , 0, 20);
      njets20->GetXaxis()->SetTitle("N_{jets}(p_{T} > 20 GeV)");
      FormatHist(njets20, 1, 2, false );

      sprintf(name, " nJets200");  
      TH1F *njets200 = new TH1F(name, name, 7 , 0, 7);
      njets200->GetXaxis()->SetTitle("N_{jets}(p_{T} > 200 GeV)");
      FormatHist(njets200, 1, 2, false );

      sprintf(name, "nTaggedchi");      
      TH1F *histTagged_chi = new TH1F(name, name, 7 , 0, 7);
      histTagged_chi->GetXaxis()->SetTitle("nTagged(#chi^{2})");
      FormatHist(histTagged_chi, 1, 2, false );

      sprintf(name, "nTaggedhep");      
      TH1F *histTagged_hep = new TH1F(name, name, 7 , 0, 7);
      histTagged_hep->GetXaxis()->SetTitle("nTagged(HepTT)");
      FormatHist(histTagged_hep, 1, 2, false );

      vector<PseudoJet> particles, jets, fatjets, djets; 

      int nEvent = bpx->GetEntries()-1;
     nEvent = 1000;
          
      cout<<"TOTAL EVENTS IN FILE : "<<bpx->GetEntries()-1<<endl;
      cout<<"TOTAL EVENTS TO READ : "<<nEvent<<endl;
          
      if(nEvent > bpx->GetEntries()){
        cout<<"More than you can handle. "<<endl;
        return -777;
      }
                
              
      int nJet=0, nFatjet=0, nThreeJetEvent=0, nFatJetEvent=0;
      int n2tag=0, n4toptag; 
      int iTagged_kin=0, iTagged_chi=0, iTagged_hep=0;
      int lastsize_chi, lastsize_hep; 
      TLorentzVector v_temp;

      for(int iEvent =0; iEvent < nEvent ; iEvent++){

        bpx->GetEntry(iEvent);       bpy->GetEntry(iEvent);       bpz->GetEntry(iEvent);   be->GetEntry(iEvent); 
                
        if(px->size() != py->size() || px->size() != pz->size() || px->size() != e->size()  ) return -666;
      
        for(int i=0; i< px->size() ; i++)
		    	particles.push_back(  PseudoJet(px->at(i), py->at(i), pz->at(i), e->at(i))  );
		  
		     
		    JetDefinition hep_def(cambridge_algorithm, fatR);
		                 
        
        ClusterSequence cs_hep(particles, hep_def);
        cout<<particles.size()<<endl;
        particles.clear();
                  
       // jets = UseSoftDrop(jets, R);
                  
         
        fatjets = sorted_by_pt(cs_hep.inclusive_jets( 300 ));

        njets200->Fill(fatjets.size());
       
       if(fatjets.size() > 1){
                   
          nFatjet+=fatjets.size();
          nFatJetEvent++;          
          iTagged_hep=0;
          lastsize_hep = Tophep.prop[0].size();
cout<<"hep tagged : "<<fatjets.size()<<endl;          
          iTagged_hep =topjetreco_hep(&fatjets, &Tophep, &Whep, &Bhep, &Bhep) - lastsize_hep;
cout<<"hep tagged : "<<fatjets.size()<<endl;                              
          if(iTagged_hep != 0){ 
            histTagged_hep->Fill(iTagged_hep);

          }
          if( iTagged_hep >=2){
                    
            top1.push_back_momenta(Tophep.GetLorentzVector(lastsize_hep) );
            top2.push_back_momenta(Tophep.GetLorentzVector(lastsize_hep+1) );

            v_temp = Tophep.GetLorentzVector(lastsize_hep) + Tophep.GetLorentzVector(lastsize_hep + 1);
            
            
            if( abs(m_Zp - v_temp.M()) <= delta_mzp ){
            
              cout<<"two tagged #################################################################"<<endl;

              Zphep.push_back_momenta(v_temp );
              
              for(int i=0; i<fatjets.size(); i++){
                for(int j=0; j<fatjets[i].constituents().size(); j++){
                  particles.push_back(PseudoJet(fatjets[i].constituents().at(i).four_mom()) );  
                }
              } 

              JetDefinition jet_def(antikt_algorithm, R); 
              ClusterSequence cs_jet(particles, jet_def);
              jets    = sorted_by_pt(cs_jet.inclusive_jets(20));

              for(int i=0; i< jets.size(); i++){        
              if(abs ( jets[i].eta()  ) < 2.5)
                djets.push_back(jets[i]);
              }
              nJet+=djets.size();
cout<<"two tagged #################################################################"<<setw(20)<<djets.size()<<setw(20)<<particles.size()<<endl;
             njets20->Fill(djets.size() );

              if(djets.size()>5){
                lastsize_chi = Topchi.prop[0].size();
                iTagged_chi=0;
                                        
                iTagged_chi =topjetreco_chi(djets,   &Topchi, &Wchi, &Bchi, &Bchi) - lastsize_chi;                       
                if(iTagged_chi != 0) 
                  histTagged_chi->Fill(iTagged_chi);

                if(iTagged_chi >=2 ){
                                  
                  top3.push_back_momenta(Topchi.GetLorentzVector(lastsize_chi));
                  top4.push_back_momenta(Topchi.GetLorentzVector(lastsize_chi+1));  

                } //2 slow tagged
                              
              } //5 tags                                                                             
            } //z mass tagged

          } // 2 top tagged jets 
        } //two fatjet event
        
        fatjets.clear(); 
        jets.clear(); 
        djets.clear(); 
                  
cout<<"end of event : "<<iEvent<<endl;
        //cout<<"Processing data......"<<100.0*float(iEvent+1)/float(nEvent+1)<<" % "<<"\r"; 	
      } //loop over all events 

      Zphep.make_properties();
      top1.make_properties();
      top2.make_properties();
      top3.make_properties();
      top4.make_properties();

      cout<<endl;
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
      canvas2->Divide(2,1);
      canvas2->cd(1);
      histTagged_chi->Draw("hist");
      canvas2->cd(2);
      histTagged_hep->Draw("hist");

	    sprintf(name, "./%s/jets/nJets.png", FolderName().c_str());
	    canvas->SaveAs(name);
      sprintf(name, "./%s/jets/nTaggedJets.png", FolderName().c_str());
	    canvas2->SaveAs(name);
	    
      cout<<"############################ Analysis Summary #############################"<<endl;

      cout<<"Number of events analysed : "<<nEvent<<endl;

      cout<<"Number of Input Jets : "<<nJet<<setw(20)<<nJet/3<<endl;
      cout<<"Number of Input fATJets : "<<nFatjet<<endl;	    
    
      cout<<"Number of events with atleast two fat jet with pt >" << m_Zp/2 - 300 <<" :  "<<nFatJetEvent<<endl;  
      
      cout<<"Number of events with atleast two top tagged fat jet : "<<top1.prop[0].size()<<endl;
      cout<<"Two Top tag fraction : "<<float(top1.prop[0].size())/float(nFatJetEvent)<<endl;
      
      cout<<"Number of input jets for chi : "<<nJet<<setw(20)<<nJet/3<<endl;
      cout<<"Two top tagging fraction chi : "<<3.0*float(top3.prop[0].size())/float(nJet)<<setw(20)<<top3.prop[0].size()<<endl;
      
      cout<<"Number of events with zp tagged  : "<<Zphep.prop[0].size()<<endl;
      cout<<"Overall Zp tagging fraction : "<<float(Zphep.prop[0].size())/float(nEvent)<<endl;
      
      cout<<"Number of reconstructed events : "<<top3.prop[0].size()<<endl;          
      cout<<" Reconstruction efficiency :"<<float(top3.prop[0].size())/float(nEvent)<<endl;
      	    
	    return 0;    
    } //jetreco  
  














































