
#include "Analyser.h"

using namespace std;
using namespace fastjet;

float top_mass = 172.0;
float w_mass = 80;

class Analyser{

  private:
    //outputs
     //inputs
    string filename;
    bool all;     

  public:
 
    Analyser(){}
    ~Analyser(){}

    Analyser(string file_name){
      filename = file_name;
      char name[100];
      sprintf(name,"%s/partons", filename.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);
      sprintf(name,"%s/jets", filename.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);   
    } 

     void SetFile(string file_name){
      filename = file_name;
      char name[100];
      sprintf(name,"%s/partons", filename.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);
      sprintf(name,"%s/jets", filename.c_str());
      mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR);   
     } 

    int ReadPartons( vector<ParticleProperty> *pp){

      TFile f(filename.c_str(), "READ");    
      TTree* tree;
      f.GetObject("partons",tree);
                
      TBranch *PartonBranch = 0;
      
      char PartonBranchName[20];
      TLorentzVector *PartonVector=0;    
      
      int nPartons = tree->GetNbranches();    
      if(nPartons == 0 || nPartons > 20)  return -1;
     
      for(int part=0; part< nPartons; part++){
      
        ParticleProperty pp_temp;      
        sprintf(PartonBranchName, "parton[%d]", part);
        tree->SetBranchAddress(PartonBranchName,&PartonBranch,&PartonVector);
        
        Long64_t tentry  = tree->LoadTree(0);
        PartonBranch->GetEntry(tentry); 
 
        for(int iEvent=0; iEvent<PartonBranchName->GetEntries()-1; iEvent++ ){
          PartonBranch->GetEntry(iEvent);
          pp_temp.push_back_momenta(*PartonVector);        
        }
        pp_temp.make_properties(); 
        pp->push_back(pp_temp);
      } // for all partons
    
      return 0;  
    } //read partons

    //algorithm  = 0 for CA, 1 for AKT
    //moethod = 0,1,2 for kin chi and  hep 
    int RecoTopJets(float R, int algorithm, int method, ParticleProperty *Top, ParticleProperty *W, ParticleProperty *B ){
    
      TFile f(filename.c_str(), "READ");    
      TTree* tree;
      f.GetObject("hadrons", tree);
      
      TBranch  *bpx=0,*bpy=0,*bpz=0,*be=0;
      vector<double>  *px=0,*py=0,*pz=0,*e=0;
    
      tree->SetBranchAddress("px",&px,&bpx);
      tree->SetBranchAddress("py",&py,&bpy);
      tree->SetBranchAddress("pz",&pz,&bpz);
      tree->SetBranchAddress("e",&e ,&be);   
   
      Long64_t tentry  = tree->LoadTree(0);
      bpx->GetEntry(tentry);       bpy->GetEntry(tentry);       bpz->GetEntry(tentry);   be->GetEntry(tentry);

      vector<PseudoJet> particles, jets, fatjets; 
      int nEvent = bpx->GetEntries()-1;
      nEvent = 100;

      char name[50], name_taggedjet[50];
      sprintf(name, "nFatJets(p_T^{jet}>200)");  
      TH1F *njets200 = new TH1F(name, name, 20 , 0,20);

      sprintf(name, "nJets(p_T^{jet}>20)");  
      TH1F *njets20 = new TH1F(name, name, 20 , 0,20);

      string methodname[3] = {"kin","chi", "hepTT"};

      sprintf(name, "nTaggedTopJets_%s", methodname[method].c_str() );      
      TH1F *histTagged = new TH1F(name, name, 20 , 0,20);

      sprintf(name, "nTopJets_%s", methodname[method].c_str() );      
      Top->setnames(name,"");
      sprintf(name, "nWJets_%s", methodname[method].c_str() );      
      W->setnames(name,"");
      sprintf(name, "nBJets_%s", methodname[method].c_str() );      
      B->setnames(name,"");

      int iTagged, nTagged=0;
          
      for(int iEvent =0; iEvent < nEvent ; iEvent++){
        bpx->GetEntry(iEvent);       bpy->GetEntry(iEvent);       bpz->GetEntry(iEvent);   be->GetEntry(iEvent);
        iTagged =0;
        if(px->size() != py->size() || px->size() != pz->size() || px->size() != e->size()  ) return 666;
      
        for(int i=0; i< px->size() ; i++)
		    	particles.push_back(  PseudoJet(px->at(i), py->at(i), pz->at(i), e->at(i))  );
		  
		    switch (algorithm){
		      case 1: JetDefinition jet_def(cambridge_algorithm, R);
		      case 2: JetDefinition jet_def(antikt_algorithm, R);  
		      default: return -202;
		    }
		                 
        ClusterSequence cs(particles, jet_def);
        jets = sorted_by_pt(cs.inclusive_jets(20));           
        fatjets = sorted_by_pt(cs.inclusive_jets(200));
  
        njets20->Fill(jets.size() );
        njets200->Fill(fatjets.size());
          
        // make it global  
        switch (method){
        
          case 1: iTagged =topjetreco_kin(jets, *Top, *W, *B)
                  histTagged->Fill(iTagged);
          case 2: iTagged =topjetreco_chi(jets, *Top, *W, *B)
                  histTagged->Fill(iTagged);
          case 3: iTagged =topjetreco_hepTT(jets, *Top, *W, *B)
                  histTagged->Fill(iTagged); 
        } 
          
        nTagged += iTagged;  
        particles.clear();
        jets.clear();
        fatjets.clear(); 
        cout<<"Processing data......"<<100.0*float(iEvent+1)/float(n2-n1+1)<<" % "<<"\r"; 
      } //loop over all events 

      Top->make_properties();
      W->make_properties();
      B->make_properties();
          
      
      TCanvas *canvas, canvas2;
	    sprintf(name,"nJetsCanvas");
	
	    canvas=  new TCanvas(name,name,1800,1000);
	    canvas2= new TCanvas(name,name,1800,1000);

	    canvas->Divide(2,1);

	    canvas->cd(1);
	    njets200->SetLineWidth(2);
	    njets200->Draw("hist");
	    canvas->cd(2);
	    n20->Draw("hist");
	    n20->SetLineWidth(2);
	    sprintf(name, "./%s/jets/njets.png", filename.c_str());
	    canvas->SaveAs(name);
	 
	    canvas2->cd();
      sprintf(name, "./%s/jets/%s.png",filename.c_str(), histTagged->GetName() );
	    histTagged->SetLineWidth(2);
	    histTagged->Draw("hist");
	    canvas2->SaveAs(name);
	    
	    return nTagged;    
    } //jetreco  
  
};  //analyser class





























