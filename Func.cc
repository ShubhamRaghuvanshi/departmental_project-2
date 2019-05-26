
#include "Func.h"
#include "particleproperty.h"

using namespace std;

using namespace fastjet;


const double m_w = 80.0;
const double m_top =172.0;
const double sigma_w =17.7;
const double sigma_top = 14.335;


  double delR( TLorentzVector v1, TLorentzVector v2 ){

   double Yttbar = v1.Rapidity() - v2.Rapidity(); 
   double Phittbar = v1.Phi() - v2.Phi() ;

    while(Phittbar> M_PI) Phittbar -= 2.0*M_PI;
    while(Phittbar <= -M_PI) Phittbar += 2.0*M_PI;
    return sqrt( Yttbar*Yttbar + Phittbar*Phittbar );    
  }

  //returns 1 if top is tagged else 0  
  int topjetreco_chi( vector<PseudoJet> jets, ParticleProperty &T, ParticleProperty &W, ParticleProperty &B ){

    //chi square minimization
    int w1, w2, b, nTagged=0; 
    double chisq_min, chisq_ijk; 
    PseudoJet t_temp, w_temp;
    
    if(jets.size() < 3) {cout<<"vector<pseudojet> cannot be considered to be top candidate"<<endl; return -666;}  
     
    while(jets.size()>2){
      chisq_min=9999999, chisq_ijk=0
      for(int i=0;i<jets.size();i++){
        for(int j=i+1; j<jets.size(); j++){  
          for(int k=j+1; k<jets.size(); k++){
      
            chisq_ijk = pow( (m_top - (jets[i]+jets[j]+jets[k]).m()),2)/sigma_top + pow( (m_w - (jets[i]+jets[j]).m()),2)/sigma_w ;  
            if(chisq_ijk < chisq_min){
              chisq_min = chisq_ijk;  
              w1 = i; w2 = j; b=k;
            }                                      
          }
        }
      }

      t_temp = jets[w1] + jets[w2] + jets[b];
      w_temp = jets[w1] + jets[w2];
      
     T->push_back_momenta( t_temp.px(), t_temp.py(), t_temp.pz(), t_temp.e()  ); 
     W->push_back_momenta( w_temp.px(), w_temp.py(), w_temp.pz(), w_temp.e()  );  
     B->push_back_momenta(b.px(), b.py(), b.pz(), b.e()  );
     
     if( abs( t_temp.m() - m_top ) < 25.0 )
       nTagged++;
     jets.erase(jets.begin() + w1 );
     jets.erase(jets.begin() + w2 );
     jets.erase(jets.begin() + b );     
  
    }
    return nTagged;
  }

  //returns 1 if top is tagged else 0
  int topjetreco_kin( vector<PseudoJet> jets, ParticleProperty &T, ParticleProperty &W, ParticleProperty &B ){

    //kinematic mass cut minimization
    PseudoJet t_temp, w_temp;
    double  dm; 
    int w1, w2, b, nTagged=0;           

    while(jets.size()>2){
  
      dm=9999999
      for(int i=0; i<jets.size(); i++){
        for(int j=i+1 ; j<jets.size() ; j++){
          if ( dm >  abs((jets[i]+jets[j]).m() - m_w )  ){
            dm =  abs((jets[i]+jets[j]).m() - m_w ) ;
            w1 = i;      w2 = j;    
          }   
        }
      }  //w jet loop
    
      dm=9999999;       
      for(int k=0; k<jets.size(); k++){
        if( k != w1 && k != w2 ){
          if( dm > abs((jets[k] + jets[w1] + jets[w2]).m() - m_top)   ){     
            dm = abs((jets[k] + jets[index_wjet1] + jets[index_wjet2]).m() - m_top) ;    
            b = k;                 
          }
        }           
      } //top jet loop           
 
      t_temp = jets[w1] + jets[w2] + jets[b];
      w_temp = jets[w1] + jets[w2];
      
      T->push_back_momenta( t_temp.px(), t_temp.py(), t_temp.pz(), t_temp.e()  ); 
      W->push_back_momenta( w_temp.px(), w_temp.py(), w_temp.pz(), w_temp.e()  );  
      B->push_back_momenta(b.px(), b.py(), b.pz(), b.e()  );
     
      if( abs( t_temp.m() - m_top ) < 25.0 )
       nTagged++;
      jets.erase(jets.begin() + w1 );
      jets.erase(jets.begin() + w2 );
      jets.erase(jets.begin() + b );     
  
    }
    return nTagged;
  }

  int topjetreco_hepTT( vector<PseudoJet> jets, ParticleProperty &Tjet, ParticleProperty &Wjet, ParticleProperty &Bjet ){
  
    int nTagged=0;
  
    for(unsigned ijet=0; ijet<jets.size(); ijet++){      
    
     	HEPTopTagger::HEPTopTagger tagger(jets[ijet]);

      // Unclustering, Filtering & Subjet Settings
      tagger.set_max_subjet_mass(30.);
      tagger.set_mass_drop_threshold(0.8);
      tagger.set_filtering_R(0.3);
      tagger.set_filtering_n(5);
      tagger.set_filtering_minpt_subjet(30.); 
      tagger.set_mode(HEPTopTagger::TWO_STEP_FILTER);  
      tagger.set_top_minpt(200); 
      tagger.set_top_mass_range(120, 240.); 
      tagger.set_fw(0.15); 

      // Run the tagger
      tagger.run();      
			if (tagger.is_tagged()){
			  Tjet->push_back_momenta( tagger.t().px(), tagger.t().py(), tagger.t().pz(), tagger.t().e()  );
				Wjet->push_back_momenta( tagger.W().px(), tagger.W().py(), tagger.W().pz(), tagger.W().e()  );
				Bjet->push_back_momenta( tagger.B().px(), tagger.B().py(), tagger.B().pz(), tagger.B().e()  );
			  if( abs(tagger.t().m() -m_top) < 25.0  )
			    nTagged++;
			}
		}	 //hep top tagger 		
    return nTagged;
  }


  TH1F* HistdelR( vector<TLorentzVector> &v1, vector<TLorentzVector> &v2, string name1, string name2  ){
    
    TLorentzVector p1,p2;  
    char name[50];
    sprintf(name,"delR_%s%s",name1.c_str(), name2.c_str());  
    TH1F *hist  = new TH1F(name,name, 200, 0, 15 );

    for(int i=0; i<v1.size(); i++){
      p1 = v1[i];
      p2 = v2[i];
      hist->Fill(delR(p1,p2));
    }
    return hist;
  }

  void FormatHist(TH1F *hist, int linecolor, int linewidth, bool stat){

    float scale = 1.0/(hist->Integral());
    hist->SetStats(stat);
    hist->SetLineColor(linecolor);  
    hist->SetLineWidth(linewidth);
    hist->Scale(scale);
  }




void hist_properties(TCanvas *c[], TH1F *h1[], TH1F *h2[], THStack *stack1[], THStack *stack2[], TLegend *leg1[], TLegend *leg2[], int mass_index, float coneradius ){
  
  int z_mass[3] = {500,1000,2000};
  char headername[50], canvasname[50];
  char outfoldername[200], outfilename[200];
  TLegendEntry *header[n_prop];

  for(int iprop =0; iprop<n_prop; iprop++){
    FormatHist(h1[iprop], mass_index+1, 2, false);
    stack1[iprop]->Add( h1[iprop] );
    leg1[iprop]->SetTextSize(0.03);
    sprintf(canvasname, "mean[mZ=%d]=%0.2f", z_mass[mass_index], h1[iprop]->GetMean());
    leg1[iprop]->AddEntry(h1[iprop],canvasname,"l" );
    sprintf(canvasname, "Entries[mZ=%d]=%0.0f", z_mass[mass_index], h1[iprop]->GetEntries()-1);
    leg1[iprop]->AddEntry(h1[iprop],canvasname,"l" );


    FormatHist(h2[iprop],mass_index+1,2,false);
    stack2[iprop]->Add( h2[iprop] );
    sprintf(canvasname, "mean[mZ=%d]=%0.2f", z_mass[mass_index], h2[iprop]->GetMean());      
    leg2[iprop]->SetTextSize(0.03);
    leg2[iprop]->AddEntry(h2[iprop],canvasname,"l" );
    sprintf(canvasname, "Entries[mZ=%d]=%0.0f", z_mass[mass_index], h2[iprop]->GetEntries()-1);
    leg2[iprop]->AddEntry(h2[iprop],canvasname,"l" );

    if(mass_index==2){
      c[iprop]->Divide(2,1);    
      c[iprop]->cd(1);     
      leg1[iprop]->SetHeader(h1[iprop]->GetName());  
      header[iprop] = (TLegendEntry*)leg1[iprop]->GetListOfPrimitives()->First();
      header[iprop]->SetTextSize(.03);  
      header[iprop]->SetTextAlign(22); 
      stack1[iprop]->Draw("hist nostack"); 
      leg1[iprop]->Draw();
    }
   
    if(mass_index==2){
      c[iprop]->cd(2);         
      leg2[iprop]->SetHeader(h2[iprop]->GetName());  
      header[iprop] = (TLegendEntry*)leg2[iprop]->GetListOfPrimitives()->First();
      header[iprop]->SetTextSize(.03);  
      header[iprop]->SetTextAlign(22); 
      stack2[iprop]->Draw("hist nostack"); 
      leg2[iprop]->Draw();
    }
  } 

  if(mass_index==2){
    sprintf(outfoldername, "./Z2ttbar_EventAnalysis");  
    int event_directory = mkdir(outfoldername, S_IRUSR | S_IWUSR | S_IXUSR);
    sprintf(outfoldername, "%s/AR_%f",outfoldername,coneradius);
    int cone_directory = mkdir(outfoldername, S_IRUSR | S_IWUSR | S_IXUSR);

 //   cout<<setw(50)<<"directory creation status : "<<event_directory<<setw(20)<<cone_directory<<endl;
    for(int iprop=0; iprop<n_prop; iprop++){
      sprintf(outfilename, "%s/%s.png",outfoldername,h1[iprop]->GetName());
      c[iprop]->SaveAs(outfilename);
    }
  }

  for(int iprop=0;iprop<n_prop;iprop++){
    h1[iprop]->SetName("hist");

    h2[iprop]->SetName("hist");  
  }
}


void hist_properties(TCanvas *c[], TH1F *h1[], THStack *stack1[], TLegend *leg1[], int mass_index, float coneradius,int n_hadronictop,int n_tagged  ){
  int z_mass[3] = {500,1000,2000};
  char headername[50], canvasname[50];
  char outfoldername[200], outfilename[200];
  TLegendEntry *header[n_prop];

  for(int iprop =0; iprop<n_prop; iprop++){
//    FormatHist(h1[iprop], mass_index+1, 2, false);
    FormatHist(h1[iprop], 1, 2, false);
    stack1[iprop]->Add( h1[iprop] );

    //sprintf(canvasname, "mean [mZ = %d] = %0.2f", z_mass[mass_index], h1[iprop]->GetMean());
    sprintf(canvasname, "mean = %0.2f",h1[iprop]->GetMean());

    leg1[iprop]->SetTextSize(0.03);
    leg1[iprop]->AddEntry(h1[iprop],canvasname,"l" );

    if(n_hadronictop==0)
 //   sprintf(canvasname, "Entries [mZ = %d] = %0.0f", z_mass[mass_index], h1[iprop]->GetEntries()-1);
    sprintf(canvasname, "Entries = %0.0f", h1[iprop]->GetEntries());
    else
 //   sprintf(canvasname, "#epsilon_{top reco}[mZ = %d] = %0.4f", z_mass[mass_index], float(n_tagged)/float(n_hadronictop)  );
    sprintf(canvasname, "#epsilon_{top reco} = %0.4f", float(n_tagged)/float(n_hadronictop)  );


    leg1[iprop]->AddEntry(h1[iprop],canvasname,"l" );
    sprintf(canvasname, "Entries = %0.0f", h1[iprop]->GetEntries());
    leg1[iprop]->AddEntry(h1[iprop],canvasname,"l" );

    if(mass_index==2){
      c[iprop]->cd();     
      leg1[iprop]->SetHeader(h1[iprop]->GetName());  
      header[iprop] = (TLegendEntry*)leg1[iprop]->GetListOfPrimitives()->First();
      header[iprop]->SetTextSize(.03);  
      header[iprop]->SetTextAlign(22); 
      stack1[iprop]->Draw("hist nostack"); 
      leg1[iprop]->Draw();
    }
  } 

  if(mass_index==2){
 //   sprintf(outfoldername, "./Z2ttbar_EventAnalysis");

    sprintf(outfoldername, "./qcd_jetanalysis");  
    int event_directory = mkdir(outfoldername, S_IRUSR | S_IWUSR | S_IXUSR);
    sprintf(outfoldername, "%s/AR_%f",outfoldername,coneradius);
    int cone_directory = mkdir(outfoldername, S_IRUSR | S_IWUSR | S_IXUSR);

  //  cout<<setw(50)<<"directory creation status : "<<event_directory<<setw(20)<<cone_directory<<endl;
    for(int iprop=0; iprop<n_prop; iprop++){
      sprintf(outfilename, "%s/%s.png",outfoldername,h1[iprop]->GetName());
      c[iprop]->SaveAs(outfilename);
    }
  }
  for(int iprop=0;iprop<n_prop;iprop++){
    h1[iprop]->SetName("hist");
  }

}

     void HistAngularDistance(ParticleProperty part1, ParticleProperty part2, TCanvas *canvas_dist, TH1F *hist_dist[], THStack *stack[], TLegend *leg[], int mass_index, float coneradius ){

      char name[50], outfilename[200], outfoldername[200];
      int z_mass[3] = {500,1000,2000};
     
      sprintf(name,"#Delta R_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[0]  = new TH1F(name,name, 50, 0, 10 );

      sprintf(name,"#Delta #Phi_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[1]  = new TH1F(name,name, 100, -6, 6 );

      sprintf(name,"#Delta y_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[2]  = new TH1F(name,name, 200, -10, 10 );

      sprintf(name,"#Delta #eta_{%s %s}",part1.GetPartName().c_str(), part2.GetPartName().c_str() );  
      hist_dist[3]  = new TH1F(name,name, 200, -10, 10 );

      float delR, delPhi, delRap, delEta;  

      for(int i=0;i<part1.prop[7].size();i++){          

        delRap = part1.prop[7][i] - part2.prop[7][i];
        delEta = part1.prop[6][i] - part2.prop[6][i];
        delPhi = part1.prop[9][i] - part2.prop[9][i];

        while(delPhi > M_PI) delPhi -= 2.0*M_PI;
        while(delPhi <= -M_PI) delPhi += 2.0*M_PI; 

        delR = sqrt( delRap*delRap + delPhi*delPhi );

        hist_dist[0]->Fill( delR ) ;
        hist_dist[1]->Fill( delPhi ) ;   
        hist_dist[2]->Fill( delRap ) ;
        hist_dist[3]->Fill( delEta ) ;

      }

      for(int iprop=0; iprop<4; iprop++){


        sprintf(name, "mean [mZ = %d] = %0.2f", z_mass[mass_index], hist_dist[iprop]->GetMean());
        leg[iprop]->SetTextSize(0.045);
        leg[iprop]->AddEntry(hist_dist[iprop],name,"l" );
        sprintf(name, "Entries [mZ = %d] = %0.0f", z_mass[mass_index], hist_dist[iprop]->GetEntries()-1);
        leg[iprop]->AddEntry(hist_dist[iprop],name,"l" );

      }

      for(int iprop=0; iprop<4; iprop++){
        FormatHist(hist_dist[iprop], mass_index+1, 2, false);
        stack[iprop]->Add(hist_dist[iprop]);    
      }

      if(mass_index==2){
        TLegendEntry *header[n_prop];
        canvas_dist->Divide(2,2);
  
  
        //for(int iprop=0; iprop<4; iprop++)
        for(int iprop=0; iprop<3; iprop++){
        canvas_dist->cd(iprop+1);
        leg[iprop]->SetHeader(hist_dist[iprop]->GetName());
        header[iprop] = (TLegendEntry*)leg[iprop]->GetListOfPrimitives()->First(); 
        header[iprop]->SetTextSize(.05);  
        header[iprop]->SetTextAlign(22);         
        stack[iprop]->Draw("hist nostack");
        leg[iprop]->Draw();
       }
                 
       sprintf(outfoldername, "./Z2ttbar_EventAnalysis");   
       sprintf(outfoldername, "%s/AR_%f",outfoldername,coneradius);
       sprintf(outfilename, "%s/%s.png",outfoldername,canvas_dist->GetName()); 
       canvas_dist->SaveAs(outfilename);  
     }

     for(int iprop=0;iprop<4;iprop++){
      hist_dist[iprop]->SetName("hist");      
     }
    } 




































