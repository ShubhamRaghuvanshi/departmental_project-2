
#include "Func.h"
#include "particleproperty.h"

using namespace std;

using namespace fastjet;

  double delR( TLorentzVector v1, TLorentzVector v2 ){

   double Yttbar = v1.Rapidity() - v2.Rapidity(); 
   double Phittbar = v1.Phi() - v2.Phi() ;

    while(Phittbar> M_PI) Phittbar -= 2.0*M_PI;
    while(Phittbar <= -M_PI) Phittbar += 2.0*M_PI;
    return sqrt( Yttbar*Yttbar + Phittbar*Phittbar );    
  }

  int topjetreco_kin( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B, ParticleProperty *TaggedTop ){

    PseudoJet t_temp, w_temp;
    double  dm; 
    int w1, w2, b;           

    while(jets.size()>2){
  
      dm=9999999;
      for(int i=0; i<jets.size(); i++){
        for(int j=i+1 ; j<jets.size() ; j++){
          if ( dm >  abs((jets[i]+jets[j]).m() - m_w ) ){
            dm =  abs((jets[i]+jets[j]).m() - m_w ) ;
            w1 = i;      w2 = j;    
          }   
        }
      }  //w jet loop
    
      dm=9999999;       
      for(int k=0; k<jets.size(); k++){
        if( dm > abs((jets[k] + jets[w1] + jets[w2]).m() - m_top) && k != w1 && k != w2   ){     
          dm = abs((jets[k] + jets[w1] + jets[w2]).m() - m_top) ;    
          b = k;                 
        }                 
      } //top jet loop           
 
      t_temp = jets[w1] + jets[w2] + jets[b];
      w_temp = jets[w1] + jets[w2];
      
      T->push_back_momenta( t_temp.px(), t_temp.py(), t_temp.pz(), t_temp.e()  ); 
      W->push_back_momenta( w_temp.px(), w_temp.py(), w_temp.pz(), w_temp.e()  );  
      B->push_back_momenta(jets[b].px(), jets[b].py(), jets[b].pz(), jets[b].e()  );
     
      if( abs( t_temp.m() - m_top ) < delta_mtop && abs( w_temp.m() - m_w) < delta_mw )
        TaggedTop->push_back_momenta( t_temp.px(), t_temp.py(), t_temp.pz(), t_temp.e()  );
        jets.erase(jets.begin() + w1 );
        jets.erase(jets.begin() + w2 );
        jets.erase(jets.begin() + b );       
    }
    return TaggedTop->prop[0].size();
  }

  //returns 1 if top is tagged else 0  
  int topjetreco_chi( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B, ParticleProperty *TaggedTop  ){

    //chi square minimization
    int w1, w2, b; 
    double chisq_min, chisq_ijk; 
    PseudoJet t_temp, w_temp;
    
    if(jets.size() < 3) {cout<<"vector<pseudojet> cannot be considered to be top candidate"<<endl; return -666;}  
     
    while(jets.size()>2){
      chisq_min=9999999; chisq_ijk=0;
      for(int i=0;i<jets.size();i++){
        for(int j=i+1; j<jets.size(); j++){  
          for(int k=0; k<jets.size(); k++){
      
            chisq_ijk = pow( (m_top - (jets[i]+jets[j]+jets[k]).m()),2)/sigma_top + pow( (m_w - (jets[i]+jets[j]).m()),2)/sigma_w ;  
            if(chisq_ijk < chisq_min && k !=i && k!=j){
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
     B->push_back_momenta( jets[b].px(), jets[b].py(), jets[b].pz(), jets[b].e()  );
     
     if( abs( t_temp.m() - m_top ) < delta_mtop && abs( w_temp.m() - m_w) < delta_mw )
        TaggedTop->push_back_momenta( t_temp.px(), t_temp.py(), t_temp.pz(), t_temp.e()  );
        jets.erase(jets.begin() + w1 );
        jets.erase(jets.begin() + w2 );
        jets.erase(jets.begin() + b );       
    }
    return TaggedTop->prop[0].size();
  }


  int topjetreco_hep( vector<PseudoJet> jets, ParticleProperty *Tjet, ParticleProperty *Wjet, ParticleProperty *Bjet, ParticleProperty *TaggedTop  ){
  
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
				Bjet->push_back_momenta( tagger.b().px(), tagger.b().py(), tagger.b().pz(), tagger.b().e()  );
			  if( abs(tagger.t().m() -m_top) < delta_mtop && abs(tagger.W().m() - m_w) < delta_mw  )
			    TaggedTop->push_back_momenta( tagger.t().px(), tagger.t().py(), tagger.t().pz(), tagger.t().e()  );
			}
		}	 //hep top tagger 		
    return TaggedTop->prop[0].size();
  }

// histogram draw methods

  void DrawdelR( ParticleProperty p1, ParticleProperty p2, string foldername ){

    char canvasname[200],pngname[200];
    TH1F *hist;
          
    sprintf(canvasname, "#Delta R_{%s_%s}", p1.GetPartName().c_str(), p2.GetPartName().c_str() );
    hist  = new TH1F(canvasname, canvasname, 50, 0, 6);
    hist->GetXaxis()->SetTitle(canvasname);
    FormatHist(hist, 2, 3, false);
    
    TLorentzVector v1,v2;  
    
    for(int i=0; i<p1.prop[7].size(); i++){
      v1 = p1.GetLorentzVector(i);
      v2 = p2.GetLorentzVector(i);
      hist->Fill( delR(v1, v2));    
    }
    
    TCanvas *canvas;  
    sprintf(canvasname, "%s_%s_canvas", p1.GetPartName().c_str(), p2.GetPartName().c_str() );
    canvas = new TCanvas(canvasname, canvasname, 1600, 1000);
    sprintf(pngname,"./%s/delR_%s_%s.png", foldername.c_str(), p1.GetPartName().c_str(), p2.GetPartName().c_str() );
    
    canvas->cd();
    hist->Draw("e1");
    hist->Draw("hist same");
    
    delete hist;
    delete canvas;
  }

  void FormatHist(TH1F *hist, int linecolor, int linewidth, bool stat){
    
    hist->SetLineColor(linecolor);  
    hist->SetLineWidth(linewidth);    
    hist->SetStats(false);    
    
    if(stat){
      float scale = 1.0/(hist->Integral());
      hist->Scale(scale);
      hist->GetYaxis()->SetTitle("Normalised units");      
      //hist->Sumw2();
    }
    else{
      hist->GetYaxis()->SetTitle("Number of events");
                     
    } 
    hist->SetName("");
    hist->SetTitle("");
  }


  int DrawHistograms(vector<ParticleProperty> particle, int drawoption, int mass_index, string foldername){
    char canvasname[200],pngname[200];
              
    int xdiv, ydiv;
    float x1=0.60, y1 = 0.7, x2=0.9, y2=0.9;
    string legentry[3] = {"M_{Z'} = 500 GeV","M_{Z'} = 1 TeV", "M_{Z'} = 2 TeV"};

    if( particle.size() % 2 ==0){
      xdiv = particle.size()/2;
      ydiv = particle.size()/2;      
    }  
    else{
      ydiv = particle.size()/2;
      xdiv = ydiv+1;
    }

    TCanvas *canvas;
    THStack *stack;
    TLegend *legend;    

    for(int iprop=0; iprop<n_prop; iprop++){     
    
      if(drawoption ==0 ) {
        legend =  new TLegend(0.65, 0.8, 0.9, 0.9);
        for(int ipart=0; ipart < particle.size(); ipart++){      
        
          sprintf(canvasname, "%s_%s", particle[ipart].GetPartName().c_str(), particle[ipart].GetPropName(iprop).c_str() );
          canvas = new TCanvas(canvasname, canvasname, 1600, 1000);               
        
          TH1F *hist;           
          hist = particle[ipart].HistProp(iprop, -1); 
          FormatHist(hist, mass_index+1, 4, false);          
          legend->AddEntry(hist, legentry[mass_index].c_str(), "l");    
          canvas->cd();
            
          hist->Draw("e1");
          hist->Draw("hist same");  
          legend->Draw("same");
          sprintf(pngname,"./%s/%s.png", foldername.c_str(), canvasname);          
          canvas->SaveAs(pngname);      
          delete hist;  
          delete canvas; 
          delete legend;         
        }     
      } // 0 : draw on separate canvas

      if(drawoption ==1 ) {
        
        TCanvas *canvas;    
        sprintf(canvasname, "%s_%s_%lu", particle[0].GetPartName().c_str(), particle[0].GetPropName(iprop).c_str(), particle.size() );
        canvas = new TCanvas(canvasname, canvasname, 1600, 1000);
        canvas->Divide(xdiv, ydiv);
        for(int ipart=0; ipart < particle.size(); ipart++){     
           
          TH1F *hist;           
          hist = particle[ipart].HistProp(iprop, -1); 
          FormatHist(hist, 1, 4, false);          
            
          canvas->cd(ipart+1);
            
          hist->Draw("e1");
          hist->Draw("hist same");  
          delete hist;                      
        } 
       sprintf(pngname,"./%s/%s.png", foldername.c_str(), canvasname );         
        canvas->SaveAs(pngname);      
        
        delete canvas;    
      } // 1 : draw on divided canvas

     
      if(drawoption == 2){
      
        if(particle.size() !=2 ){
          cout<<"He's got so much in his heart, but he doesn't know what to do."<<endl;
          return -666;
        }

  
        sprintf(canvasname, "%s_%s_div_%s_%s", particle[0].GetPartName().c_str(), particle[0].GetPropName(iprop).c_str(), 
        particle[1].GetPartName().c_str(), particle[1].GetPropName(iprop).c_str() );
        
        sprintf(pngname,"./%s/%s.png", foldername.c_str(), canvasname );      

        TH1F *hist1, *hist2 ;
        char histname[200];
        hist1 = particle[0].HistProp(iprop, 1); 
        hist2 = particle[1].HistProp(iprop, 1); 
        hist1->Divide(hist2);

        sprintf(histname, "#epsilon_%s", particle[1].GetPartName().c_str() );
        hist1->GetYaxis()->SetTitle(histname);        
        hist1->SetLineWidth(3);
        hist1->SetLineColor(2);
        
        hist1->Draw("e");
        hist1->Draw("hist same");  

        canvas->SaveAs(pngname);      
                    
        delete hist1;
        delete hist2;
        delete canvas;
      } // 2: draw the ratio of two 
     
     
      if(drawoption == 3 ) {
      
        stack  = new THStack("stack","");
        legend =  new TLegend(x1,y1,x2,y2);    

        sprintf(canvasname, "%s_%s", particle[0].GetPartName().c_str(), particle[0].GetPropName(iprop).c_str() );
        canvas = new TCanvas(canvasname, canvasname, 1600, 1000);               
              
        canvas->cd();                  
        for(int ipart=0; ipart < particle.size(); ipart++){      
   
          TH1F *hist;           
          hist = particle[ipart].HistProp(iprop, -1); 
          FormatHist(hist, ipart+1, 4, true);          
            
          hist->Draw("e1 same"); 
          stack->Add(hist);  
          legend->AddEntry(hist, legentry[ipart].c_str(), "l");
          
          if(ipart == particle.size()-1){
            stack->Draw("hist nostack same");            
            stack->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());          
            legend->Draw("same");
            sprintf(pngname,"./%s/%s_stack.png", foldername.c_str(), canvasname);          
            canvas->SaveAs(pngname);      
          }
        }
                                               
        delete stack;
        delete canvas;
        delete legend;
       // delete hist;     
      } //3 : draw on top

     
     
     
      } //prop
  return 0;         
  }
  
  

















