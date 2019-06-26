
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


  double delR( PseudoJet v1, PseudoJet v2 ){

   double Yttbar = v1.rap() - v2.rap(); 
   double Phittbar = v1.phi() - v2.phi() ;

    while(Phittbar> M_PI) Phittbar -= 2.0*M_PI;
    while(Phittbar <= -M_PI) Phittbar += 2.0*M_PI;
    return sqrt( Yttbar*Yttbar + Phittbar*Phittbar );    
  }

  vector<PseudoJet> UseSoftDrop( vector<PseudoJet> jets, double R){
  
    vector<PseudoJet> groomedjets;

    double z_cut = 0.10;
    double beta = 2.0;

    contrib::SoftDrop sd(beta, z_cut, R);
    
    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
      // Run SoftDrop and examine the output
      PseudoJet sd_jet = sd(jets[ijet]);
    
      //  cout << "original    jet: " << jets[ijet] << endl;
      //  cout << "SoftDropped jet: " << sd_jet << endl;
      groomedjets.push_back(sd_jet);
      assert(sd_jet != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
   
//      cout << " masses before and after grooming "<<setw(20)<<jets[ijet].constituents().size()<<setw(20)<<sd_jet.constituents().size()<<endl;
//      cout << " masses before and after grooming "<<setw(20)<<jets[ijet].m()<<setw(20)<<sd_jet.m()<<endl;
    }
 
    return groomedjets;
  }

   float sizeofjet( TLorentzVector parent, vector<double> *h_px, vector<double> *h_py, vector<double> *h_pz, vector<double> *h_e  ){
  
    float R;
    TLorentzVector jetmomenta, hadronmomenta;  
   
    for(R=0; R < 10.0; R = R+0.1 ){
      jetmomenta.SetPxPyPzE(0,0,0,0);
  
      for(int i=0; i< h_px->size();  i++){ 
        hadronmomenta.SetPxPyPzE( h_px->at(i), h_py->at(i), h_pz->at(i), h_e->at(i)  );      

        if(  parent.DeltaR(hadronmomenta) < R ) 
          jetmomenta = jetmomenta + hadronmomenta;          
      }
      if(jetmomenta.M() >= parent.M()) break;
   }//R loop
   return R;
 }


  int topjetreco_kin( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B, ParticleProperty *TopTagged ){

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
      if( abs( t_temp.m() - m_top ) < delta_mtop && abs( w_temp.m() - m_w) < delta_mw ){
      
        
      
        T->push_back_momenta( t_temp.px(), t_temp.py(), t_temp.pz(), t_temp.e()  ); 
        W->push_back_momenta( w_temp.px(), w_temp.py(), w_temp.pz(), w_temp.e()  );  
        B->push_back_momenta(jets[b].px(), jets[b].py(), jets[b].pz(), jets[b].e()  );

       // break;
      }  
        jets.erase(jets.begin() + w1 );
        jets.erase(jets.begin() + w2 );
        jets.erase(jets.begin() + b );       
    }
    return T->prop[0].size();
  }

  //returns 1 if top is tagged else 0  
  int topjetreco_chi( vector<PseudoJet> jets, ParticleProperty *T, ParticleProperty *W, ParticleProperty *B, ParticleProperty *TopTagged){

    //chi square minimization
    int w1, w2, b; 
    double chisq_min, chisq_ijk; 
    PseudoJet t_temp, w_temp;
    vector<TLorentzVector> v_sort;
   //if(jets.size() < 3) {cout<<"vector<pseudojet> cannot be considered to be top candidate : "<<jets.size()<<endl; return -666;}  
     
    while(jets.size()>2){
      chisq_min=9999999; chisq_ijk=0;
      for(int i=0;i<jets.size();i++){
        for(int j=i+1; j<jets.size(); j++){  
          for(int k=0; k<jets.size(); k++){
      
            chisq_ijk = pow( (m_top - (jets[i]+jets[j]+jets[k]).m()),2)/(sigma_top*sigma_top) + pow( (m_w - (jets[i]+jets[j]).m()),2)/(sigma_w*sigma_w) ;  
            if(chisq_ijk < chisq_min && k !=i && k!=j){
              chisq_min = chisq_ijk;  
              w1 = i; w2 = j; b=k;
            }                                      
          }
        }
      }
           
     if( chisq_min < 2 ){

      t_temp = jets[w1] + jets[w2] + jets[b];
      w_temp = jets[w1] + jets[w2];
  
       T->push_back_momenta( t_temp.px(), t_temp.py(), t_temp.pz(), t_temp.e()  ); 
       W->push_back_momenta( w_temp.px(), w_temp.py(), w_temp.pz(), w_temp.e()  );  
       B->push_back_momenta( jets[b].px(), jets[b].py(), jets[b].pz(), jets[b].e()  );
     // break;
     }
    
    jets.erase(jets.begin() + w1 );
    jets.erase(jets.begin() + w2 );
    jets.erase(jets.begin() + b );       
 
    }
    return T->prop[0].size();
  }


  int topjetreco_hep( vector<PseudoJet> *jets, ParticleProperty *Tjet, ParticleProperty *Wjet, ParticleProperty *Bjet, ParticleProperty *TopTagged ){
  
    for(unsigned ijet=0; ijet<jets->size(); ijet++){      
     	HEPTopTagger::HEPTopTagger tagger(jets->at(ijet));

      // Unclustering, Filtering & Subjet Settings
      tagger.set_max_subjet_mass(30.);
      tagger.set_mass_drop_threshold(0.8);
      tagger.set_filtering_R(0.3);
      tagger.set_filtering_n(5);
      tagger.set_filtering_minpt_subjet(30.); 
      tagger.set_mode(HEPTopTagger::TWO_STEP_FILTER);  
      tagger.set_top_minpt(200); 
      tagger.set_top_mass_range(m_top - 2.5*delta_mtop, m_top + 2.5*delta_mtop); 
      tagger.set_fw(0.15); 

      // Run the tagger
      tagger.run();      
			if (tagger.is_tagged()){
			  Tjet->push_back_momenta( tagger.t().px(), tagger.t().py(), tagger.t().pz(), tagger.t().e()  );
				Wjet->push_back_momenta( tagger.W().px(), tagger.W().py(), tagger.W().pz(), tagger.W().e()  );
				Bjet->push_back_momenta( tagger.b().px(), tagger.b().py(), tagger.b().pz(), tagger.b().e()  );	
				jets->erase(jets->begin() + ijet);			      	      
			}
		}	 //hep top tagger 		
    return Tjet->prop[0].size();
  }

// histogram draw methods

  int DrawdelR( vector<ParticleProperty> particle, int mass_index, string foldername ){


    TCanvas *canvas;  
    TLegend *legend;
    TH1F *hist;
    TLine *line;
    THStack *stack;         

    string legentry[3] = {"M_{Z'} = 500 GeV","M_{Z'} = 1000 GeV", "M_{Z'} = 2000 GeV" };
    string algoentry[3] = {"kinematic mass cut", "#chi^{2} minimization", "HepTT"};

    char canvasname[200],pngname[200];


    sprintf(canvasname, "#DeltaR(top_{parton}, top_{tagged} ) for %s", legentry[mass_index].c_str() );
    canvas = new TCanvas(canvasname, canvasname, 1600, 1000);
    canvas->cd();

    //hist->GetXaxis()->SetTitle(canvasname);
    //FormatHist(hist, 2, 3, false);
    
    float x1=0.65, y1 = 0.7, x2=0.88, y2=0.88;  
    legend =  new TLegend(x1,y1,x2,y2);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03); 
    stack  = new THStack("stack","");
    
    TLorentzVector v1,v2;        
    
    for(int ipart =0 ; ipart < particle.size(); ipart = ipart + 2){
    
      hist  = new TH1F(canvasname, canvasname, 50, 0, 4);
      FormatHist(hist, ipart/2 + 2, 3, false);
      hist->SetFillColor(ipart/2+2);
    
      if(particle[ipart].prop[7].size() != particle[ipart+1].prop[7].size()) {
        cout<<"Tumse na ho paega beta, Tumhare lachhan hame bilkul thik nhi lag rhe beta tumse na ho paega. "<<endl;
        cout<<"sizes : "<<particle[ipart].prop[7].size()<<setw(20)<<particle[ipart+1].prop[7].size()<<endl;
        return -666;
      }      
      for(int i=0; i<particle[ipart].prop[7].size(); i++){
        v1 = particle[ipart].GetLorentzVector(i);
        v2 = particle[ipart+1].GetLorentzVector(i);
        hist->Fill( delR(v1, v2));    
      }
      stack->Add(hist);
      legend->AddEntry(hist, algoentry[ipart/2].c_str(), "f");   
    }      
     
    stack->Draw(); 
    stack->GetXaxis()->SetTitle("#DeltaR(top_{parton}, top_{tagged})"); 
    stack->GetYaxis()->SetTitle("Number of events"); 
    sprintf(canvasname,"#DeltaR(top_{parton}, top_{tagged}) for %s", legentry[mass_index].c_str() );
    stack->SetTitle(canvasname); 
    legend->Draw("same");    
    canvas->Update();
    line = new TLine(0.3,0,0.3 , canvas->GetUymax());
    line->SetLineColor(kBlack);
    line->SetLineWidth(1);
    line->Draw();             
    sprintf(pngname,"./%s/delR.png", foldername.c_str());
       
    canvas->SaveAs(pngname);
    delete hist;
    delete canvas;
    delete stack;
    delete legend;
    delete line;

    return 0;
  }

  void FormatHist(TH1F *hist, int linecolor, int linewidth, bool stat){
    
    hist->SetLineColor(linecolor);  
    hist->SetLineWidth(linewidth);    
      
    
    if(stat){
      float scale = 1.0/(hist->Integral());
      hist->Scale(scale);
      hist->GetYaxis()->SetTitle("Normalised units"); 
      hist->SetStats(false) ;    
      //hist->Sumw2();
    }
    else{
      hist->GetYaxis()->SetTitle("Number of events");
      hist->SetStats(false);                     
    } 
    
    hist->GetYaxis()->SetTitleSize(.04); 
    hist->GetXaxis()->SetTitleSize(.04);
   // hist->SetTitleSize(0.02);
    hist->SetTitle("");
  }


  int DrawHistograms(vector<ParticleProperty> particle, int drawoption, int mass_index, string foldername){
    char canvasname[200],pngname[200];
              
    int xdiv, ydiv;
    float x1=0.65, y1 = 0.7, x2=0.88, y2=0.88;
    string legentry[3] = {"M_{Z'}=500","M_{Z'}=1000", "M_{Z'}=2000"};
    string algoentry[4] = {"kinematic mass cut", "#chi^{2} method", "HepTopTagger","Original Top Quark"};


    TCanvas *canvas;
    THStack *stack;
    TLegend *legend;    
    
    //cout<<"particles to be plotted : "<<particle.size()<<endl;

    for(int iprop=0; iprop<10; iprop++){  
    
      if(drawoption ==0 ) {
        
        for(int ipart=0; ipart < particle.size(); ipart++){      
        
          sprintf(canvasname, "%s_%s", particle[ipart].GetPartName().c_str(), particle[ipart].GetPropName(iprop).c_str() );
          canvas = new TCanvas(canvasname, canvasname, 1600, 1000);               
          
          TH1F *hist;           
          hist = particle[ipart].HistProp(iprop, 1); 
               
          FormatHist(hist, 1 , 2, true);          
          
          canvas->cd();

        //  hist->Draw("e1");
          hist->Draw("hist"); 
          
          sprintf(pngname,"./%s/%s.png", foldername.c_str(), canvasname);          
          canvas->SaveAs(pngname);    
                      
          delete hist;  
          delete canvas; 
            
        }     
      } // 0 : draw on separate canvas

      else if(drawoption ==1 ) {
      
        int batchsize =1;
        float recoeff;
        if(particle.size() == 2 ){
          xdiv = 2;
          ydiv = 1;
        }  
        else if(particle.size() == 3 ){
          xdiv = 3;
          ydiv = 1;
        }  
        else if(particle.size() == 4 ){
          xdiv = 2;
          ydiv = 2;
        }  
        
        else if(particle.size() == 12 ||particle.size() == 15){
          xdiv = 2;
          ydiv = 2;                  
        }
        else{ cout<<"define new canvas divisions (xdiv, ydiv)"<<endl;    return -888; }
        
        
        
        sprintf(canvasname,"%s", particle[0].propXaxis[iprop].c_str());
        canvas = new TCanvas(canvasname, canvasname, 1600, 1000);
        canvas->Divide(xdiv, ydiv);
        
        TH1F *hist;                         
        for(int iCanvas=0; iCanvas< xdiv*ydiv; iCanvas++ ){
          canvas->cd(iCanvas+1);
          stack  = new THStack("stack","");  
     //     legend =  new TLegend( 0.59, y1 - 0.03 , 0.9 , y2);
          
     //     legend->SetTextSize(0.045); 
      //    for(int iBatch= batchsize-1; iBatch>=0; iBatch--){
        
          for(int iBatch = 0; iBatch<batchsize; iBatch++){             
           // hist = particle[ (batchsize)*iCanvas + iBatch + batchsize ].HistProp(iprop, 3);

            hist = particle[iCanvas*(batchsize) + iBatch].HistProp(iprop, 1);
            FormatHist(hist, iBatch+1, 2, true);
            stack->Add(hist);
  
       //   recoeff =  float( hist->GetEntries() ) /  float( particle[iBatch +batchsize].prop[0].size()) ;
      //    if(iCanvas == 3)
      //     recoeff =  float( hist->GetEntries() ) /  float( particle[iBatch].prop[0].size()) ;
  

     //     if(mass_index == 0)
     //       sprintf(canvasname, "#epsilon_{tag}(%s)=%.2f " , legentry[iBatch].c_str(), recoeff  );            
     //     else
     //       sprintf(canvasname, "#epsilon_{reco}(%s)=%.2f" , legentry[iBatch].c_str(), recoeff  );            

      
    //      if(iCanvas == 0 )
    //        legend->AddEntry(hist, legentry[iBatch].c_str(), "l");
    //      else       
    //        legend->AddEntry(hist, canvasname, "l"); 
            
     //       hist->Draw("e1");
            hist->SetName("hist");
           // hist->Draw("e1 same"); 
            //if(iBatch == batchsize-1){ 
            if(iBatch == batchsize-1){                            
                          
             // sprintf(pngname, "%s", particle[ (batchsize)*iCanvas + iBatch + batchsize ].GetPartName().c_str());
             sprintf(pngname, "%s", particle[ iCanvas*(batchsize) + iBatch ].GetPartName().c_str());
              stack->SetTitle(pngname);
//              stack->SetTitle(hist->GetTitle());
              stack->Draw("hist nostack");
              sprintf(pngname,"%s^{top}",particle[0].propXaxis[iprop].c_str());
              
              stack->GetXaxis()->SetTitleSize(.04);
              stack->GetYaxis()->SetTitleSize(.04);
              stack->GetXaxis()->SetTitle(pngname);
              sprintf(pngname,"%s", hist->GetYaxis()->GetTitle());
              stack->GetYaxis()->SetTitle(pngname);              
     //         legend->SetBorderSize(0); 

//if(iprop != 4)
       //       legend->Draw("same");            

            }
         }  //batch          
        }  //canvas
        if(mass_index == 0 )
          sprintf(canvasname, "TagParticles%lu_%s", particle.size(), particle[0].GetPropName(iprop).c_str());
        else
          sprintf(canvasname, "RecoParticles%lu_%s", particle.size(), particle[0].GetPropName(iprop).c_str());

        sprintf(pngname,"./%s/%s.png", foldername.c_str(), canvasname );         
        canvas->SaveAs(pngname);      
        delete hist;    
        delete canvas;
    //    delete legend;  
        delete stack;
      } // 1 : draw on divided canvas
      
      else if(drawoption == 2){
      

        TH1F *hist2, *hist1, *hist3;

        hist2 = particle[particle.size()-2].HistProp(iprop, mass_index+1);
        hist3 = particle[particle.size()-1].HistProp(iprop, mass_index+1);
        
 //   cout<<"enties : "<<hist2->GetEntries()<<setw(20)<<hist3->GetEntries()<<endl;
        
        sprintf(canvasname, " %s_%s_DIV_%s_%s", particle[0].GetPartName().c_str(), particle[0].GetPropName(iprop).c_str(), 
        particle[particle.size()-1].GetPartName().c_str(), particle[particle.size()-1].GetPropName(iprop).c_str() );
         
        canvas = new TCanvas(canvasname, canvasname, 1600, 1000);
        canvas->cd();     

        stack  = new THStack("stack","");
        legend =  new TLegend(x1-0.53 ,y1+0.01, x2 -0.53,y2+0.01);    
        
        for(int ipart= 0; ipart< particle.size()-2 ; ipart = ipart+1){
          char histname[200];

          hist1 = particle[ipart].HistProp(iprop, mass_index+1);
          hist1->SetName(""); 
          if(ipart == particle.size() -3)      
            hist1->Divide(hist3);
          else 
            hist1->Divide(hist2);
            
          FormatHist(hist1, ipart+2, 3, false);
         // hist1->SetFillColor(ipart +2);
          stack->Add(hist1);
          legend->AddEntry(hist1, algoentry[ipart].c_str(), "l");
                   
          if(ipart == particle.size() -3){
            stack->Draw("hist nostack");                        

            sprintf(pngname,"top reconstruction efficiency vs %s for %s", particle[particle.size()-2].propXaxis[iprop].c_str(), legentry[mass_index].c_str());
            stack->SetTitle(pngname);
            sprintf(pngname,"%s^{%s}", particle[0].propXaxis[iprop].c_str(), particle[particle.size()-2].GetPartName().c_str() );
            stack->GetXaxis()->SetTitle(pngname);
            sprintf(pngname, "#epsilon_{%sreco}", particle[particle.size()-2].GetPartName().c_str() );
            stack->GetYaxis()->SetTitle(pngname);        

            legend->SetBorderSize(0);
            legend->Draw("same");                  
            sprintf(pngname,"./%s/Recoeff%lu_1.0_%s_%s.png", foldername.c_str(), particle.size(), 
            particle[particle.size()-1].GetPartName().c_str(),particle[particle.size()-2].GetPropName(iprop).c_str() );               
            canvas->SaveAs(pngname);  
            
          }                                              
        }
        
        delete canvas;  
        delete hist1;  
        delete stack;
        delete hist2;
        delete hist3;
        delete legend;  
      } // 2: draw the ratio of two 
     
     
      else if(drawoption == 3 ) {
      
        stack  = new THStack("stack","");
        legend =  new TLegend(x1,y1,x2,y2);    

        sprintf(canvasname, "%s_%s", particle[0].GetPartName().c_str(), particle[0].GetPropName(iprop).c_str() );
        canvas = new TCanvas(canvasname, canvasname, 1600, 1000);               
              
        canvas->cd();                  
        for(int ipart= particle.size()-1; ipart >=0 ; ipart--){      
   
          TH1F *hist;           
          hist = particle[ipart].HistProp(iprop, -1); 
          FormatHist(hist, ipart+1, 4, false);          

          if(ipart == particle.size() -1 ){
            sprintf(pngname,"Reconstructed  %s^{top} for %s ", particle[0].propXaxis[iprop].c_str(), legentry[mass_index].c_str() );
            hist->SetTitle(pngname);
            sprintf(pngname,"%s^{top}", particle[0].propXaxis[iprop].c_str() );
            hist->GetXaxis()->SetTitle(pngname); 
          }
 
 
          hist->Draw("e1 same"); 
          stack->Add(hist);  
          legend->AddEntry(hist, particle[ipart].GetPartName().c_str(), "l");
          
          if(ipart == 0){

            stack->Draw("hist nostack same"); 
                       
//            stack->GetXaxis()->SetTitle(pngname);
//            stack->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle()); 
            legend->SetBorderSize(0);         
            legend->Draw("same");
            sprintf(pngname,"./%s/%s_stack.png", foldername.c_str(), canvasname);          
            canvas->SaveAs(pngname);      
          }
        }
                                               
        delete stack;
        delete canvas;
        delete legend;
       // delete hist;     
      } //3 : draw on same canvas

     else{}
     
     
      } //prop
  return 0;         
  }
  
 
 void Drawone( ParticleProperty p, int iprop , string foldername){

  TCanvas *canvas;
  TLegend *legend;    
  TH1F *hist;

  float x1=0.60, y1 = 0.7, x2=0.85, y2=0.85;

  legend =  new TLegend(x1,y1,x2,y2);    
  
  
  hist = p.HistProp(iprop, -1); 
 
  FormatHist(hist, 1, 4, true);          

  char canvasname[100], legentry[100];

  sprintf(canvasname, "%s_%s", p.GetPartName().c_str(), p.GetPropName(iprop).c_str() );
  canvas = new TCanvas(canvasname, canvasname, 1600, 1000);               
  sprintf(legentry, "%s", hist->GetTitle() );

  legend->AddEntry(hist,legentry, "l");
    
  sprintf(legentry, "./%s/%s_%s.png", foldername.c_str(), p.GetPartName().c_str(), p.GetPropName(iprop).c_str());
  hist->Draw("hist");
 // legend->Draw("same");
  canvas->SaveAs(legentry);
 delete canvas;
 delete legend;
 delete hist;
 
 }
 
 void Draw1v2( ParticleProperty p1, ParticleProperty p2, int iprop1 , int iprop2 ,string foldername){

  char canvasname[100], legentry[100];
  
  TCanvas *canvas;
  TH2F *hist;

  float x1=0.60, y1 = 0.7, x2=0.85, y2=0.85;

  char xname[50], yname[50];
  sprintf(xname, "%s^{%s}", p1.propXaxis[iprop1].c_str(), p1.GetPartName().c_str() );
  sprintf(yname, "%s^{%s}", p2.propXaxis[iprop2].c_str(), p2.GetPartName().c_str() );
  sprintf(canvasname, "%s vs %s", xname, yname );
      
  int propsize;
  
  if(p1.prop[iprop1].size() <  p2.prop[iprop2].size() )
    propsize = p1.prop[iprop1].size();
  else    
    propsize = p2.prop[iprop2].size();

  float min1 = 9999999, max1 = -9999999,min2 = 9999999, max2 = -9999999;
  for(int i=0; i< propsize; i++){
  
    if(min1 > p1.prop[iprop1][i]) min1 = p1.prop[iprop1][i];
    if(max1 < p1.prop[iprop1][i]) max1 = p1.prop[iprop1][i];
    if(min2 > p2.prop[iprop2][i]) min2 = p2.prop[iprop2][i];
    if(max2 < p2.prop[iprop2][i]) max2 = p2.prop[iprop2][i]; 
  }
  
  hist = new TH2F(canvasname, canvasname, 50, 0, 800, 50, 0, 5 );
  
  for(int i=0; i< propsize; i++){
    hist->Fill(p1.prop[iprop1][i], p2.prop[iprop2][i]);  
  }

  hist->GetXaxis()->SetTitle(xname);
  hist->GetYaxis()->SetTitle(yname);
//  hist->SetTitle("");

  sprintf(canvasname, "%s_%s", p1.GetPartName().c_str(), p1.GetPropName(iprop1).c_str() );
  canvas = new TCanvas(canvasname, canvasname, 1600, 1000);               
  
  gStyle->SetPalette(kRainBow);
  canvas->cd();
  sprintf(legentry, "./%s/%s_vs_%s.png", foldername.c_str(), xname, yname);
  hist->SetStats(false);
  hist->Draw("COLZ");
  //legend->Draw("same");
  canvas->SaveAs(legentry);
  
  delete canvas;
  delete hist;
 }

 int DrawdelRvspT(ParticleProperty top, ParticleProperty w, ParticleProperty b, string foldername ){
 
  if(top.prop[5].size() != w.prop[5].size() || top.prop[5].size() != b.prop[5].size()){
    cout<<"Equal sized particles needed"<<endl;
    return -999;
  }
 
  char canvasname[100], legentry[100];
  
  TCanvas *canvas;
  TH2F *hist;

  char xname[50], yname[50];
  sprintf(xname, "%s^{%s}", top.propXaxis[5].c_str(), top.GetPartName().c_str() );
  sprintf(yname, " #DeltaR(%s %s)", w.GetPartName().c_str(), b.GetPartName().c_str() );
  sprintf(canvasname, "%s vs %s", xname, yname );
      
  float min1 = 9999999, max1 = -9999999;
  
  for(int i=0; i< top.prop[5].size(); i++){  
    if(min1 > top.prop[5][i]) min1 = top.prop[5][i];
    if(max1 < top.prop[5][i]) max1 = top.prop[5][i];
  }
  
  hist = new TH2F(canvasname, canvasname, 50, 0, 800, 50, 0, 5 );
  
  for(int i=0; i< top.prop[5].size(); i++){
    hist->Fill(  top.prop[5][i], delR(w.GetLorentzVector(i), b.GetLorentzVector(i) )  );  
  }

  hist->GetXaxis()->SetTitle(xname);
  hist->GetYaxis()->SetTitle(yname);
 //  hist->SetTitle("");

 
  canvas = new TCanvas(canvasname, canvasname, 1600, 1000);                
  gStyle->SetPalette(kRainBow);
  canvas->cd();
  sprintf(legentry, "./%s/%s_vs_%s.png", foldername.c_str(), xname, yname);
  hist->SetStats(false);
  hist->Draw("COLZ");
  //legend->Draw("same");
  canvas->SaveAs(legentry);
  
  delete canvas;
  delete hist;
  
 
 } 
















