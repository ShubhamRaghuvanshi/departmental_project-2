

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>  
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLorentzVector.h"
//  #include "particleproperty.cc"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
//#include "./HEPTopTagger.hh"
//#include "Func.cc"

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
  int topjetreco_chi( vector<PseudoJet> jets, TLorentzVector &top, TLorentzVector &w, TLorentzVector &b ){

        //chi square minimization
         int index_wjet1, index_wjet2, index_bjet; 
         double chisq_min=9999999, chisq_ijk=0; 
          for(int i=0;i<jets.size();i++){
            for(int j=i+1; j<jets.size(); j++){  
              for(int k=0; k<jets.size(); k++){
                if( k!=i && k!=j ){
      
                  chisq_ijk = (m_top - (jets[i]+jets[j]+jets[k]).m())*(m_top - (jets[i]+jets[j]+jets[k]).m() )/sigma_top + (m_w - (jets[i]+jets[j]).m())*(m_w - (jets[i]+jets[j]).m())/sigma_w ;  

//cout<<setw(20)<<"ijk : "<<setw(20)<<i<<setw(20)<<j<<setw(20)<<k<<setw(20)<<chisq_ijk<<setw(20)<<chisq_min<<setw(20)<<(m_top - (jets[i]+jets[j]+jets[k]).m())*(m_top - (jets[i]+jets[j]+jets[k]).m() )/sigma_w<<setw(20)<<(m_w - (jets[i]+jets[j]).m())*(m_w - (jets[i]+jets[j]).m())/sigma_w<<endl;
                  if(chisq_ijk < chisq_min){
                    chisq_min = chisq_ijk;  
                    index_wjet1 = i; index_wjet2 = j; index_bjet=k;
                  }                      
                }
              }
            }
          }
// cout<<"chi : "<<setw(20)<<setw(20)<<index_wjet1<<setw(20)<<index_wjet2<<endl;
 


          double w_px =  (jets[index_wjet1] + jets[index_wjet2]).px() ;
          double w_py =  (jets[index_wjet1] + jets[index_wjet2]).py() ;
          double w_pz =  (jets[index_wjet1] + jets[index_wjet2]).pz() ;
          double w_e  =  (jets[index_wjet1] + jets[index_wjet2]).e() ;
 
          double b_px =  (jets[index_bjet]).px() ;
          double b_py =  (jets[index_bjet]).py() ;
          double b_pz =  (jets[index_bjet]).pz() ;
          double b_e  =  (jets[index_bjet]).e() ;
 
          double top_px = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).px() ;
          double top_py = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).py() ;
          double top_pz = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).pz() ;
          double top_e  = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).e() ;
 
          w.SetPxPyPzE( w_px, w_py, w_pz, w_e );
          b.SetPxPyPzE( b_px, b_py, b_pz, b_e );
          top.SetPxPyPzE( top_px, top_py, top_pz, top_e ); 
          if(abs ( m_top - top.M()) < 25  ) return 1;  
          else return 0;       
  }

  //returns 1 if top is tagged else 0
  int topjetreco_kin( vector<PseudoJet> jets, TLorentzVector &top, TLorentzVector &w, TLorentzVector &b ){

        //kinematic mass cut minimization
         double  dmw_kin=9999999; 
         int index_wjet1, index_wjet2, index_bjet;           
          //search through event for two jets from w
          for(int i=0; i<jets.size(); i++){
            for(int j=i+1 ; j<jets.size() ; j++){
              if(i!=j){    
                if ( dmw_kin >  abs((jets[i]+jets[j]).m() - m_w )  ){
                  dmw_kin =  abs((jets[i]+jets[j]).m() - m_w ) ;
                  index_wjet1 = i;      index_wjet2 = j;    
                }   
               // cout<<setw(20)<<"WJET : "<<setw(20)<<i<<setw(20)<<j<<setw(20)<<dm_wjet<<endl; 
              }
            }
          }  //w jet loop
  
          if( dmw_kin < 99999 ) {
          
             double dmt_kin=9999999; 
            // search from rest of the jets, a jet from top     
            for(int k=0; k<jets.size(); k++){
              if( k != index_wjet1 && k != index_wjet2 ){
                if( dmt_kin > abs((jets[k] + jets[index_wjet1] + jets[index_wjet2]).m() - m_top)   ){     
                  dmt_kin = abs((jets[k] + jets[index_wjet1] + jets[index_wjet2]).m() - m_top) ;    
                  index_bjet = k;                 
                }
               // cout<<setw(20)<<"TOP JET : "<<k<<setw(20)<<index_wjet1<<setw(20)<<index_wjet2<<setw(20)<<dm_topjet<<endl;
              }           
            } //top jet loop           
          }


          double w_px =  (jets[index_wjet1] + jets[index_wjet2]).px() ;
          double w_py =  (jets[index_wjet1] + jets[index_wjet2]).py() ;
          double w_pz =  (jets[index_wjet1] + jets[index_wjet2]).pz() ;
          double w_e  =  (jets[index_wjet1] + jets[index_wjet2]).e() ;
  
          double b_px =  (jets[index_bjet]).px() ;
          double b_py =  (jets[index_bjet]).py() ;
          double b_pz =  (jets[index_bjet]).pz() ;
          double b_e  =  (jets[index_bjet]).e() ;

          double top_px = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).px() ;
          double top_py = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).py() ;
          double top_pz = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).pz() ;
          double top_e  = (jets[index_bjet] + jets[index_wjet1] + jets[index_wjet2]).e() ;
 
          w.SetPxPyPzE( w_px, w_py, w_pz, w_e );
          b.SetPxPyPzE( b_px, b_py, b_pz, b_e );
          top.SetPxPyPzE( top_px, top_py, top_pz, top_e );    
        
          if(abs ( m_top - top.M()) < 25  ) return 1;  
          else return 0;  
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



/*
int main(){

  int z_mass[3] = {500,1000,2000};
  char filename[50],treename[50];

  sprintf(filename, "../Zp2ttbar_1l_Event%d.root", z_mass[0]);
  sprintf(treename, "z2ttbar_1l_4j");  
  
  vector<double> px[2],py[2],pz[2],e[2];
  int status = AnalyseEvents(filename, treename, 0, 99, 0);

  return 0;
} 


*/


































