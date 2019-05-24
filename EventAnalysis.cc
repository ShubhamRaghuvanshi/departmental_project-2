
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
#include "TLegend.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TLegendEntry.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "./HEPTopTagger.hh"
#include "Analyser.cc"

using namespace std;

using namespace fastjet;

  const int nprop=10;
  const float coneradius = 0.5;  
  int z_mass[3] = {500,1000,2000};

    
int main(){
 
  char filename[50],treename[50],canvasname[50];

  TCanvas *c_z[nprop];
  TLegend *leg_z[nprop];  
  THStack *zstack[nprop];    
  TH1F *hist_z[nprop]; 
  
  TCanvas *c_top[nprop];
  TLegend *leg_top[nprop], *leg_topbar[nprop];  
  THStack *topstack[nprop], *topbarstack[nprop];    
  TH1F *hist_top[nprop], *hist_topbar[nprop] ;

  TCanvas *c_w[nprop];
  TLegend *leg_wplus[nprop], *leg_wminus[nprop];  
  THStack *wplusstack[nprop], *wminusstack[nprop];    
  TH1F *hist_wplus[nprop], *hist_wminus[nprop] ;

  TCanvas *c_b[nprop];
  TLegend *leg_b[nprop], *leg_bbar[nprop];  
  THStack *bstack[nprop], *bbarstack[nprop];    
  TH1F *hist_b[nprop], *hist_bbar[nprop] ;

  TCanvas *c_quark[nprop];
  TLegend *leg_quark[nprop];  
  THStack *quarkstack[nprop];    
  TH1F *hist_quark[nprop];

  TCanvas *c_antiquark[nprop];
  TLegend *leg_antiquark[nprop];  
  THStack *antiquarkstack[nprop];    
  TH1F *hist_antiquark[nprop];


  TCanvas *c_hadron[nprop];
  TLegend *leg_hadron[nprop];  
  THStack *hadronstack[nprop];    
  TH1F *hist_hadron[nprop];

  TCanvas *c_hadronictop[nprop];
  TLegend *leg_hadronictop[nprop];  
  THStack *hadronictopstack[nprop];    
  TH1F *hist_hadronictop[nprop];

  TCanvas *c_hadronicb[nprop];
  TLegend *leg_hadronicb[nprop];  
  THStack *hadronicbstack[nprop];    
  TH1F *hist_hadronicb[nprop];

  TCanvas *c_hadronicw[nprop];
  TLegend *leg_hadronicw[nprop];  
  THStack *hadronicwstack[nprop];    
  TH1F *hist_hadronicw[nprop];

  TCanvas *c_topchi[nprop];
  TLegend *leg_topchi[nprop];  
  THStack *topchistack[nprop];    
  TH1F *hist_topchi[nprop];

  TCanvas *c_wchi[nprop];
  TLegend *leg_wchi[nprop];  
  THStack *wchistack[nprop];    
  TH1F *hist_wchi[nprop];

  TCanvas *c_bchi[nprop];
  TLegend *leg_bchi[nprop];  
  THStack *bchistack[nprop];    
  TH1F *hist_bchi[nprop]; 

  TCanvas *c_tophep[nprop];
  TLegend *leg_tophep[nprop];  
  THStack *tophepstack[nprop];    
  TH1F *hist_tophep[nprop];

  TCanvas *c_whep[nprop];
  TLegend *leg_whep[nprop];  
  THStack *whepstack[nprop];    
  TH1F *hist_whep[nprop];

  TCanvas *c_bhep[nprop];
  TLegend *leg_bhep[nprop];  
  THStack *bhepstack[nprop];    
  TH1F *hist_bhep[nprop];


  TCanvas *c_topkin[nprop];
  TLegend *leg_topkin[nprop];      
  THStack *topkinstack[nprop];    
  TH1F *hist_topkin[nprop];

  TCanvas *c_wkin[nprop];
  TLegend *leg_wkin[nprop];  
  THStack *wkinstack[nprop];    
  TH1F *hist_wkin[nprop];

  TCanvas *c_bkin[nprop];
  TLegend *leg_bkin[nprop];  
  THStack *bkinstack[nprop];    
  TH1F *hist_bkin[nprop];

  TCanvas *c_dtt;
  TLegend *leg_dtt[4];  
  THStack *dttstack[4];    
  TH1F *hist_dtt[4];

  TCanvas *c_dwb;
  TLegend *leg_dwb[4];  
  THStack *dwbstack[4];    
  TH1F *hist_dwb[4];

  TCanvas *c_dwbchi;
  TLegend *leg_dwbchi[4];  
  THStack *dwbchistack[4];    
  TH1F *hist_dwbchi[4];

  TCanvas *c_dttrchi;
  TLegend *leg_dttrchi[4];  
  THStack *dttrchistack[4];    
  TH1F *hist_dttrchi[4];

  TCanvas *c_dwwrchi;
  TLegend *leg_dwwrchi[4];  
  THStack *dwwrchistack[4];    
  TH1F *hist_dwwrchi[4];

  TCanvas *c_dbbrchi;
  TLegend *leg_dbbrchi[4];  
  THStack *dbbrchistack[4];    
  TH1F *hist_dbbrchi[4];

  TCanvas *c_dttrkin;
  TLegend *leg_dttrkin[4];  
  THStack *dttrkinstack[4];    
  TH1F *hist_dttrkin[4];

  TCanvas *c_dwwrkin;
  TLegend *leg_dwwrkin[4];  
  THStack *dwwrkinstack[4];    
  TH1F *hist_dwwrkin[4];

  TCanvas *c_dbbrkin;
  TLegend *leg_dbbrkin[4];  
  THStack *dbbrkinstack[4];    
  TH1F *hist_dbbrkin[4];
//
  TCanvas *c_dttrhep;
  TLegend *leg_dttrhep[4];  
  THStack *dttrhepstack[4];    
  TH1F *hist_dttrhep[4];

  TCanvas *c_dwwrhep;
  TLegend *leg_dwwrhep[4];  
  THStack *dwwrhepstack[4];    
  TH1F *hist_dwwrhep[4];

  TCanvas *c_dbbrhep;
  TLegend *leg_dbbrhep[4];  
  THStack *dbbrhepstack[4];    
  TH1F *hist_dbbrhep[4];


//
  TCanvas *c_dwbkin;
  TLegend *leg_dwbkin[4];  
  THStack *dwbkinstack[4];    
  TH1F *hist_dwbkin[4];

  TCanvas *c_dqq;
  TLegend *leg_dqq[4];  
  THStack *dqqstack[4];    
  TH1F *hist_dqq[4];

  TCanvas *c_dwbh;
  TLegend *leg_dwbh[4];  
  THStack *dwbhstack[4];    
  TH1F *hist_dwbh[4];

  TCanvas *c_njets;
  TLegend *leg_njets[2];      
  TH1F *hist_njets20, *hist_njets150;
 
  gStyle->SetOptTitle(0);

    c_dtt = new TCanvas("dtt_canvas", "dtt_canvas",1800,1000);
    c_dwb = new TCanvas("dwb_canvas", "dwb_canvas",1800,1000);
    c_dwbchi = new TCanvas("dwbchi_canvas", "dwbchi_canvas",1800,1000);

    c_dttrchi = new TCanvas("dttrchi_canvas", "dttrchi_canvas",1800,1000);
    c_dwwrchi = new TCanvas("dwwrchi_canvas", "dwwrchi_canvas",1800,1000);
    c_dbbrchi = new TCanvas("dbbrchi_canvas", "dbbrchi_canvas",1800,1000);

    c_dttrkin = new TCanvas("dttrkin_canvas", "dttrkin_canvas",1800,1000);
    c_dwwrkin = new TCanvas("dwwrkin_canvas", "dwwrkin_canvas",1800,1000);
    c_dbbrkin = new TCanvas("dbbrkin_canvas", "dbbrkin_canvas",1800,1000);

    c_dttrhep = new TCanvas("dttrhep_canvas", "dttrhep_canvas",1800,1000);
    c_dwwrhep = new TCanvas("dwwrhep_canvas", "dwwrhep_canvas",1800,1000);
    c_dbbrhep = new TCanvas("dbbrhep_canvas", "dbbrhep_canvas",1800,1000);


    c_dwbkin = new TCanvas("dwbkin_canvas", "dwbkin_canvas",1800,1000);
    c_dwbh = new TCanvas("dwbh_canvas", "dwb_canvas",1800,1000);
    c_dqq = new TCanvas("dqq_canvas", "dqq_canvas",1800,1000);
  
    float x1=0.57, y1 = 0.35, x2=0.9, y2=0.9;
    for(int iprop=0; iprop<4; iprop++ ){
      sprintf(canvasname,"dtt canvas_%d", iprop);  

      dttstack[iprop] = new THStack(canvasname,canvasname);
      leg_dtt[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dwb canvas_%d", iprop);   
      dwbstack[iprop] = new THStack(canvasname,canvasname);
      leg_dwb[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dwbchi canvas_%d", iprop);   
      dwbchistack[iprop] = new THStack(canvasname,canvasname);
      leg_dwbchi[iprop] = new TLegend(x1,y1,x2,y2);
//
      sprintf(canvasname,"dttrchi canvas_%d", iprop);  
      dttrchistack[iprop] = new THStack(canvasname,canvasname);
      leg_dttrchi[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dwwrchi canvas_%d", iprop);  
      dwwrchistack[iprop] = new THStack(canvasname,canvasname);
      leg_dwwrchi[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dbbrchi canvas_%d", iprop);  
      dbbrchistack[iprop] = new THStack(canvasname,canvasname);
      leg_dbbrchi[iprop] = new TLegend(x1,y1,x2,y2);
//

      sprintf(canvasname,"dttrkin canvas_%d", iprop);  
      dttrkinstack[iprop] = new THStack(canvasname,canvasname);
      leg_dttrkin[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dwwrkin canvas_%d", iprop);  
      dwwrkinstack[iprop] = new THStack(canvasname,canvasname);
      leg_dwwrkin[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dbbrkin canvas_%d", iprop);  
      dbbrkinstack[iprop] = new THStack(canvasname,canvasname);
      leg_dbbrkin[iprop] = new TLegend(x1,y1,x2,y2);
//

      sprintf(canvasname,"dttrhep canvas_%d", iprop);  
      dttrhepstack[iprop] = new THStack(canvasname,canvasname);
      leg_dttrhep[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dwwrhep canvas_%d", iprop);  
      dwwrhepstack[iprop] = new THStack(canvasname,canvasname);
      leg_dwwrhep[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dbbrhep canvas_%d", iprop);  
      dbbrhepstack[iprop] = new THStack(canvasname,canvasname);
      leg_dbbrhep[iprop] = new TLegend(x1,y1,x2,y2);
//

      sprintf(canvasname,"dwbchi canvas_%d", iprop);   
      dwbkinstack[iprop] = new THStack(canvasname,canvasname);
      leg_dwbkin[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dwbh canvas_%d", iprop);   
      dwbhstack[iprop] = new THStack(canvasname,canvasname);
      leg_dwbh[iprop] = new TLegend(x1,y1,x2,y2);

      sprintf(canvasname,"dqq canvas_%d", iprop);   
      dqqstack[iprop] = new THStack(canvasname,canvasname);
      leg_dqq[iprop] = new TLegend(x1,y1,x2,y2);

    }
    
x1=0.65; y1 = 0.6; x2=0.9; y2=0.9;

float xx1=0.57, yy1 = 0.6, xx2=0.9, yy2=0.9;
    for(int iprop =0; iprop<nprop; iprop++){

     sprintf(canvasname,"z canvas_%d", iprop);  
     c_z[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_z[iprop] = new TLegend(x1,y1,x2,y2);
     zstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"top canvas_%d", iprop);  
     c_top[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_top[iprop] = new TLegend(xx1,yy1,xx2,yy2);
     leg_topbar[iprop] = new TLegend(xx1,yy1,xx2,yy2); 
     topstack[iprop] = new THStack(canvasname,canvasname);
     topbarstack[iprop] = new THStack(canvasname,canvasname);       

     sprintf(canvasname,"w canvas_%d", iprop);  
     c_w[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_wplus[iprop] = new TLegend(xx1,yy1,xx2,yy2);
     leg_wminus[iprop] = new TLegend(xx1,yy1,xx2,yy2);
     wplusstack[iprop] = new THStack(canvasname,canvasname);
     wminusstack[iprop] = new THStack(canvasname,canvasname);       

     sprintf(canvasname,"b canvas_%d", iprop);  
     c_b[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_b[iprop] = new TLegend(xx1,yy1,xx2,yy2);
     leg_bbar[iprop] = new TLegend(xx1,yy1,xx2,yy2); 
     bstack[iprop] = new THStack(canvasname,canvasname);
     bbarstack[iprop] = new THStack(canvasname,canvasname);       

     sprintf(canvasname,"quark canvas_%d", iprop);  
     c_quark[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_quark[iprop] = new TLegend(xx1,yy1,xx2,yy2);
     quarkstack[iprop] = new THStack(canvasname,canvasname); 
   
     sprintf(canvasname,"antiquark canvas_%d", iprop);  
     c_antiquark[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_antiquark[iprop] = new TLegend(xx1,yy1,xx2,yy2);
     antiquarkstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"hadron canvas_%d", iprop);  
     c_hadron[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_hadron[iprop] = new TLegend(x1,y1,x2,y2);
     hadronstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"hadronictop canvas_%d", iprop);  
     c_hadronictop[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_hadronictop[iprop] = new TLegend(x1,y1,x2,y2);
     hadronictopstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"hadronicb canvas_%d", iprop);  
     c_hadronicb[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_hadronicb[iprop] = new TLegend(x1,y1,x2,y2);
     hadronicbstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"hadronicw canvas_%d", iprop);  
     c_hadronicw[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_hadronicw[iprop] = new TLegend(x1,y1,x2,y2);
     hadronicwstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"topchi canvas_%d", iprop);  
     c_topchi[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_topchi[iprop] = new TLegend(x1,y1,x2,y2);
     topchistack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"wchi canvas_%d", iprop);  
     c_wchi[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_wchi[iprop] = new TLegend(x1,y1,x2,y2);
     wchistack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"bchi canvas_%d", iprop);  
     c_bchi[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_bchi[iprop] = new TLegend(x1,y1,x2,y2);
     bchistack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"tophep canvas_%d", iprop);  
     c_tophep[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_tophep[iprop] = new TLegend(x1,y1,x2,y2);
     tophepstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"whep canvas_%d", iprop);  
     c_whep[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_whep[iprop] = new TLegend(x1,y1,x2,y2);
     whepstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"bhep canvas_%d", iprop);  
     c_bhep[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_bhep[iprop] = new TLegend(x1,y1,x2,y2);
     bhepstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"topkin canvas_%d", iprop);  
     c_topkin[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_topkin[iprop] = new TLegend(x1,y1,x2,y2);
     topkinstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"wkin canvas_%d", iprop);  
     c_wkin[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_wkin[iprop] = new TLegend(x1,y1,x2,y2);
     wkinstack[iprop] = new THStack(canvasname,canvasname);

     sprintf(canvasname,"bkin canvas_%d", iprop);  
     c_bkin[iprop] =  new TCanvas(canvasname,canvasname,1800,1000);
     leg_bkin[iprop] = new TLegend(x1,y1,x2,y2);
     bkinstack[iprop] = new THStack(canvasname,canvasname);

  }


    c_njets = new TCanvas("njets","njets",1800,1000);
    c_njets->Divide(2,3);
    leg_njets[0]  = new TLegend(x1,y1,x2,y2);
    leg_njets[1]  = new TLegend(x1,y1,x2,y2);

  for(int i=0;i<3;i++){ 

    cout<<"######## Z mass ="<< z_mass[i]<< "############"<<endl;
    sprintf(filename, "../Zp2ttbar_1l_Event%d.root", z_mass[i]);
    sprintf(treename, "z2ttbar_1l_4j");     
   // Analyser analyser(filename, treename, 0, 100, true ); 
    Analyser analyser(filename, treename, 0, 5000, true);
    analyser.RunAnalyser(coneradius,i);
                    
//prop loop 
    for(int iprop =0; iprop<nprop; iprop++){
      hist_z[iprop]      = analyser.Z.HistProp(iprop,0);
      hist_top[iprop] = analyser.Top.HistProp(iprop,1);
     hist_topbar[iprop] = analyser.Topbar.HistProp(iprop,1);  
      hist_wplus[iprop]  = analyser.Wplus.HistProp(iprop,2);
      hist_wminus[iprop] = analyser.Wminus.HistProp(iprop,2);
      hist_b[iprop] = analyser.B.HistProp(iprop,3);
      hist_bbar[iprop] = analyser.Bbar.HistProp(iprop,3);
      hist_quark[iprop]      = analyser.Quark.HistProp(iprop,6);
      hist_antiquark[iprop]      = analyser.AntiQuark.HistProp(iprop,6);
    //  hist_hadron[iprop]      = analyser.Final_State_Hardons.HistProp(iprop,5);

      hist_hadronictop[iprop]      = analyser.Hadronic_Top.HistProp(iprop,1);
      hist_hadronicb[iprop]      = analyser.Hadronic_B.HistProp(iprop,3);
      hist_hadronicw[iprop]      = analyser.Hadronic_W.HistProp(iprop,2);

      hist_topchi[iprop] = analyser.Top_Chi.HistProp(iprop,4);
      hist_wchi[iprop] = analyser.W_Chi.HistProp(iprop,2);
      hist_bchi[iprop] = analyser.B_Chi.HistProp(iprop,3);

      hist_tophep[iprop] = analyser.Top_Hep.HistProp(iprop,4);
      hist_whep[iprop] = analyser.W_Hep.HistProp(iprop,2);
      hist_bhep[iprop] = analyser.B_Hep.HistProp(iprop,3);


      hist_topkin[iprop] = analyser.Top_Kin.HistProp(iprop,4);
      hist_wkin[iprop] = analyser.W_Kin.HistProp(iprop,2);
      hist_bkin[iprop] = analyser.B_Kin.HistProp(iprop,3); 

//cout<<setw(20)<<iprop<<endl;
 
    } //prop loop
    hist_properties(c_z, hist_z, zstack, leg_z, i, coneradius, 0 ,0);
    hist_properties(c_top, hist_top,   hist_topbar, topstack,   topbarstack, leg_top, leg_topbar, i, coneradius );  
    hist_properties(c_w,   hist_wplus, hist_wminus, wplusstack, wminusstack, leg_wplus, leg_wminus, i, coneradius );
    hist_properties(c_b, hist_b,   hist_bbar, bstack,   bbarstack, leg_b, leg_bbar, i, coneradius );     
  
  //  hist_properties(c_hadron, hist_hadron, hadronstack, leg_hadron, i, coneradius,0 );
    hist_properties(c_quark, hist_quark, quarkstack, leg_quark, i, coneradius,0,0 );
    hist_properties(c_antiquark, hist_antiquark, antiquarkstack, leg_antiquark, i, coneradius,0,0 );
    hist_properties(c_hadronictop, hist_hadronictop, hadronictopstack, leg_hadronictop, i, coneradius,0,0 );
    hist_properties(c_hadronicw, hist_hadronicw, hadronicwstack, leg_hadronicw, i, coneradius,0,0 );
    hist_properties(c_hadronicb, hist_hadronicb, hadronicbstack, leg_hadronicb, i, coneradius,0 ,0);

    hist_properties(c_topchi, hist_topchi, topchistack, leg_topchi, i, coneradius, hist_z[0]->GetEntries(), analyser.n_tagged_chi );
    hist_properties(c_wchi, hist_wchi, wchistack, leg_wchi, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_chi  );
    hist_properties(c_bchi, hist_bchi, bchistack, leg_bchi, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_chi  );

    hist_properties(c_tophep, hist_tophep, tophepstack, leg_tophep, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_hep );
    hist_properties(c_whep, hist_whep, whepstack, leg_whep, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_hep  );
    hist_properties(c_bhep, hist_bhep, bhepstack, leg_bhep, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_hep );

    hist_properties(c_topkin, hist_topkin, topkinstack, leg_topkin, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_kin );
    hist_properties(c_wkin, hist_wkin, wkinstack, leg_wkin, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_kin  );
    hist_properties(c_bkin, hist_bkin, bkinstack, leg_bkin, i, coneradius,  hist_z[0]->GetEntries(), analyser.n_tagged_kin  );
 
    HistAngularDistance(analyser.Top ,analyser.Topbar, c_dtt, hist_dtt, dttstack,leg_dtt, i, coneradius );
    HistAngularDistance(analyser.Hadronic_W ,analyser.Hadronic_B, c_dwbh, hist_dwbh, dwbhstack,leg_dwbh, i, coneradius );
    HistAngularDistance(analyser.Wplus ,analyser.B, c_dwb, hist_dwb, dwbstack, leg_dwb, i, coneradius );
    HistAngularDistance(analyser.W_Chi ,analyser.B_Chi, c_dwbchi, hist_dwbchi, dwbchistack, leg_dwbchi, i, coneradius );
    HistAngularDistance(analyser.W_Kin ,analyser.B_Kin, c_dwbkin, hist_dwbkin, dwbkinstack, leg_dwbkin, i, coneradius );
    HistAngularDistance(analyser.Quark ,analyser.AntiQuark, c_dqq, hist_dqq, dqqstack,leg_dqq, i, coneradius );


    HistAngularDistance(analyser.Top_Chireco_tagged ,analyser.Top_Chireco, c_dttrchi, hist_dttrchi, dttrchistack,leg_dttrchi, i, coneradius );
    HistAngularDistance(analyser.W_Chireco_tagged ,analyser.W_Chireco, c_dwwrchi, hist_dwwrchi, dwwrchistack, leg_dwwrchi, i, coneradius );
    HistAngularDistance(analyser.B_Chireco_tagged ,analyser.B_Chireco, c_dbbrchi, hist_dbbrchi, dbbrchistack, leg_dbbrchi, i, coneradius );

    HistAngularDistance(analyser.Top_Kinreco_tagged ,analyser.Top_Kinreco, c_dttrkin, hist_dttrkin, dttrkinstack,leg_dttrkin, i, coneradius );
    HistAngularDistance(analyser.W_Kinreco_tagged ,analyser.W_Kinreco, c_dwwrkin, hist_dwwrkin, dwwrkinstack, leg_dwwrkin, i, coneradius );
    HistAngularDistance(analyser.B_Kinreco ,analyser.B_Kinreco_tagged, c_dbbrkin, hist_dbbrkin, dbbrkinstack, leg_dbbrkin, i, coneradius );

    HistAngularDistance(analyser.Top_Hepreco_tagged ,analyser.Top_Hepreco, c_dttrhep, hist_dttrhep, dttrhepstack,leg_dttrhep, i, coneradius );
    HistAngularDistance(analyser.W_Hepreco_tagged ,analyser.W_Hepreco, c_dwwrhep, hist_dwwrhep, dwwrhepstack, leg_dwwrhep, i, coneradius );
    HistAngularDistance(analyser.B_Hepreco_tagged ,analyser.B_Hepreco, c_dbbrhep, hist_dbbrhep, dbbrhepstack, leg_dbbrhep, i, coneradius );
  
  }

  return 0;
} 






