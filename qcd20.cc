
  #include "particleproperty.h"
  #include "Analyser.h"
  #include "Func.h" 
  #include <time.h>

  using namespace std;
  using namespace fastjet;

  int main(){
    clock_t t;
    //######################################## Initialize Analyser object #####################

    Analyser analyser;
    vector<ParticleProperty> jet;
    TH1F *hist1, *hist2;
    float x1=0.65, y1 = 0.7, x2=0.88, y2=0.88;   
    char filename[100],foldername[200],pngname[200];
    sprintf(filename, "QCD20.root" );
    analyser.SetFile(filename);    

    analyser.RecoJets(0.5, 1.0, false);

    sprintf(foldername, "./%s/jets", analyser.FolderName().c_str());

////

    jet.push_back(analyser.Tophep);
    jet.push_back(analyser.Topkin);
    jet.push_back(analyser.Topchi);
    jet.push_back(analyser.FatTop);
    DrawHistograms(jet, 0 , -1, foldername);	  
    jet.clear();
////

    hist1 = analyser.Tophep.HistProp(5, 1);
    hist2 = analyser.FatTop.HistProp(5, 1);
    hist1->Divide(hist2);
    FormatHist(hist1, 2, 2, true);


    TCanvas *canvas = new TCanvas("canvas", "canvas", 1600,1000);
    TLegend *leg = new TLegend(x1 -0.5,y1,x2-0.5,y2);
    leg->SetBorderSize(0);
     
    double e = double(analyser.Tophep.prop[0].size())/double(analyser.FatTop.prop[0].size());
    sprintf(pngname, "#epsilon_{mistag}^{HEPTT} = %.4f", e);
    leg->AddEntry(hist1, pngname, "l");  
    canvas->cd();
    hist1->Draw("hist");
    leg->Draw("same");
    sprintf(pngname, "%s/misreco_pt.png", foldername);
    canvas->SaveAs(pngname);
////    
    t = clock() - t;
    cout<<"Time of execution (seconds ): "<<t/float(CLOCKS_PER_SEC)<<endl;
    return 0;
  }
  
  
  // /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/SoftDrop.cc 
  // /home/ehep/Downloads/products/fjcontrib-1.041/RecursiveTools/RecursiveSymmetryCutBase.cc
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
