

#include "particleproperty.h"
#include "Func.h"

using namespace std;


ParticleProperty::ParticleProperty(vector<double> &px, vector<double> &py, vector<double> &pz, vector<double> &e, string name, string append){

  particle_name = name;
  particle_surname = append;

  prop[0] = px; 
  prop[1] = py;
  prop[2] = pz;
  prop[3]  = e; 

}


void ParticleProperty::push_back_momenta( TLorentzVector v){
    prop[0].push_back(v.Px());            
    prop[1].push_back(v.Py());            
    prop[2].push_back(v.Pz());            
    prop[3].push_back(v.E());            
}

void ParticleProperty::push_back_momenta( vector<TLorentzVector> v){

  for(int i=0; i< v.size(); i++){
    prop[0].push_back(v[i].Px());            
    prop[1].push_back(v[i].Py());            
    prop[2].push_back(v[i].Pz());            
    prop[3].push_back(v[i].E());            
  }

}

int ParticleProperty::push_back_momenta(vector<double> px, vector<double> py, vector<double> pz, vector<double> e){

  if(px.size() != py.size() || px.size() != py.size() ||px.size() != pz.size() ||px.size() != e.size())
    return -666;
  for(int i=0; i<px.size(); i++){
    prop[0].push_back(px[i]);            
    prop[1].push_back(py[i]);            
    prop[2].push_back(pz[i]);            
    prop[3].push_back(e[i]);
  }
  return 0;              
}
    
void ParticleProperty::push_back_momenta(double px, double py, double pz, double e){

    prop[0].push_back(px);            
    prop[1].push_back(py);            
    prop[2].push_back(pz);            
    prop[3].push_back(e);
                
}
    
int ParticleProperty::add_momenta(ParticleProperty p1, ParticleProperty p2){
 
  if(p1.prop[0].size() != p2.prop[0].size()){
    cout<<"unequal sizes : "<<p1.prop[0].size()<<setw(20)<<p2.prop[0].size()<<endl;
    return -333;
  }
  for(int i=0; i<p1.prop[0].size(); i++){
    prop[0].push_back(p1.prop[0][i] + p2.prop[0][i]); 
    prop[1].push_back(p1.prop[1][i] + p2.prop[1][i]);
    prop[2].push_back(p1.prop[2][i] + p2.prop[2][i]); 
    prop[3].push_back(p1.prop[3][i] + p2.prop[3][i]);
  
  }
  make_properties();
  return 0;
}    
    
void ParticleProperty::clear_properties(){

  for(int iprop=0; iprop<n_prop; iprop++ )
    prop[iprop].clear();
}    

void ParticleProperty::push_property(int iprop, float val){

  prop[iprop].push_back(val);
}    
    
// constructs all properties given 4 momenta vectors
void ParticleProperty::make_properties(){
  vector<double> px,py,pz,e;
  px = prop[0];
  py = prop[1];
  pz = prop[2];
  e  = prop[3];    

  TLorentzVector v;

  for(int i=0;i<prop[0].size();i++){
    v.SetPxPyPzE( px[i],py[i],pz[i],e[i] );
    prop[4].push_back(v.M());            
    prop[5].push_back(v.Perp());
    prop[6].push_back(v.Eta());
    prop[7].push_back(v.Rapidity());    
    prop[8].push_back(v.Phi());
    prop[9].push_back(v.Theta());
    
  } 
}



void ParticleProperty::setnames( string name, string append ){
  particle_name = name;
  particle_surname = append;
  
}


void ParticleProperty::SetFromFile(string filename, string treename, int treeindex){

  TFile f(filename.c_str(),"READ");
  TTree *tree; 
  f.GetObject(treename.c_str(),tree);
  
  vector<double> *p[4]={0,0,0,0};
  TBranch *b_p[4] = {0,0,0,0};

  char name[50];
  for(int iprop=0; iprop<4;iprop++){
    sprintf(name,"%s_%s%s",particle_name.c_str(), property[iprop].c_str(), particle_surname.c_str());  
    tree->SetBranchAddress(name,&p[iprop], &b_p[iprop]) ;
  }
    
  Long64_t tentry  = tree->LoadTree(treeindex);

  for(int iprop=0; iprop<4;iprop++)
    b_p[iprop]->GetEntry(tentry);

  for(int iprop=0; iprop<4;iprop++)
    b_p[iprop]->GetEntry(0);

  prop[0] = *p[0];
  prop[1] = *p[1];
  prop[2] = *p[2];
  prop[3] = *p[3];

  TLorentzVector v;

  for(int i=0;i<prop[0].size();i++){
    v.SetPxPyPzE( prop[0][i], prop[1][i],  prop[2][i], prop[3][i] );
    prop[4].push_back(v.M());            
    prop[5].push_back(v.Perp());
    prop[6].push_back(v.Eta());
    prop[7].push_back(v.Rapidity());    
    prop[8].push_back(v.Phi());
    prop[9].push_back(v.Theta());
  }
  f.Close();
}
        
void ParticleProperty::PrintVectorSizes(){
  cout<<setw(20)<<particle_name<<setw(20)<<particle_surname;
  for(int i=0;i<n_prop;i++)
    cout<<setw(8)<<prop[i].size();
  cout<<endl;
}

void ParticleProperty::PrintEntry(int n){
  for(int i=0;i<n_prop;i++)
    cout<<setw(15)<<prop[i][n]  ;
  cout<<endl;
} 

 
TH1F* ParticleProperty::HistProp(int iprop, int particle_id){

  double  prop_low[n_prop]  ;
  double  prop_high[n_prop]  ;

    double temp_high = -9999999, temp_low = 9999999; 
    for(int i=0; i < prop[iprop].size(); i++){
      
      if(temp_high <  prop[iprop][i] )     
        temp_high = prop[iprop][i] ;
    
      if(temp_low >  prop[iprop][i] )     
        temp_low = prop[iprop][i] ;   
    }
    prop_low[iprop] = temp_low;
    prop_high[iprop] = temp_high;
    
//cout<<setw(20)<<propXaxis[iprop]<<setw(20)<<prop_low[iprop]<<setw(20)<<prop_high[iprop]<<setw(20)<<prop[0].size()<<endl;
    
  if(particle_id == 0 ){

        prop_low[0] = -2000.0;   prop_high[0] = 2000.0;
        prop_low[1] = -2000.0;   prop_high[1] = 2000.0;
        prop_low[2] = -2500.0;   prop_high[2] = 2500.0;
        prop_low[3] = 0.0;       prop_high[3] = 2500.0;
        prop_low[4] = 400.0;     prop_high[4] = 2100.0;
        prop_low[5] = 0.0;       prop_high[5] = 1500.0;
        prop_low[6] = -6.0;     prop_high[6] = 6.0;
        prop_low[7] = -6.0;     prop_high[7] = 6.0;
        prop_low[8] = -6.0;      prop_high[8] = 6.0;
        prop_low[9] = 0.0;       prop_high[9] = 6.0;         
  }
  if(particle_id == 1 ){

        prop_low[0] = -2000.0;   prop_high[0] = 2000.0;
        prop_low[1] = -2000.0;   prop_high[1] = 2000.0;
        prop_low[2] = -2500.0;   prop_high[2] = 2500.0;
        prop_low[3] = 0.0;       prop_high[3] = 2500.0;
     
       // prop_low[4] = m_top - delta_mtop -5.0;       prop_high[4] =  m_top + delta_mtop + 5.0;
        
        prop_low[4] = 0;         prop_high[4] = 300.0; 
        prop_low[5] = 200.0;       prop_high[5] = 850.0;
        prop_low[6] = -6.0;     prop_high[6] = 6.0;
        prop_low[7] = -6.0;     prop_high[7] = 6.0;
        prop_low[8] = -6.0;      prop_high[8] = 6.0;
        prop_low[9] = 0.0;       prop_high[9] = 6.0;         
  }
  if(particle_id == 2 ){

        prop_low[0] = -2000.0;   prop_high[0] = 2000.0;
        prop_low[1] = -2000.0;   prop_high[1] = 2000.0;
        prop_low[2] = -2500.0;   prop_high[2] = 2500.0;
        prop_low[3] = 0.0;       prop_high[3] = 2500.0;
     
        prop_low[4] = m_top - delta_mtop -2.0;       prop_high[4] =  m_top + delta_mtop + 2.0;
        
        prop_low[5] = 0.0;       prop_high[5] = 800.0;
        prop_low[6] = -6.0;     prop_high[6] = 6.0;
        prop_low[7] = -6.0;     prop_high[7] = 6.0;
        prop_low[8] = -6.0;      prop_high[8] = 6.0;
        prop_low[9] = 0.0;       prop_high[9] = 6.0;         
  }
  if(particle_id == 3 ){

        prop_low[0] = -2000.0;   prop_high[0] = 2000.0;
        prop_low[1] = -2000.0;   prop_high[1] = 2000.0;
        prop_low[2] = -2500.0;   prop_high[2] = 2500.0;
        prop_low[3] = 0.0;       prop_high[3] = 2500.0;
     
        prop_low[4] = m_top - delta_mtop -5.0;       prop_high[4] =  m_top + delta_mtop + 5.0;
        
        prop_low[5] = 0.0;       prop_high[5] = 1300.0;
        prop_low[6] = -6.0;     prop_high[6] = 6.0;
        prop_low[7] = -6.0;     prop_high[7] = 6.0;
        prop_low[8] = -6.0;      prop_high[8] = 6.0;
        prop_low[9] = 0.0;       prop_high[9] = 6.0;         
  }


  if(particle_id == 4 ){

        prop_low[0] = -2000.0;   prop_high[0] = 2000.0;
        prop_low[1] = -2000.0;   prop_high[1] = 2000.0;
        prop_low[2] = -2500.0;   prop_high[2] = 2500.0;
        prop_low[3] = 0.0;       prop_high[3] = 2500.0;
     
        prop_low[4] = 100;       prop_high[4] =  500;
        prop_low[5] = 0.0;       prop_high[5] = 1300.0;
        prop_low[6] = -6.0;     prop_high[6] = 6.0;
        prop_low[7] = -6.0;     prop_high[7] = 6.0;
        prop_low[8] = -6.0;      prop_high[8] = 6.0;
        prop_low[9] = 0.0;       prop_high[9] = 6.0;         
  }


  
  int prop_bin[n_prop]  =  {100,  100,  100, 100, 100, 50 , 100, 100, 50, 50, 50 };

  char name[100];
  sprintf(name,"%s_%s%s",particle_name.c_str(), property[iprop].c_str(), particle_surname.c_str());  
  TH1F *hist  = new TH1F(name,name, prop_bin[iprop], prop_low[iprop], prop_high[iprop] );

  for(int i=0;i<prop[iprop].size();i++){
    hist->Fill( prop[iprop][i] ) ;
  }

  //sprintf( name, "%s", propXaxis[iprop].c_str() );   
  sprintf( name, "%s^{%s}", propXaxis[iprop].c_str(), particle_name.c_str() ); 
  hist->GetXaxis()->SetTitle(name);
  
  return hist;
}

  
string ParticleProperty::GetPropName(int iprop){ return property[iprop]; }
string ParticleProperty::GetPartName(){ return particle_name; }
string ParticleProperty::GetPartSurname(){ return particle_surname; }

TLorentzVector ParticleProperty::GetLorentzVector(int i){
  TLorentzVector v;
  v.SetPxPyPzE( prop[0][i], prop[1][i], prop[2][i], prop[3][i]);
  return v;
}

vector<TLorentzVector> ParticleProperty::GetLorentzVector(){

  vector<TLorentzVector> v_vec;
  TLorentzVector v;
  
  for(int i=0; i<prop[0].size(); i++){
    v.SetPxPyPzE( prop[0][i], prop[1][i], prop[2][i], prop[3][i]);
    v_vec.push_back(v);
  }
  return v_vec;
}












