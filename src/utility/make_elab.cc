#include "make_elab.h"


int main(int argc, char** argv)
{


  string error_msg = "make_ecm <Inifile> \n ";
  int n_args = 1;
  
  if( argc < n_args + 1 ){ cerr << error_msg; exit(1); }
  
  string xmlini;  {istringstream a(argv[1]); a >> xmlini;};  
  
  //============================
  //===== READ THE XML =====



  int L;
  double xi;
  double PI = (atan(double(1)) * double(4.0));
  string name, Elab;

  Array1dO<Real> momsq;
  Array1dO<string> irrep;  

  XMLReader xml_in(xmlini); 

  cout << "reading: " << xmlini << endl;

  try{
    read(xml_in,"/MakeELab/momO",Elab);
    read(xml_in,"/MakeELab/L",L);
    read(xml_in,"/MakeELab/Xi",xi);
    read(xml_in,"/MakeELab/name",name); 
    read(xml_in,"/MakeELab/irrep",irrep); 
    read(xml_in,"/MakeELab/momSq",momsq); 
  }

  catch( const string& error ){
    cerr << "Error reading input file : " << error << endl;
  }

  //============================
  //===== MAKE THE ELAB =====

  EnsemReal ElabO;

  read(Elab, ElabO);

  for(int i = 1; i <= irrep.size(); i++){
    double mom_sq;
    
    mom_sq = pow((2*PI)/(L*xi), 2) * SEMBLE::toScalar(momsq[i]);

    EnsemReal Elab_this_mom = make_elab(ElabO,mom_sq);

    ostringstream outfile; outfile << name << irrep[i] << ".jack";  write(outfile.str(), Elab_this_mom );

    cout << "done writing: " << outfile.str() << endl;

  }

}

EnsemReal make_elab(EnsemReal ElabO,double momsq ){

  EnsemReal Elab; Elab.resize(ElabO.size());

  for(int bin = 0; bin < Elab.size(); bin++){
    Elab.elem(bin) = sqrt(pow(SEMBLE::toScalar(ElabO.elem(bin)),2) + momsq);
  }

  return Elab;

}

 

