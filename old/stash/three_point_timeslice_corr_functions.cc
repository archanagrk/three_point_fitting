#include "three_point_timeslice_corr_functions.h"

/* pass in the control objects */
three_point_timeslice_corr_nexp_params::three_point_timeslice_corr_nexp_params(XMLReader& xml, const string& path){
  try {
    XMLReader paramtop(xml, path); 
    read(paramtop, "src_exp", src_exp);
    read(paramtop, "snk_exp", snk_exp);
  }
  catch(const std::string& e) 
    {  std::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;  exit(1);   }
};


/* construct using struct */
ThreePointtimesliceCorrNExp::ThreePointtimesliceCorrNExp( three_point_timeslice_corr_nexp_params p ) 
  :  src_exp(p.src_exp), snk_exp(p.snk_exp){ initialize(); }; 


/* construct using integers */
ThreePointtimesliceCorrNExp::ThreePointtimesliceCorrNExp( int src_exp_, int snk_exp_ ) 
  : src_exp(src_exp_), snk_exp(snk_exp_){ initialize(); };  


/* make the function name & log the param names */
void ThreePointtimesliceCorrNExp::initialize(){ 
  {stringstream ss; ss << src_exp << "src_exp" << snk_exp << "exp_pc"; the_name = ss.str();}
  
  par_names.push_back( "F" );

  if(snk_exp == 1){par_names.push_back( "dmf" ); par_names.push_back( "Ff" );}

  if(src_exp == 1) {par_names.push_back( "dmi" ); par_names.push_back( "Fi" );}

  if(src_exp >= 2){
    for(int i = 1; i <= snk_exp; i++){
      { stringstream ss; ss << "dFi" << i;  par_names.push_back( ss.str() ); }
      { stringstream ss; ss << "dmi" << i; par_names.push_back( ss.str() ); }
    }
    
  }


  if(snk_exp >= 2){
    par_names.push_back( "F" );
    for(int i = 1; i < snk_exp; i++){
      { stringstream ss; ss << "dFf" << i;  par_names.push_back( ss.str() ); }
      { stringstream ss; ss << "dmf" << i; par_names.push_back( ss.str() ); }
    }

  }

}


/* return the function value taking a double as the timeslice argument */
double ThreePointtimesliceCorrNExp::operator()( std::pair<double,double> t, const mapstringdouble& pars ) const {
  
  /* mass for each new exp is the previous exp mass plus a dm */
  /* can limit the dm to positive values to ensure ordering of masses */
  

  double t_snk = t.first;
  double t_curr = t.second;

  double F = pars.find("F")->second;

    
  if(src_exp == 0){

    if(snk_exp == 0){
      return F;}
    
    else if(snk_exp == 1){

      double Ff = pars.find("Ff")->second;
      double dmf = pars.find("dmf")->second; 

      return F + Ff * exp( - dmf * (t_snk - t_curr) );
    }
    
  }



  else if(src_exp == 1){

      double Fi = pars.find("Fi")->second;
      double dmi = pars.find("dmi")->second;

      if(snk_exp == 0){

        return F + Fi * exp( - dmi * t_curr );
      }

      else if(snk_exp == 1){

        double Ff = pars.find("Ff")->second;
        double dmf = pars.find("dmf")->second; 

        return F + Fi*exp( - dmi * t_curr ) + Ff*exp( - dmf * (t_snk - t_curr) ) ;
      }
  }


  
  else{
    
    double out = 0.0;
  
    for(int i = 1; i <= src_exp; i++){
      double dFi, dmi;
      { stringstream ss; ss << "dFi" << i;   dFi = pars.find(ss.str())->second; }
      { stringstream ss; ss << "dmi" << i; dmi = pars.find(ss.str())->second; }
      out += dFi*exp( - dmi *  t_curr );
    }

    for(int i = 1; i <= snk_exp; i++){
      double dFf, dmf;
      { stringstream ss; ss << "dFf" << i;   dFf = pars.find(ss.str())->second; }
      { stringstream ss; ss << "dmf" << i; dmf = pars.find(ss.str())->second; }
      out += dFf*exp( - dmf *  t_snk );
    }
  
    out += F ;
    return out;
  }
}


/* return the function value taking an Abscissa as the timeslice argument */
double ThreePointtimesliceCorrNExp::operator()(const Abscissa& x, const mapstringdouble& pars ) const{

  if( x.abscissa_type() != "pair_int_abscissa" )
    { cerr << "ThreePointtimesliceCorrNExp:: supplied " << x.abscissa_type() << " when pair_int_abscissa required, exiting" << endl; exit(1); }
  
  /* t value must actually be an PairIntAbscissa */
  std::pair<double,double> t = (static_cast<const PairIntAbscissa&>(x)).get_x();

  
  check_pars( pars ); // 1 microsec
  
  return operator()( t, pars );
}
