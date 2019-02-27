#include "z_timeslice_functions.h"

/* pass in the control objects */
z_timeslice_nexp_params::z_timeslice_nexp_params(XMLReader& xml, const string& path){
  try {
    XMLReader paramtop(xml, path); 
    read(paramtop, "n_exp", n_exp);
    read(paramtop, "t0", t0);
  }
  catch(const std::string& e) 
    {  std::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;  exit(1);   }
};


/* construct using struct */
ZtimesliceNExp::ZtimesliceNExp( z_timeslice_nexp_params p ) 
  : n_exp(p.n_exp), t0(p.t0){ initialize(); }; 


/* construct using integers */
ZtimesliceNExp::ZtimesliceNExp( int n_exp_, int t0_ ) 
  : n_exp(n_exp_), t0(t0_){ initialize(); };  


/* make the function name & log the param names */
void ZtimesliceNExp::initialize(){ 
  { stringstream ss; ss << n_exp << "exp_pc"; the_name = ss.str();} 

  par_names.push_back( "A" );
  if( n_exp == 1 ){ par_names.push_back( "m" ); par_names.push_back( "B" ); }
  else if( n_exp == 2 ){ par_names.push_back( "m" ); par_names.push_back( "B" );  par_names.push_back( "dm" ); par_names.push_back( "C" ); }
  else if( n_exp > 2){
    for(int i = 1; i <= n_exp; i++){
      { stringstream ss; ss << "A" << i;  par_names.push_back( ss.str() ); }
      { stringstream ss; ss << "dm" << i; par_names.push_back( ss.str() ); }
    }
  }
}


/* return the function value taking a double as the timeslice argument */
double ZtimesliceNExp::operator()( double t, const mapstringdouble& pars ) const {
  
  /* mass for each new exp is the previous exp mass plus a dm */
  /* can limit the dm to positive values to ensure ordering of masses */
  
  double t_min_t0 = t - t0;
  

  if(n_exp == 0){

    double A = pars.find("A")->second;
    
    return A;     
  }
  if(n_exp == 1){ 

    double m = pars.find("m")->second;
    double A = pars.find("A")->second;
    double B = pars.find("B")->second;

    return A + B * exp( - m * t_min_t0 );
  }
  else if(n_exp == 2){

    double m = pars.find("m")->second; 
    double dm = pars.find("dm")->second; 
    double A = pars.find("A")->second; 
    double B = pars.find("B")->second; 
    double C = pars.find("C")->second; 

    return A + B*exp( - m * t_min_t0 ) + C*exp( - (m + dm) * t_min_t0 );
  }
  else{
    double m = pars.find("m")->second;  double mm = m;
    double AA = pars.find("A")->second;
    double out = 0.0; 
  
    for(int i = 1; i <= n_exp; i++){
      double A, dm;
      { stringstream ss; ss << "A" << i;   A = pars.find(ss.str())->second; }
      { stringstream ss; ss << "dm" << i; dm = pars.find(ss.str())->second; }
      mm += dm;
      out += A*exp( - mm *  t_min_t0 );
    }
  
    out += AA;
    return out;
  }
}


/* return the function value taking an Abscissa as the timeslice argument */
double ZtimesliceNExp::operator()(const Abscissa& x, const mapstringdouble& pars ) const{

  if( x.abscissa_type() != "int_abscissa" )
    { cerr << "ZtimesliceNExp:: supplied " << x.abscissa_type() << " when int_abscissa required, exiting" << endl; exit(1); }
  
  /* t value must actually be an IntAbscissa */
  int t = (static_cast<const IntAbscissa&>(x)).get_x();
  
  check_pars( pars ); // 1 microsec
  
  return operator()( double(t), pars );
}

