#include "formfac_qsq_functions.h"

/* pass in the control objects */
ff_qsq_params::ff_qsq_params(XMLReader& xml, const string& path){
  try {
    XMLReader paramtop(xml, path); 
    read(paramtop, "m_p", m_p);
    read(paramtop, "qsq_max", qsq_max);
    read(paramtop, "topt", topt);
    read(paramtop, "tcut", tcut);
    read(paramtop, "n_max", n_max);

  }
  catch(const std::string& e) 
    {  std::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;  exit(1);   }
};

  //************************************************************************************************

/* construct using struct */
FFqsqVMD::FFqsqVMD( ff_qsq_params p ) 
  : m_p(p.m_p){ initialize(); }; 


/* construct using double */
FFqsqVMD::FFqsqVMD( double m_p_ ) 
  : m_p(m_p_){ initialize(); };  


/* make the function name & log the param names */
void FFqsqVMD::initialize(){ 
  { stringstream ss; ss << "m_p_" << m_p ; the_name = ss.str();} 

  par_names.push_back( "f0" );
  par_names.push_back( "msq" );

}


/* return the function value taking a double as the argument */
double FFqsqVMD::operator()( double qsq, const mapstringdouble& pars ) const {
  
  double msq = pars.find("msq")->second;
  double f0 = pars.find("f0")->second;

  return f0/ (1 + (qsq/msq));     
  
}


/* return the function value taking an Abscissa as the argument */
double FFqsqVMD::operator()(const Abscissa& x, const mapstringdouble& pars ) const{

  if( x.abscissa_type() != "double_abscissa" )
    { cerr << "FFqsqVMD:: supplied " << x.abscissa_type() << " when double_abscissa required, exiting" << endl; exit(1); }
  
  /* qsq value must actually be an DoubleAbscissa */
  double qsq = (static_cast<const DoubleAbscissa&>(x)).get_x();
  
  check_pars( pars ); 
  
  return operator()( qsq, pars );
}

  //************************************************************************************************


/* construct using struct */
FFqsqGauss::FFqsqGauss( ff_qsq_params p ) 
  : m_p(p.m_p){ initialize(); }; 


/* construct using double */
FFqsqGauss::FFqsqGauss( double m_p_ ) 
  : m_p(m_p_){ initialize(); };  


/* make the function name & log the param names */
void FFqsqGauss::initialize(){ 
  { stringstream ss; ss << "m_p_" << m_p; the_name = ss.str();} 

  par_names.push_back( "f0" );
  par_names.push_back( "betasq" );

}


/* return the function value taking a double as the argument */
double FFqsqGauss::operator()( double qsq, const mapstringdouble& pars ) const {
  
  double betasq = pars.find("betasq")->second;
  double f0 = pars.find("f0")->second;

  return f0 * exp(-qsq/(16*betasq));     
  
}


/* return the function value taking an Abscissa as the argument */
double FFqsqGauss::operator()(const Abscissa& x, const mapstringdouble& pars ) const{

  if( x.abscissa_type() != "double_abscissa" )
    { cerr << "FFqsqGauss:: supplied " << x.abscissa_type() << " when double_abscissa required, exiting" << endl; exit(1); }
  
  /* qsq value must actually be an DoubleAbscissa */
  double qsq = (static_cast<const DoubleAbscissa&>(x)).get_x();
  
  check_pars( pars );
  
  return operator()( qsq, pars );
}

  //************************************************************************************************


/* construct using struct */
FFqsqZExp::FFqsqZExp( ff_qsq_params p ) 
  : m_p(p.m_p), qsq_max(p.qsq_max), topt(p.topt), tcut(p.tcut), n_max(p.n_max){ initialize(); }; 


/* construct using double */
FFqsqZExp::FFqsqZExp( double m_p_, double qsq_max_, bool topt_, double tcut_, int n_max_ ) 
  : m_p(m_p_), qsq_max(qsq_max_), topt(topt_), tcut(tcut_), n_max(n_max_){ initialize(); };  


/* make the function name & log the param names */
void FFqsqZExp::initialize(){ 
  { stringstream ss; ss << "m_p_" << m_p << "_tcut_" << tcut << "_n_max_" << n_max; the_name = ss.str();} 

  par_names.push_back( "f0" );

  // par_names.push_back( "msq" );

  for(int i = 1; i < n_max; i++ ){ stringstream ss; ss << "f" << i;  par_names.push_back( ss.str() ); }

}


/* return the function value taking a double as the argument */
double FFqsqZExp::operator()( double qsq, const mapstringdouble& pars ) const {

  
  // double msq = pars.find("msq")->second;
  double f0 = pars.find("f0")->second;

  double t0 = 0;
  if(topt){t0 = tcut * (1 - sqrt(1 + (qsq_max/tcut)) );}

  double one_term = sqrt(tcut + qsq);
  double two_term = sqrt(tcut - t0);
  double z = ( one_term  - two_term )/( one_term  + two_term );

  //double coeff = 1/(1 + (qsq/msq) );

  //double out = coeff * f0 ; 
  double out = f0 ; 

  for(int i = 1; i < n_max; i++){
    stringstream f; f << "f" << i; 
    double f_i = pars.find(f.str())->second;

    // out += coeff * f_i * pow(z,i);
    out += f_i * pow(z,i);

    }

  return out;     
  
}


/* return the function value taking an Abscissa as the argument */
double FFqsqZExp::operator()(const Abscissa& x, const mapstringdouble& pars ) const{

  if( x.abscissa_type() != "double_abscissa" )
    { cerr << "FFqsqZExp:: supplied " << x.abscissa_type() << " when double_abscissa required, exiting" << endl; exit(1); }
  
  /* qsq value must actually be an DoubleAbscissa */
  double qsq = (static_cast<const DoubleAbscissa&>(x)).get_x();
  
  check_pars( pars ); 
  
  return operator()( qsq, pars );
}


