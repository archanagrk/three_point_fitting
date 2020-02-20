#include "three_point_timeslice_corr_functions.h"

/* pass in the control objects */
three_point_timeslice_corr_nexp_params::three_point_timeslice_corr_nexp_params(XMLReader& xml, const string& path){
  try {
    XMLReader paramtop(xml, path); 
    read(paramtop, "src_exp", src_exp);
    read(paramtop, "snk_exp", snk_exp);
    read(paramtop, "dt", dt);
  }
  catch(const std::string& e) 
    {  std::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;  exit(1);   }
};


/* construct using struct */
ThreePointtimesliceCorrNExp::ThreePointtimesliceCorrNExp( three_point_timeslice_corr_nexp_params p ) 
  :  src_exp(p.src_exp), snk_exp(p.snk_exp), dt(p.dt){ initialize(); }; 


/* construct using vector of integers */
ThreePointtimesliceCorrNExp::ThreePointtimesliceCorrNExp( vector<int> src_exp_, vector<int> snk_exp_, vector<int> dt_ ) 
  : src_exp(src_exp_), snk_exp(snk_exp_), dt(dt_){ initialize(); };  


/* make the function name & log the param names */
void ThreePointtimesliceCorrNExp::initialize(){

  {stringstream ss; 
    for(size_t t = 0; t < dt.size(); t++){
    ss << "src_exp_" << src_exp[t] << "_snk_exp_" << snk_exp[t] << "_dt_" << dt[t]; the_name = ss.str();
    }}
  
  /* all functions have F */
  par_names.push_back( "F" );


  /* loop over the number of source-sink seperations and independently choose the number of source exp, sink exp in case */
  for(size_t t = 0; t < dt.size(); t++){

    if(snk_exp[t] == 1){    
      { stringstream ss; ss << "Ff" << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
      { stringstream ss; ss << "Ef" << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }
    }

    if(src_exp[t] == 1) {
      { stringstream ss; ss << "Fi" << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
      { stringstream ss; ss << "Ei" << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }      
    }

    if(src_exp[t] >= 2){
      for(size_t i = 1; i <= snk_exp[t]; i++){
        { stringstream ss; ss << "Fi" << i << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
        { stringstream ss; ss << "Ei" << i << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }
      }
      
    }


    if(snk_exp[t] >= 2){
      par_names.push_back( "F" );
      for(size_t i = 1; i < snk_exp[t]; i++){
        { stringstream ss; ss << "Ff" << i << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
        { stringstream ss; ss << "Ef" << i << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }
      }

    }
  }


}


/* return the function value taking a pair of double as the timeslice argument */
double ThreePointtimesliceCorrNExp::operator()( std::pair<double,double> t, const mapstringdouble& pars ) const {
  
/* can limit the Ei and Ef to positive values to ensure ordering of masses */
  

  double t_snk = t.first;
  double t_curr = t.second;
  int Dt = 0;

  double F = pars.find("F")->second;

  /* all functions have F */
  double out = F;

  /* loop over the number of source-sink seperations and independently choose the function based on number of source exp, sink exp in case */
  for(size_t i = 0; i < dt.size(); i++){ if(dt[i] == t_snk){ Dt = i; break;} }

  if(src_exp[Dt] == 0){

    
    if(snk_exp[Dt] == 1){

      stringstream ff; ff << "Ff" << "_Dt" << dt[Dt];
      stringstream ef; ef << "Ef" << "_Dt" << dt[Dt];

      double Ff = pars.find(ff.str())->second;
      double Ef = pars.find(ef.str())->second; 

      out += Ff * exp( - Ef * (t_snk - t_curr) ); 
    }
    
  }



  else if(src_exp[Dt] == 1){

    stringstream fi; fi << "Fi" << "_Dt" << dt[Dt];
    stringstream ei; ei << "Ei" << "_Dt" << dt[Dt];

    double Fi = pars.find(fi.str())->second;
    double Ei = pars.find(ei.str())->second; 

    out += Fi*exp( - Ei * t_curr );

    if(snk_exp[Dt] == 1){

      stringstream ff; ff << "Ff" << "_Dt" << dt[Dt];
      stringstream ef; ef << "Ef" << "_Dt" << dt[Dt];

      double Ff = pars.find(ff.str())->second;
      double Ef = pars.find(ef.str())->second; 

      out += Ff*exp( - Ef * (t_snk - t_curr) ) ;
    }

  }


  
  else{

    for(int i = 1; i <= src_exp[Dt]; i++){
      double Fi, Ei;
      { stringstream ss; ss << "Fi" << i << "_Dt" << dt[Dt];   Fi = pars.find(ss.str())->second; }
      { stringstream ss; ss << "Ei" << i << "_Dt" << dt[Dt]; Ei = pars.find(ss.str())->second; }
      out += Fi*exp( - Ei *  t_curr );
    }

    for(int i = 1; i <= snk_exp[Dt]; i++){
      double Ff, Ef;
      { stringstream ss; ss << "Ff" << i << "_Dt" << dt[Dt];   Ff = pars.find(ss.str())->second; }
      { stringstream ss; ss << "Ef" << i << "_Dt" << dt[Dt]; Ef = pars.find(ss.str())->second; }
      out += Ff*exp( - Ef *  t_snk );
    }

  }


  return out;
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
