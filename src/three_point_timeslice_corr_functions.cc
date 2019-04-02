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


/* construct using integers */
ThreePointtimesliceCorrNExp::ThreePointtimesliceCorrNExp( int src_exp_, int snk_exp_, vector<int> dt_ ) 
  : src_exp(src_exp_), snk_exp(snk_exp_), dt(dt_){ initialize(); };  


/* make the function name & log the param names */
void ThreePointtimesliceCorrNExp::initialize(){ 
  {stringstream ss; ss << src_exp << "src_exp" << snk_exp << "snk_exp"; the_name = ss.str();}
  
  par_names.push_back( "F" );

  for(int t = 0; t < dt.size(); t++){

    if(snk_exp == 1){    
      { stringstream ss; ss << "Ff" << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
      { stringstream ss; ss << "Ef" << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }
    }

    if(src_exp == 1) {
      { stringstream ss; ss << "Fi" << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
      { stringstream ss; ss << "Ei" << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }      
    }

    if(src_exp >= 2){
      for(int i = 1; i <= snk_exp; i++){
        { stringstream ss; ss << "Fi" << i << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
        { stringstream ss; ss << "Ei" << i << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }
      }
      
    }


    if(snk_exp >= 2){
      par_names.push_back( "F" );
      for(int i = 1; i < snk_exp; i++){
        { stringstream ss; ss << "Ff" << i << "_Dt" << dt[t];  par_names.push_back( ss.str() ); }
        { stringstream ss; ss << "Ef" << i << "_Dt" << dt[t]; par_names.push_back( ss.str() ); }
      }

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

  double out = F;
  
  if(src_exp == 0){

    if(snk_exp == 0){
      return out;}
    
    else if(snk_exp == 1){
      for(int t = 0; t < dt.size(); t++){

        stringstream ff; ff << "Ff" << "_Dt" << dt[t];
        stringstream ef; ef << "Ef" << "_Dt" << dt[t];

        double Ff = pars.find(ff.str())->second;
        double Ef = pars.find(ef.str())->second; 

        out += Ff * exp( - Ef * (t_snk - t_curr) ); 

      }
      return out;
    }
    
  }



  else if(src_exp == 1){

    for(int t = 0; t < dt.size(); t++){

      stringstream fi; fi << "Fi" << "_Dt" << dt[t];
      stringstream ei; ei << "Ei" << "_Dt" << dt[t];

      double Fi = pars.find(fi.str())->second;
      double Ei = pars.find(ei.str())->second; 

      out += Fi*exp( - Ei * t_curr );

      if(snk_exp == 0){continue;}

      else if(snk_exp == 1){

        stringstream ff; ff << "Ff" << "_Dt" << dt[t];
        stringstream ef; ef << "Ef" << "_Dt" << dt[t];

        double Ff = pars.find(ff.str())->second;
        double Ef = pars.find(ef.str())->second; 

        out += Ff*exp( - Ef * (t_snk - t_curr) ) ;
      }
    }
    return out;
  }


  
  else{
    
    double out = 0.0;

    for(int t = 0; t < dt.size(); t++){

      for(int i = 1; i <= src_exp; i++){
        double Fi, Ei;
        { stringstream ss; ss << "Fi" << i << "_Dt" << dt[t];   Fi = pars.find(ss.str())->second; }
        { stringstream ss; ss << "Ei" << i << "_Dt" << dt[t]; Ei = pars.find(ss.str())->second; }
        out += Fi*exp( - Ei *  t_curr );
      }

      for(int i = 1; i <= snk_exp; i++){
        double Ff, Ef;
        { stringstream ss; ss << "Ff" << i << "_Dt" << dt[t];   Ff = pars.find(ss.str())->second; }
        { stringstream ss; ss << "Ef" << i << "_Dt" << dt[t]; Ef = pars.find(ss.str())->second; }
        out += Ff*exp( - Ef *  t_snk );
      }
    }

    out += F ;
    return out;
  }

  return 0;
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
