#ifndef __Z_TIMESLICE_FUNCTIONS_H__
#define __Z_TIMESLICE_FUNCTIONS_H__

#include "fitting_lib/functions.h"
#include "fitting_lib/tools.h"
#include "fitting_lib/params.h"



//**************************************************************************************
/* class derived from Function
   only difference is the addition of 
   an operator()(double t, ...)
   which can be called to plot continuous curves
   static_cast<ZtimesliceFunction*>(function) can be used to convert Function pointers
*/
//***************************************************************************************

class ZtimesliceFunction : public Function {
 public:
  virtual ~ZtimesliceFunction(){};

  /* still not specified */
  virtual string name() const = 0;                                                        /* you must provide a function name string */
  virtual double operator()( const Abscissa& x, const mapstringdouble& pars ) const = 0;  /* you must provide a function returing a double */

  /* the actual function that does the work -- can take a double as the time -- this is useful for plotting continuous curves */
  virtual double operator()( double t, const mapstringdouble& pars ) const = 0;

};




//~~~~~~~~~~~~~~~~~~~~~~~~~~~
// sum of exp for prin corrs
//~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* an xml control struct for factory construction */
struct z_timeslice_nexp_params{
  z_timeslice_nexp_params(){};
  z_timeslice_nexp_params(XMLReader& xml_in, const string& path);
  
  int n_exp, t0;
};


/* fit n_exp exponentials to a prin_corr, enforcing c(t0) = 1 */
class ZtimesliceNExp : public ZtimesliceFunction{
 public:
  ZtimesliceNExp( int n_exp_, int t0_ );                  /* construct with explicit integers */
  ZtimesliceNExp( z_timeslice_nexp_params p );    /* construct by factory approach    */
  
  string name() const { return the_name; };
  double operator()( const Abscissa& x, const mapstringdouble& pars ) const;

  /* the actual function that does the work -- can take a double as the time -- this is useful for plotting */
  double operator()( double t, const mapstringdouble& pars ) const ;
  
 private:
  string the_name; 
  int n_exp, t0;
  void initialize(); /* builds the name & logs the param names */
};







#endif
