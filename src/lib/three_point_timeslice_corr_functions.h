#ifndef __THREE_POINT_TIMESLICE_CORR_FUNCTIONS_H__
#define __THREE_POINT_TIMESLICE_CORR_FUNCTIONS_H__

#include "fitting_lib/functions.h"
#include "fitting_lib/tools.h"
#include "fitting_lib/params.h"
#include "pair_int_abscissa.h"
#include <numeric>


//**************************************************************************************
/* class derived from Function
   only difference is the addition of 
   an operator()(std::pair<double ,double> t, ...)
   which can be called to plot continuous curves
   static_cast<ThreePointtimesliceCorrFunction*>(function) can be used to convert Function pointers

   Takes in vector<int> dt of different Dt (source-sink seperations) and a vector<int> src_exp and snk_exp which
   correspond to the number of src_exp and snk_exp for each Dt. 
   
   The order of elements is important here: the vectors are sorted in ascending order of Dt.
*/
//***************************************************************************************

class ThreePointtimesliceCorrFunction : public Function {
 public:
  virtual ~ThreePointtimesliceCorrFunction(){};

  /* still not specified */
  virtual string name() const = 0;                                                        /* you must provide a function name string */
  virtual double operator()( const Abscissa& x, const mapstringdouble& pars ) const = 0;  /* you must provide a function returing a double */

  /* the actual function that does the work -- can take a double as the time -- this is useful for plotting continuous curves */
  virtual double operator()( std::pair<double,double> t, const mapstringdouble& pars ) const = 0;

};




//~~~~~~~~~~~~~~~~~~~~~~~~~~~
// sum of exp for 3 pt
//~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* an xml control struct for factory construction */
struct three_point_timeslice_corr_nexp_params{
  three_point_timeslice_corr_nexp_params(){};
  three_point_timeslice_corr_nexp_params(XMLReader& xml_in, const string& path);
  
  vector<int> snk_exp;
  vector<int> src_exp;
  vector<int> dt;
};


/* fit snk_exp exponentials to a 3 pt */
class ThreePointtimesliceCorrNExp : public ThreePointtimesliceCorrFunction{
 public:
  ThreePointtimesliceCorrNExp( vector<int> src_exp_, vector<int> snk_exp_, vector<int> dt_);                  /* construct with explicit integers */
  ThreePointtimesliceCorrNExp( three_point_timeslice_corr_nexp_params p );    /* construct by factory approach    */
  
  string name() const { return the_name; };
  double operator()( const Abscissa& x, const mapstringdouble& pars ) const;

  /* the actual function that does the work -- can take a double as the time -- this is useful for plotting */
  double operator()( std::pair<double,double> t, const mapstringdouble& pars ) const ;
  
 private:
  string the_name; 
  vector<int> snk_exp;
  vector<int> src_exp;
  vector<int> dt;
  
  void initialize(); /* builds the name & logs the param names */
};







#endif
