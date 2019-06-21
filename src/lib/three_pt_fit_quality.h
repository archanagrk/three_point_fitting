#ifndef __THREE_PT_FIT_QUALITY_H__
#define __THREE_PT_FIT_QUALITY_H__

#include "fitting_lib/avg_fitter.h"
#include "fitting_lib/fit_quality.h"
#include "pair_int_abscissa.h"

#include "three_pt_fit_quality_factory.h"
#include "three_point_timeslice_corr_functions_factory.h"

/* build a factory of ThreePtFitQuality types */
/* many will be specific to particular fit functions */
/* the AvgFit passed in contains the data and the function, 
   so the quality can compute all kinds of stuff */

/*
  the fit quality constructor can in principle depend upon some fixed parameters
  so probably want an xml control or similar to smuggle them in via the factory
*/

bool cmp(double a, double b);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Estimate the errors in the fits better
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* an xml control struct for factory construction */
struct gen3pt_params{
  gen3pt_params(){};
  gen3pt_params(XMLReader& xml_in, const string& path);
  
  double power;
  double power_exp;
  double power_time;
  bool multiply_exp;
};



class Generic3pt : public FitQuality{

 public:
 Generic3pt( double power_, double power_exp_, double power_time_, bool multiply_exp_ ) : power(power_), power_exp(power_exp_), power_time(power_time_), multiply_exp(multiply_exp_) {};

 Generic3pt( gen3pt_params p) : power( p.power ), power_exp(p.power_exp), power_time(p.power_time), multiply_exp( p.multiply_exp )
    {
      cout << "constructed a Generic_3pt with power = " << power;
      if(multiply_exp){ cout << ", multiplying by mass"; }
      cout << endl;
    };
  
  double operator()( const AvgFit& fit ) const;
  string name() const { return "generic"; } 
  
 private:
  double power;
  double power_exp;
  double power_time;
  bool multiply_exp;

};


class QNGen : public FitQuality{

 public:
 QNGen( double power_, double power_exp_, double power_time_, bool multiply_exp_ ) : power(power_), power_exp(power_exp_), power_time(power_time_), multiply_exp(multiply_exp_) {};

 QNGen( gen3pt_params p) : power( p.power ), power_exp(p.power_exp), power_time(p.power_time), multiply_exp( p.multiply_exp )
    {
      cout << "constructed a QNGen with power = " << power;
      if(multiply_exp){ cout << ", multiplying by mass"; }
      cout << endl;
    };
  
  double operator()( const AvgFit& fit ) const;
  string name() const { return "qn_gen"; } 
  
 private:
  double power;
  double power_exp;
  double power_time;
  bool multiply_exp;

};

#endif

