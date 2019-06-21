#ifndef __THREE_PT_FIT_QUALITY_H__
#define __THREE_PT_FIT_QUALITY_H__

#include "avg_fitter.h"

/* build a factory of ThreePtFitQuality types */
/* many will be specific to particular fit functions */
/* the AvgFit passed in contains the data and the function, 
   so the quality can compute all kinds of stuff */

/*
  the fit quality constructor can in principle depend upon some fixed parameters
  so probably want an xml control or similar to smuggle them in via the factory
*/

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Estimate the errors in the fits better
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* an xml control struct for factory construction */
struct gen3pt_params{
  gen3pt_params(){};
  gen3pt_params(XMLReader& xml_in, const string& path);
  
  int power;
  int power_mass;
  int power_F;
  int power_time;
  bool multiply_mass;
};



class Generic3pt : public FitQuality{
 public:
 Generic3pt( int power_, int power_mass_,int power_F_, int power_time_, bool multiply_mass_ ) : power(power_), power_mass(power_mass_),  power_F(power_F_), power_time(power_time_), multiply_mass(multiply_mass_) {};
 Generic3pt( gen3pt_params p) : power( p.power ), power_mass(p.power_mass), power_F(p.power_F), power_time(p.power_time), multiply_mass( p.multiply_mass )
    {
      cout << "constructed a Generic_3pt with power = " << power;
      if(multiply_mass){ cout << ", multiplying by mass"; }
      cout << endl;
    };
  
  double operator()( const AvgFit& fit ) const;
  string name() const { return "generic"; } 
  
 private:
  int power;
  int power_mass;
  int power_F;
  int power_time;
  bool multiply_mass;
};

#endif

