#ifndef __THREE_POINT_TIMESLICE_FITTING_H__
#define __THREE_POINT_TIMESLICE_FITTING_H__

#include "data_selectors.h"
#include "fit_quality_factory.h"
#include "functions_factory.h"
#include "fit_selector.h"

/*
   Three-point functions -- with different source sink seperations - abscissa pair for delta_t and t
*/

//************************************************************************
// FITTING THREE POINT FUNCTIONS
//************************************************************************

struct fit_three_point_control{
  vector<int> tmin, tmax; 
  vector<int> dt;
  int Nt_min = 1;

  vector<int> tsrc, tsnk;

  double F_start = 0.005;
  double F_err   = 0.7;
  
  bool only_cnst = false;
  bool src_exp = false;
  bool only_one_exp = false;
  bool minos = false;
  bool correlated = true;
  
  bool long_log = false;
  bool plots = false;
};

//************************************************************************

struct fit_three_point_output{
  bool success;
  string fit_summary;
  
  EnsemReal F;
  string F_fit_variation;
  
  string fit_long_log;
  string plot_data;

  map<pair<int,int>, pair<double,double> > best_data_description; // can use this for outlier elimination
}; 

//************************************************************************
fit_three_point_output fit_three_point_corr( const Data& data,                        /* the three_point data */
				    const vector<bool> accepted_data,        /* times which passed noise cut etc... */
				    fit_three_point_control control,
				    FitQuality* fit_qual,                    /* how to compare different fits */
				    double chisq_ndof_cutoff                 /* accept these fits */
				    );
//************************************************************************

//************************************************************************
// PLOTTING UTILITIES
//************************************************************************

struct plot_three_point_timeslice_function_data{
  vector<int> n_points; vector<double> tmin, tmax;
  vector<double> dt;
  vector< pair<pair<double,double>,double> > central, upper, lower;
};

plot_three_point_timeslice_function_data plot_three_point_timeslice_ensem_function( ThreePointtimesliceCorrFunction* fn,
							    const map<string, ensem_param_value>& ensem_pars,
							    const vector<std::pair<double,double>>& t_values, const vector<double>& delt );

/* central value of average fits */
vector< pair<pair<double,double>,double> >  plot_three_point_timeslice_function( ThreePointtimesliceCorrFunction* fn,
						       const map<string, param_value>& pars,
						       const vector<std::pair<double,double>>& t_values);




#endif
