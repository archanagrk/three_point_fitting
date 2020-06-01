#ifndef __Z_TIMESLICE_FITTING_H__
#define __Z_TIMESLICE_FITTING_H__

#include "fitting_lib/data_selectors.h"
#include "fitting_lib/fit_quality_factory.h"
#include "z_timeslice_functions_factory.h"
#include "fitting_lib/fit_selector.h"

/*
   Z(t) functions -- with different source sink seperations - abscissa pair for delta_t and t
*/

//************************************************************************
// FITTING Z(T) FUNCTIONS
//************************************************************************

struct fit_z_timeslice_control{
  int tmin, tmax;
  int Nt_min = 1;

  double A_start = 0.4;
  double A_err   = 0.4;
  
  bool only_cnst = false;
  bool only_one_exp = false;
  bool minos = false;
  bool correlated = true;
  
  bool long_log = false;
  bool plots = false;
};

//************************************************************************

struct fit_z_timeslice_output{
  bool success;
  string fit_summary;
  
  EnsemReal m;
  string m_fit_variation;
  
  string fit_long_log;
  string plot_data;

  map<int, pair<double,double> > best_data_description; // can use this for outlier elimination
}; 

//************************************************************************
fit_z_timeslice_output fit_z_timeslice( const Data& data,                        /* the z(t) data */
				    int t0,
				    const vector<bool> accepted_data,        /* times which passed noise cut etc... */
				    fit_z_timeslice_control control,
				    FitQuality* fit_qual,                    /* how to compare different fits */
				    double chisq_ndof_cutoff                 /* accept these fits */
				    );
//************************************************************************

//************************************************************************
// PLOTTING UTILITIES
//************************************************************************

struct plot_timeslice_function_data{
  int n_points; double tmin, tmax;
  vector< pair<double,double> > central, upper, lower;
};

plot_timeslice_function_data plot_timeslice_ensem_function( ZtimesliceFunction* fn,
							    const map<string, ensem_param_value>& ensem_pars,
							    const vector<double>& t_values );

/* central value of average fits */
vector< pair<double,double> > plot_timeslice_function( ZtimesliceFunction* fn,
						       const map<string, param_value>& pars,
						       const vector<double>& t_values );




#endif
