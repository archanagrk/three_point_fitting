#ifndef __FF_QSQ_FITTING_H__
#define __FF_QSQ_FITTING_H__

#include "fitting_lib/data_selectors.h"
#include "fitting_lib/fit_quality_factory.h"
#include "formfac_qsq_functions_factory.h"
#include "fitting_lib/fit_selector.h"

#include"semble/semble_file_management.h"
#include "semble/semble_meta.h"

using namespace SEMBLE;

/*
   Z(t) functions -- with different source sink seperations - abscissa pair for delta_t and t
*/

//************************************************************************
// FITTING Z(T) FUNCTIONS
//************************************************************************

struct fit_ff_qsq_control{
  double qsqmin, qsqmax;

  double m_p;
  int Nq_min;
  string fit_form;

  double t0 = 0;
  double tcut;
  int n_max;
  double fbound;

  double msq_start;
  double msq_err;

  double f0_start;
  double f0_err;

  bool minos = false;
  bool correlated = true;
  
  bool long_log = false;
  bool plots = true;
  bool topt = true;
};

//************************************************************************

struct fit_ff_qsq_output{
  bool success;
  string fit_summary;
  
  EnsemReal msq;
  string msq_fit_variation;
  
  string fit_long_log;
  string plot_data;

  map<double, pair<double,double> > best_data_description; // can use this for outlier elimination
}; 

//************************************************************************
fit_ff_qsq_output fit_ff_qsq( const Data& data,                        /* the z(t) data */
				    const vector<bool> accepted_data,        /* times which passed noise cut etc... */
				    fit_ff_qsq_control control,
				    FitQuality* fit_qual,                    /* how to compare different fits */
				    double chisq_ndof_cutoff                 /* accept these fits */
				    );
//************************************************************************

//************************************************************************
// PLOTTING UTILITIES
//************************************************************************

struct plot_qsq_function_data{
  int n_points; double qsqmin, qsqmax;
  vector< pair<double,double> > central, upper, lower;
};

plot_qsq_function_data plot_qsq_ensem_function( FFqsqFunction* fn,
							    const map<string, ensem_param_value>& ensem_pars,
							    const vector<double>& qsq_values );

/* central value of average fits */
vector< pair<double,double> > plot_qsq_function( FFqsqFunction* fn,
						       const map<string, param_value>& pars,
						       const vector<double>& qsq_values );




#endif
