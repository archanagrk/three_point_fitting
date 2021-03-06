#ifndef __THREE_POINT_TIMESLICE_FITTING_H__
#define __THREE_POINT_TIMESLICE_FITTING_H__

#include "fitting_lib/data_selectors.h"
#include "three_pt_fit_quality_factory.h"
#include "three_point_timeslice_corr_functions_factory.h"
#include "fitting_lib/fit_selector.h"


/*
   Three-point functions -- with different source sink seperations - abscissa pair for delta_t and t
*/

//************************************************************************
// FITTING THREE POINT FUNCTIONS
//************************************************************************


struct fit_three_point_control{
  vector<std::tuple<int,int,int>> tmin_max;
  vector<int> dt;
  vector<double> Ei_min, Ef_min, Ei_max, Ef_max, Nt_min;

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
  bool fixed_ranges = false;
  bool limit_max_fcn_calls = true;

  int max_fcn_calls = 100000; //on an average ~ 1500 fcn calls max I have seen is ~3000
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


struct single_dt_avg_fit_output{
  bool success;
  std::map<AvgFit*, pair<pair<int,int>,pair<int,int>> > avg_fits;
  int dt;
};

//************************************************************************
// DRIVER OVER DIFFERENT SOURCE-SINK SEPERATIONS
//************************************************************************


fit_three_point_output fit( const vector<Data>& data,                        /* the three_point data */
				    const vector<vector<bool>> accepted_data,                        /* times which passed noise cut etc... */
				    vector<fit_three_point_control> control,
				    FitQuality* fit_qual,                                            /* how to compare different fits */
				    double chisq_ndof_cutoff, double num                             /* accept these fits */
				    );
//************************************************************************
// FIX THE RANGES FOR EACH DT
//************************************************************************


single_dt_avg_fit_output single_dt_fit( const Data& data,                       /* the three_point data */
				    const vector<bool> accepted_data,                                   /* times which passed noise cut etc... */
				    fit_three_point_control& control,
				    FitQuality* fit_qual,                                               /* how to compare different fits */
				    double chisq_ndof_cutoff                                           /* accept these fits */
				    );
//************************************************************************
// FIT EACH SOURCE-SINK SEPERATION
//************************************************************************


fit_three_point_output fit_three_point_corr( const Data& data,                        /* the three_point data */
				    const vector<bool> accepted_data,                                         /* times which passed noise cut etc... */
				    fit_three_point_control& control,
				    FitQuality* fit_qual, vector<single_dt_avg_fit_output> single_dt_fits,    /* how to compare different fits */
				    double chisq_ndof_cutoff                                                  /* accept these fits */
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


//************************************************************************
// FUNCTION TO LOOP OVER ALL THE T_MINS AND T_MAXS
//************************************************************************


std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> get_all_t_ranges(fit_three_point_control control, bool slide, bool src, bool snk,
                      std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> range, int count_dt);

std::vector<pair<pair<int,int>,pair<int,int>>> get_single_dt_t_ranges(fit_three_point_control control, bool slide, bool src, bool snk);

 void cart_product(
    std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >>& rvvi,  // final result
    std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>>>&  rvi,   // current result 
    std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >>::const_iterator me, // current input
    std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >>::const_iterator end); // final input


#endif
