#include "z_timeslice_fitting.h"

fit_z_timeslice_output fit_z_timeslice( const Data& data,                        
				    int t0,
				    const vector<bool> accepted_data,       
				    fit_z_timeslice_control control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff               
				    )
{
  /* set minuit controls using defaults */
  MinuitControl minuit_controls; minuit_controls.minos = control.minos;
  
  /* hold the various fits */
  vector<AvgFit*> fits;

  /* the output object */
  fit_z_timeslice_output out;

  //************************************************************************************************

  //************************************************************************************************
      /*
      THIS IS EXP FITTING WITH CONSTANT 
      */
  //************************************************************************************************
 
  Function* cnst = new ZtimesliceNExp( 0 , t0 );  
  Function* cnst_one_exp = new ZtimesliceNExp( 1, t0 );  
  Function* cnst_two_exp = new ZtimesliceNExp( 2, t0 );

  /* perform constant fits */
  {

    map<string, param_value> start_params;
  
    param_value A(control.A_start, control.A_err); 
    A.minos = control.minos;
    start_params.insert( make_pair("A", A ) ); 

  
    vector<bool> previous( accepted_data.size(), false );
    for(int tlow = control.tmax - control.Nt_min; tlow >= control.tmin; tlow--){
      if(tlow == t0){continue;}
    
      vector<bool> active_data;
      {
        Abscissa* x_tlow = new IntAbscissa(tlow - 1);
        active_data = data_above_x( data, x_tlow );
        delete x_tlow;
      }
      active_data = active_data && accepted_data;
      
      if( (count_active(active_data) >= control.Nt_min) && (active_data != previous) )
      {
        stringstream name; name << "cnst_tmin" << tlow;
  
        AvgFit* this_fit = new AvgFit( data, active_data, cnst, start_params, minuit_controls, control.correlated, name.str() );
        fits.push_back( this_fit );
        previous = active_data;
  
        /*cout << "==================================================================================" << endl;
        cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
        cout << data.print_data(active_data) << endl;
        cout << this_fit->report() << endl;
        cout << "==================================================================================" << endl; */
      }
      
    } // next tlow
  } //end constant fits
  //************************************************************************************************

  /* pick the best constant fit to choose a start mass for one_exp fits */

  int tmin_cnst = 0;
  param_value A_cnst(control.A_start, control.A_err); A_cnst.minos = control.minos;

  if( fits.size() > 0 )
  {
    map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
    minuit_fit_result                           best_cnst = (ordered.begin()->second)->get_result();
  
    {
      vector<bool> active = (ordered.begin()->second)->get_active_data();
      /* in the current case, the indexing of active_data corresponds to the timeslice number */
      while( !active[tmin_cnst] ){ tmin_cnst++; }
    }
  
    A_cnst = (best_cnst.par_values.find("A"))->second;
  }
  else{
    tmin_cnst = control.tmin;
  }
  //************************************************************************************************
  /* perform constant and one exp fits */
  if(!control.only_cnst)
  {
    {
      map<string, param_value> start_params;
      {
        A_cnst.error *= 5.0; /* boost the error on the mass */
        start_params.insert( make_pair("A", A_cnst) );
        
        param_value m(2.0,2.0);
        m.minos = control.minos;
        start_params.insert( make_pair("m", m ) );

        param_value B(0.4,0.4);
        B.minos = control.minos;
        start_params.insert( make_pair("B", B ) );
        
      }

    
      vector<bool> previous( accepted_data.size(), false );
      for(int tlow = control.tmax - control.Nt_min; tlow >= control.tmin; tlow--)
      {
        if(tlow == t0){continue;}
        
        vector<bool> active_data;
        {
          Abscissa* x_tlow = new IntAbscissa(tlow - 1);
          active_data = data_above_x( data, x_tlow );
          delete x_tlow;
        }
        active_data = active_data && accepted_data;
        
        if( (count_active(active_data) >= control.Nt_min) && (active_data != previous) )
        {
          stringstream name; name << "cnst_one_exp_tmin" << tlow;
    
        AvgFit* this_fit = new AvgFit( data, active_data, cnst_one_exp, start_params, minuit_controls, control.correlated, name.str() );
        fits.push_back( this_fit );
        previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** 1 EXP , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
        
      } // next tlow
    } //end one exp fits
    //************************************************************************************************
  
    /* pick the best cnst_one_exp fit to choose a start mass for cnst_two_exp fits */

    int tmin_cnst_one_exp = 0;

    param_value m_cnst_one_exp(2.0,2.0); m_cnst_one_exp.minos = control.minos;
    param_value B_cnst_one_exp(0.4,0.4); B_cnst_one_exp.minos = control.minos;
    param_value A_cnst_one_exp(control.A_start, control.A_err); A_cnst_one_exp.minos = control.minos;


  
    if( fits.size() > 0 )
    {
      map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
      minuit_fit_result                           best_cnst_one_exp = (ordered.begin()->second)->get_result();
    
      {
        vector<bool> active = (ordered.begin()->second)->get_active_data();
        /* in the current case, the indexing of active_data corresponds to the timeslice number */
        while( !active[tmin_cnst_one_exp] ){ tmin_cnst_one_exp++; }
      }
    
      A_cnst_one_exp = (best_cnst_one_exp.par_values.find("A"))->second;
      m_cnst_one_exp = (best_cnst_one_exp.par_values.find("m"))->second;
      B_cnst_one_exp = (best_cnst_one_exp.par_values.find("B"))->second;
    }
    else{
      tmin_cnst_one_exp = control.tmin;
    }

  
    //************************************************************************************************
    /* perform two exp fits */
    if(!control.only_one_exp)
    {
      {
        map<string, param_value> start_params;
        {
          m_cnst_one_exp.error *= 5.0; /* boost the error on the mass */
          start_params.insert( make_pair("m", m_cnst_one_exp) );

          A_cnst_one_exp.error *= 5.0; /* boost the error on the cnst */
          start_params.insert( make_pair("A", A_cnst_one_exp) );
          
          B_cnst_one_exp.error *= 5.0; /* boost the error on the B */
          start_params.insert( make_pair("B", B_cnst_one_exp) );


          param_value dm = m_cnst_one_exp;
          dm.value *= 2.0; dm.error *= 5.0;
          dm.minos = false;
          dm.low_limited = true; dm.low_limit = 0.0;
          start_params.insert( make_pair("dm", dm ) );
          
          param_value C(0.4,0.4);
          start_params.insert( make_pair("C", C) ); 
        }
        
        vector<bool> previous( accepted_data.size(), false );   
        for(int tlow = tmin_cnst_one_exp - 1; tlow >= control.tmin; tlow--){
          if(tlow == t0){continue;}
        
          vector<bool> active_data;
          {
            Abscissa* x_tlow = new IntAbscissa(tlow - 1);
            active_data = data_above_x( data, x_tlow );
            delete x_tlow;
          }
          active_data = active_data && accepted_data;
        
          if( (count_active(active_data) >= control.Nt_min) && (active_data != previous) )
          {	
            stringstream name; name << "cnst_two_exp_tmin" << tlow;
    
            AvgFit* this_fit = new AvgFit( data, active_data, cnst_two_exp, start_params, minuit_controls, control.correlated, name.str() );
            fits.push_back( this_fit );
            previous = active_data;
    
            /*cout << "==================================================================================" << endl;
            cout << "*** 2 EXP , tmin = " << tlow  << " ***" << endl;
            cout << data.print_data(active_data) << endl;
            cout << this_fit->report() << endl;
            cout << "==================================================================================" << endl;*/
          }
        }//next tlow
      
      }//end two exp fits
    }//end of !only exp loop
  }//end of !only cnst loop

  //************************************************************************************************

    
    
  //************************************************************************************************
  /* done with avg fitting, select a best fit, do an ensem fit ... */
  FitSelector fit_selector( fits, fit_qual, chisq_ndof_cutoff );
  //************************************************************************************************
  //cout << fit_selector.get_summary();
  

  //************************************************************************************************
  /* build the output object */
  out.success     = fit_selector.ensem_success();
  out.fit_summary = fit_selector.get_summary();
  
  if( !out.success || control.long_log ){
    out.fit_long_log = fit_selector.get_long_log();
  }
  
  if(out.success)
  {
    param_fit_variation A_param = fit_selector.get_param( "A" );
    out.m                       = A_param.ensem;
    out.m_fit_variation         = A_param.report();

    /* make the best fit data description in case you want to look for outliers */
    vector<bool> ensem_fit_active_data = fit_selector.get_ensem_fit()->get_active_data();
    vector<Abscissa*> x;
    vector<double>    y, y_err;
    
    Function* best_fit_fn = (fit_selector.get_ensem_fit())->get_function();  
    ZtimesliceFunction* ft = static_cast<ZtimesliceFunction*>(best_fit_fn);
      
    map<string, ensem_param_value> ensem_pars = (fit_selector.get_ensem_fit())->get_ensem_param_values();
    
    {
      int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
      vector<double> active_timeslices;
      for(int i = 0; i < x.size(); i++){
        int t_i = ( static_cast<IntAbscissa*>(x[i]) )->get_x();
        active_timeslices.push_back( double(t_i) );
      }
      plot_timeslice_function_data t_y_yerr = plot_timeslice_ensem_function(ft, ensem_pars, active_timeslices );
      
      for(int i = 0; i < t_y_yerr.n_points; i++)
      {
        out.best_data_description.insert( make_pair( int( t_y_yerr.central[i].first ),
                make_pair( t_y_yerr.central[i].second, t_y_yerr.upper[i].second - t_y_yerr.central[i].second )
                )
            );
      }
    }
      
    /* make the plot data text if requested */
    if(control.plots)
    {
      stringstream plot;
      vector<AvgFit*> accepted_fits = fit_selector.get_accepted_fits();

      plot << "## tmin= " << control.tmin << " tmax= " << control.tmax << " nfits = " << accepted_fits.size() << endl;

      /* change this to a nice format write */
      {
        char buff[10]; int nn = sprintf(buff, "%7.5f +/- %7.5f", toDouble(mean(out.m)), toDouble(sqrt(variance(out.m))) );
        plot << "## A= " << buff << endl;
      }

      chisq_dof chisq_desc = (fit_selector.get_ensem_fit())->get_chisq();
      plot << "## chisq/ndof= " << chisq_desc.one_line_report() << endl;
      plot << endl << endl;
      
      { /* write active data */
        plot << "# active data" << endl;
        int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
  
        for(int i = 0; i < x.size(); i++)
        {
          int t_i = ( static_cast<IntAbscissa*>(x[i]) )->get_x(); plot << t_i << " "; 
          plot << y[i] * A_param.ensem_mean << " ";
          plot << y_err[i] * A_param.ensem_mean << endl;
        }
          plot << endl << endl;
      }
      
      
      { /* write inactive data */
        plot << "# inactive data" << endl;
        int count = data.get_active_data(!ensem_fit_active_data, x, y, y_err);
  
        for(int i = 0; i < x.size(); i++)
        {
          int t_i = ( static_cast<IntAbscissa*>(x[i]) )->get_x(); plot << t_i << " "; 
          plot << y[i] * A_param.ensem_mean << " ";
          plot << y_err[i] * A_param.ensem_mean << endl;
        }
        plot << endl << endl;
      }
      
      int n_points = (control.tmax - control.tmin)*10; double dt = 0.1;
      vector<double> plot_t; for(int i = 0; i <= n_points; i++){ plot_t.push_back( double(control.tmin) + double(i)*dt ); }
      
      { /* write ensem fit function */
        plot << "# ensem fit" << endl;
        plot_timeslice_function_data plot_fn = plot_timeslice_ensem_function(ft, ensem_pars, plot_t );
  
        for(int i = 0; i < plot_fn.n_points; i++)
        {
          double t = plot_fn.central[i].first;
          plot << t << "  "
          << plot_fn.lower[i].second   * A_param.ensem_mean << "  "
          << plot_fn.central[i].second * A_param.ensem_mean << "  "
          << plot_fn.upper[i].second   * A_param.ensem_mean << endl;
        }
        plot << endl << endl;
      }
      
      { /* write avg fit variation functions */ 
        for(vector<AvgFit*>::iterator fit = accepted_fits.begin(); fit != accepted_fits.end(); fit++)
        {
          plot << "# avg fit: " << (*fit)->get_fit_name() << endl;
    
          Function* fit_fn = (*fit)->get_function();  
          ZtimesliceFunction* ft = static_cast<ZtimesliceFunction*>(fit_fn);
    
          map<string, param_value> pars = ( (*fit)->get_result()).par_values;
    
          vector< pair<double,double> > plot_fn = plot_timeslice_function(ft, pars, plot_t );
    
          for(int i = 0; i < plot_fn.size(); i++){
            double t = plot_fn[i].first;
            plot << t << "  " << plot_fn[i].second * A_param.ensem_mean << endl;
          }
    
          plot << endl << endl;
        }//next fit variation
      }
      
      out.plot_data = plot.str();
    }// end if plots

  }//end if success
  
  /* clean up pointers */
  for(vector<AvgFit*>::iterator it = fits.begin(); it!= fits.end(); it++){ delete *it; }
  delete cnst; delete cnst_one_exp; delete cnst_two_exp;
  
  return out; 

}; 



//************************************************************************
// PLOTTING UTILITIES
//************************************************************************

plot_timeslice_function_data plot_timeslice_ensem_function( ZtimesliceFunction* fn,
							    const map<string, ensem_param_value>& ensem_pars,
							    const vector<double>& t_values )
{
  if( ensem_pars.size() == 0 ){ cerr << "plot_timeslice_ensem_function:: no parameters provided, exiting " << endl; exit(1); }
  int ncfgs = ((ensem_pars.begin())->second).ensem.size();
  
  plot_timeslice_function_data out;
  out.n_points = t_values.size(); out.tmin = t_values[0]; out.tmax = t_values[ t_values.size() - 1]; /* assuming t_values is ordered ! */
  
  vector< map<string, double> > scaled_down_pars(ncfgs);
  for(map<string, ensem_param_value>::const_iterator par = ensem_pars.begin(); par != ensem_pars.end(); par++)
    {
    EnsemReal down = rescaleEnsemDown( (par->second).ensem );
    for(int cfg = 0; cfg < ncfgs; cfg++){
      (scaled_down_pars[cfg]).insert( make_pair( par->first, toDouble( peekEnsem(down, cfg) ) ) );
    }
  }
  
  for(int i = 0; i < t_values.size(); i++){
    double t = t_values[i];
    
    EnsemReal y_ensem; y_ensem.resize(ncfgs);
    for(int cfg = 0; cfg < ncfgs; cfg++){
      double y = (*fn)( t, scaled_down_pars[cfg] ); pokeEnsem(y_ensem, Real(y), cfg);
    }//next cfg
    
    y_ensem = rescaleEnsemUp(y_ensem);
    pair<double,double> yy = mean_err(y_ensem);
    out.central.push_back( make_pair( t , yy.first ) );
    out.upper.push_back(   make_pair( t , yy.first + yy.second ) );
    out.lower.push_back(   make_pair( t , yy.first - yy.second ) );   
  }//next t
  
  return out;
}


vector< pair<double,double> > plot_timeslice_function( ZtimesliceFunction* fn,
						       const map<string, param_value>& pars,
						       const vector<double>& t_values  ){
  vector< pair<double,double> > out;

  map<string, double> central_pars;
  for(map<string, param_value>::const_iterator par = pars.begin(); par != pars.end(); par++){
    central_pars.insert( make_pair( par->first, (par->second).value ) );
  }
  
  for(int i = 0; i < t_values.size(); i++){
    double t = t_values[i];
    double y = (*fn)( t, central_pars );
    out.push_back( make_pair( t , y ) );
  }//next t
  
  return out;
}

