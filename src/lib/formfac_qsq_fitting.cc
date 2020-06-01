#include "formfac_qsq_fitting.h"

fit_ff_qsq_output fit_ff_qsq( const Data& data,                        
				    const vector<bool> accepted_data,       
				    fit_ff_qsq_control control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff               
				    )
{
  /* set minuit controls using defaults */
  MinuitControl minuit_controls; minuit_controls.minos = control.minos;
  
  /* hold the various fits */
  vector<AvgFit*> fits;

  /* the output object */
  fit_ff_qsq_output out;

  //************************************************************************************************

  //************************************************************************************************
      /*
      THIS IS FITTING WITH VMD
      */
  //************************************************************************************************

  Function* vmd = new FFqsqVMD( control.m_p );
  Function* gauss = new FFqsqGauss( control.m_p );  

  control.msq_start = control.m_p;
  control.msq_err   = data.get_all_y_err()[0];

  control.f0_start = data.get_all_y_mean()[0];
  control.f0_err   = data.get_all_y_err()[0];


  if(control.fit_form == "ff_qsq_vmd" || control.fit_form == "all" )
  { 
    /* perform vmd fit */
    {

      map<string, param_value> start_params;
    
      param_value msq(control.msq_start, control.msq_err); 
      param_value f0(control.f0_start, control.f0_err);       
      msq.minos = control.minos; f0.minos = control.minos;  

      start_params.insert( make_pair("msq", msq ) ); 
      start_params.insert( make_pair("f0", f0 ) );       

    
      vector<bool> previous( accepted_data.size(), false );
      for(double qmax = control.qsqmax; qmax >= control.qsqmin; qmax -= 0.005){

        vector<bool> active_data;
        {
          Abscissa* x_qhigh = new DoubleAbscissa(qmax+0.00001);
          active_data = data_below_x( data, x_qhigh );
          delete x_qhigh;
        }
        active_data = active_data && accepted_data;  
        
        if( (count_active(active_data) >= control.Nq_min) && (active_data != previous) )
        {
          stringstream name; name << "vmd_qmin" << control.qsqmin << "_qmax" << qmax;
    
          AvgFit* this_fit = new AvgFit( data, active_data, vmd, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** VMD , qsqmin = " << qlow  << " ***" << endl;
          cout << data.prdouble_data(accepted_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
        
      }//next qmax
        
    } //end vmd fit
    //************************************************************************************************
  }

  //************************************************************************************************
      /*
      THIS IS FITTING WITH GAUSSIAN
      */
  //************************************************************************************************
  if(control.fit_form == "ff_qsq_gauss" || control.fit_form == "all" )
  { 

    /* perform gauss fit */
    {

      map<string, param_value> start_params;
    
      param_value betasq(control.msq_start, control.msq_err); 
      param_value f0(control.f0_start, control.f0_err);       
      betasq.minos = control.minos; f0.minos = control.minos;  

      start_params.insert( make_pair("betasq", betasq ) ); 
      start_params.insert( make_pair("f0", f0 ) );   

    
      vector<bool> previous( accepted_data.size(), false );
      for(double qmax = control.qsqmax; qmax >= control.qsqmin; qmax -= 0.005){

        vector<bool> active_data;
        {
          Abscissa* x_qhigh = new DoubleAbscissa(qmax+0.00001);
          active_data = data_below_x( data, x_qhigh );
          delete x_qhigh;
        }
        active_data = active_data && accepted_data; 
    
      
        if( (count_active(active_data) >= control.Nq_min) && (active_data != previous) )
        {
          stringstream name; name << "gauss_qmin" << control.qsqmin << "_qmax" << qmax;
    
          AvgFit* this_fit = new AvgFit( data, active_data, gauss, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          previous = active_data;

          /*cout << "==================================================================================" << endl;
          cout << "*** GAUSS , qsqmin = " << qlow  << " ***" << endl;
          cout << data.prdouble_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
      }//next qmax        
    } //end gauss fit
    //************************************************************************************************

  }

  //************************************************************************************************
      /*
      THIS IS FITTING WITH Z EXPANSION
      */
  //************************************************************************************************
  if(control.fit_form == "ff_qsq_z_exp" || control.fit_form == "all" )
  { 

    /* perform z exp fit */
    {

      map<string, param_value> start_params;
    
      // param_value msq(control.msq_start, control.msq_err); 
      param_value f0(control.f0_start, control.f0_err); 

      f0.low_limited = true; f0.low_limit = -control.fbound;
      f0.high_limited = true; f0.high_limit = control.fbound;   
      // msq.minos = control.minos; f0.minos = control.minos;  f0.fixed = false;

      // start_params.insert( make_pair("msq", msq ) ); 
      start_params.insert( make_pair("f0", f0 ) );   

    
      vector<bool> previous( accepted_data.size(), false );
      for(double qmax = control.qsqmax; qmax >= control.qsqmin; qmax -= 0.005){

        vector<bool> active_data;
        {
          Abscissa* x_qhigh = new DoubleAbscissa(qmax+0.00001);
          active_data = data_below_x( data, x_qhigh );
          delete x_qhigh;
        }
        active_data = active_data && accepted_data; 
    
      
        if( (count_active(active_data) >= control.Nq_min) && (active_data != previous) )
        {
          for(int n_max = 1; n_max <= control.n_max; n_max++){
            stringstream name; name << "z_exp" << n_max << "_qmin" << control.qsqmin << "_qmax" << qmax;

            map<string, param_value> start_params_nmax = start_params;

            for(int pars = 1; pars < n_max; pars++){
              stringstream f; f << "f" << pars; 
              start_params_nmax.insert( make_pair(f.str(), f0 ) );   
            }
          
            Function* z_exp = new FFqsqZExp( control.m_p, qmax, control.topt, control.tcut, n_max );   
      
            AvgFit* this_fit = new AvgFit( data, active_data, z_exp, start_params_nmax, minuit_controls, control.correlated, name.str() );
            fits.push_back( this_fit );
          }
          previous = active_data;

          /*cout << "==================================================================================" << endl;
          cout << "*** GAUSS , qsqmin = " << qlow  << " ***" << endl;
          cout << data.prdouble_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
      }//next qmax        
    } //end z exp fit
    //************************************************************************************************

  }

  //************************************************************************************************
  /* done with avg fitting, select a best fit, do an ensem fit ... */
  FitSelector fit_selector( fits, fit_qual, chisq_ndof_cutoff );
  //************************************************************************************************
  //cout << fit_selector.get_summary();
  

  //************************************************************************************************
  /* build the output object */
  out.success     = fit_selector.ensem_success();
  
  if( !out.success || control.long_log ){
    out.fit_long_log = fit_selector.get_long_log();
    cout << fit_selector.get_summary();
  }

  out.fit_summary = fit_selector.get_summary();

  if(out.success)
  {

    param_fit_variation msq_param;
    string fit_name  = fit_selector.get_ensem_fit()->get_fit_name();

    // if( (fit_name.find("vmd")!=std::string::npos) || (fit_name.find("z_exp")!=std::string::npos) ){msq_param = fit_selector.get_param( "msq" );}
    if( (fit_name.find("vmd")!=std::string::npos)){msq_param = fit_selector.get_param( "msq" );}
    if(fit_name.find("gauss")!=std::string::npos){msq_param = fit_selector.get_param( "betasq" );}
    if(fit_name.find("z_exp")!=std::string::npos){msq_param = fit_selector.get_param( "f0" );}

    out.msq                       = msq_param.ensem;
    out.msq_fit_variation         = msq_param.report();


    /* make the best fit data description in case you want to look for outliers */
    vector<bool> ensem_fit_active_data = fit_selector.get_ensem_fit()->get_active_data();
    vector<Abscissa*> x;
    vector<double>    y, y_err;
    
    Function* best_fit_fn = (fit_selector.get_ensem_fit())->get_function();  
    FFqsqFunction* ft = static_cast<FFqsqFunction*>(best_fit_fn);

    map<string, ensem_param_value> ensem_pars = (fit_selector.get_ensem_fit())->get_ensem_param_values();

    /* get the radius */
    EnsemReal radius; radius.resize(((ensem_pars.begin())->second).ensem.size());
    // EnsemReal norm; norm.resize(((ensem_pars.begin())->second).ensem.size());


    if( (fit_name.find("vmd")!=std::string::npos) ){ 
      //radius = rescaleEnsemUp(sqrt(toScalar(6.0)*rescaleEnsemDown(ensem_pars.find("f0")->second.ensem)/rescaleEnsemDown(ensem_pars.find("msq")->second.ensem)));}
      radius = rescaleEnsemUp(sqrt(toScalar(6.0)/rescaleEnsemDown(ensem_pars.find("msq")->second.ensem)));}  

    if(fit_name.find("gauss")!=std::string::npos){
      radius = rescaleEnsemUp(sqrt((toScalar(3.0))/(toScalar(8.0)*rescaleEnsemDown(ensem_pars.find("betasq")->second.ensem)))) ;}
    
    if(fit_name.find("z_exp")!=std::string::npos){

      FFqsqZExp* fn = static_cast<FFqsqZExp*>(best_fit_fn);

      double tcut = control.tcut;

      int n_max = fn->get_nmax();

      double qsqmax = fn->get_qsqmax();
      double t0;

      if(control.topt){t0 = tcut * (1 - sqrt(1 + (control.qsqmax/tcut)) );}
      else{t0 = 0;}


      double dz = ((sqrt(tcut) - sqrt(tcut -t0))/(2*sqrt(tcut)* pow(sqrt(tcut) + sqrt(tcut - t0),2) ) ) + (1/(2*sqrt(tcut)*(sqrt(tcut) + sqrt(tcut - t0))));
      double z  = (sqrt(tcut) - sqrt(tcut -t0)) / (sqrt(tcut) + sqrt(tcut -t0));

      // EnsemReal dcoeff = toScalar(6.0)/rescaleEnsemDown(ensem_pars.find("msq")->second.ensem);
      EnsemReal f0 = rescaleEnsemDown(ensem_pars.find("f0")->second.ensem);

      // radius = (dcoeff * f0);
      EnsemReal norm = f0;

      for(int num = 1; num < n_max; num++){
        stringstream f; f << "f" << num;
        EnsemReal ff = rescaleEnsemDown(ensem_pars.find(f.str())->second.ensem);
        // radius += (ff * toScalar(-6.0 * num * pow(z,num-1) * dz)) + (dcoeff * ff * toScalar(pow(z,num)));
        radius += (ff * toScalar(-6.0 * num * pow(z,num-1) * dz));
        norm += ff * toScalar(pow(z,num));
      }

      radius = rescaleEnsemUp(sqrt(radius/norm));
    }



    {
      int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
      vector<double> active_qsq;
      for(int i = 0; i < x.size(); i++){
        double qsq_i = ( static_cast<DoubleAbscissa*>(x[i]) )->get_x();
        active_qsq.push_back( qsq_i );
      }
      plot_qsq_function_data x_y_yerr = plot_qsq_ensem_function(ft, ensem_pars, active_qsq );
      
      for(int i = 0; i < x_y_yerr.n_points; i++)
      {
        out.best_data_description.insert( make_pair( x_y_yerr.central[i].first ,
                make_pair( x_y_yerr.central[i].second, x_y_yerr.upper[i].second - x_y_yerr.central[i].second )
                )
            );
      }
    }
      
    /* make the plot data text if requested */
    if(control.plots)
    {
      stringstream plot;
      vector<AvgFit*> accepted_fits = fit_selector.get_accepted_fits();
      /* radius */
      {
        char buff[10]; int nn = sprintf(buff, "%7.5f +/- %7.5f", toDouble(mean(radius)), toDouble(sqrt(variance(radius))) );
        plot << "## radius= " << buff << endl;
      }

      plot << "## qsqmin= " << control.qsqmin << " qsqmax= " << control.qsqmax << " nfits = " << accepted_fits.size() << endl;

      /* change this to a nice format write */
      {
        char buff[10]; int nn = sprintf(buff, "%7.5f +/- %7.5f", toDouble(mean(out.msq)), toDouble(sqrt(variance(out.msq))) );
        plot << "## msq/betasq= " << buff << endl;
      }

      chisq_dof chisq_desc = (fit_selector.get_ensem_fit())->get_chisq();
      plot << "## chisq/ndof= " << chisq_desc.one_line_report() << endl;
      plot << endl << endl;
      
      { /* write active data */
        plot << "# active data" << endl;
        int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
  
        for(int i = 0; i < x.size(); i++)
        {
          double qsq_i = ( static_cast<DoubleAbscissa*>(x[i]) )->get_x(); plot << qsq_i << " "; 
          plot << y[i] << " ";
          plot << y_err[i] << endl;
        }
          plot << endl << endl;
      }
      
      
      { /* write inactive data */
        plot << "# inactive data" << endl;
        int count = data.get_active_data(!ensem_fit_active_data, x, y, y_err);
  
        for(int i = 0; i < x.size(); i++)
        {
          double qsq_i = ( static_cast<DoubleAbscissa*>(x[i]) )->get_x(); plot << qsq_i << " "; 
          plot << y[i] << " ";
          plot << y_err[i] << endl;
        }
        plot << endl << endl;
      }
      
      int n_points = (control.qsqmax - control.qsqmin)*100000; double dq = 0.00001;
      vector<double> plot_t; for(int i = 0; i <= n_points; i++){ plot_t.push_back( double(control.qsqmin) + double(i)*dq ); }
      
      { /* write ensem fit function */
        plot << "# ensem fit" << endl;
        plot_qsq_function_data plot_fn = plot_qsq_ensem_function(ft, ensem_pars, plot_t );
  
        for(int i = 0; i < plot_fn.n_points; i++)
        {
          double qsq = plot_fn.central[i].first;
          plot << qsq << "  "
          << plot_fn.lower[i].second   << "  "
          << plot_fn.central[i].second << "  "
          << plot_fn.upper[i].second   << endl;
        }
        plot << endl << endl;
      }
      
      { /* write avg fit variation functions */ 
        for(vector<AvgFit*>::iterator fit = accepted_fits.begin(); fit != accepted_fits.end(); fit++)
        {
          plot << "# avg fit: " << (*fit)->get_fit_name() << endl;
    
          Function* fit_fn = (*fit)->get_function();  
          FFqsqFunction* ft = static_cast<FFqsqFunction*>(fit_fn);
    
          map<string, param_value> pars = ( (*fit)->get_result()).par_values;
    
          vector< pair<double,double> > plot_fn = plot_qsq_function(ft, pars, plot_t );
    
          for(int i = 0; i < plot_fn.size(); i++){
            double qsq = plot_fn[i].first;
            plot << qsq << "  " << plot_fn[i].second << endl;
          }
    
          plot << endl << endl;
        }//next fit variation
      }
      
      out.plot_data = plot.str();
    }// end if plots

  }//end if success
  
  /* clean up pointers */
  for(vector<AvgFit*>::iterator it = fits.begin(); it!= fits.end(); it++){ delete *it; }

  delete vmd; delete gauss;
  
  return out; 

}; 



//************************************************************************
// PLOTTING UTILITIES
//************************************************************************

plot_qsq_function_data plot_qsq_ensem_function( FFqsqFunction* fn,
							    const map<string, ensem_param_value>& ensem_pars,
							    const vector<double>& qsq_values )
{
  if( ensem_pars.size() == 0 ){ cerr << "plot_qsq_ensem_function:: no parameters provided, exiting " << endl; exit(1); }
  int ncfgs = ((ensem_pars.begin())->second).ensem.size();
  
  plot_qsq_function_data out;
  out.n_points = qsq_values.size(); out.qsqmin = qsq_values[0]; out.qsqmax = qsq_values[ qsq_values.size() - 1]; /* assuming qsq_values is ordered ! */
  
  vector< map<string, double> > scaled_down_pars(ncfgs);
  for(map<string, ensem_param_value>::const_iterator par = ensem_pars.begin(); par != ensem_pars.end(); par++)
    {
    EnsemReal down = rescaleEnsemDown( (par->second).ensem );
    for(int cfg = 0; cfg < ncfgs; cfg++){
      (scaled_down_pars[cfg]).insert( make_pair( par->first, toDouble( peekEnsem(down, cfg) ) ) );
    }
  }
  
  for(int i = 0; i < qsq_values.size(); i++){
    double t = qsq_values[i];
    
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


vector< pair<double,double> > plot_qsq_function( FFqsqFunction* fn,
						       const map<string, param_value>& pars,
						       const vector<double>& qsq_values  ){
  vector< pair<double,double> > out;

  map<string, double> central_pars;
  for(map<string, param_value>::const_iterator par = pars.begin(); par != pars.end(); par++){
    central_pars.insert( make_pair( par->first, (par->second).value ) );
  }
  
  for(int i = 0; i < qsq_values.size(); i++){
    double q = qsq_values[i];
    double f = (*fn)( q, central_pars );
    out.push_back( make_pair( q , f ) );
  }//next t
  
  return out;
}

