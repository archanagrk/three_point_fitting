#include "old.h"

fit_three_point_output fit_three_point_corr( const Data& data,                        
				    const vector<bool> accepted_data,       
				    fit_three_point_control control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff               
				    )
{
  /* set minuit controls using defaults */
  MinuitControl minuit_controls; minuit_controls.minos = control.minos;
  
  /* hold the various fits */
  vector<AvgFit*> fits;
  vector<AvgFit*> fits_src_exp;
  vector<AvgFit*> fits_snk_exp;

  /* the output object */
  fit_three_point_output out;

  //************************************************************************************************

  //************************************************************************************************
      /*
      THIS IS EXP FITTING WITH CONSTANT FOR FORM FACTORS
      */
  //************************************************************************************************

  Function* cnst = new ThreePointtimesliceCorrNExp( 0, 0 );  
  Function* cnst_src_exp = new ThreePointtimesliceCorrNExp( 1, 0 );  
  Function* cnst_snk_exp = new ThreePointtimesliceCorrNExp( 0, 1 ); 
  Function* cnst_two_exp = new ThreePointtimesliceCorrNExp( 1, 1 );

  int ct;

  /* perform constant fits */
  {

    map<string, param_value> start_params;
  
    param_value F(control.F_start, control.F_err); //the constant that would be the ff. Constraints?
    F.minos = control.minos;
    start_params.insert( make_pair("F", F ) ); 

  
    vector<bool> previous( accepted_data.size(), false );
    ct = 0;
    for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it)
    {
      for(int tlow = control.tmax[ct] - control.Nt_min; tlow >= control.tmin[ct]; tlow--)
      {
        if(tlow == *it){continue;}
      
        vector<bool> active_data;
        {
          Abscissa* x_tlow = new PairIntAbscissa(make_pair(*it, tlow - 1));
          active_data = data_above_x( data, x_tlow );
          delete x_tlow;
        }
        active_data = active_data && accepted_data;
        
        if( (count_active(active_data) >= control.Nt_min) && (active_data != previous) )
        {
          stringstream name; name << "cnst_dt" << *it << "_tmin" << tlow;
    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          fits_src_exp.push_back( this_fit );
          fits_snk_exp.push_back( this_fit );

          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
      } // next tlow
      ct++;
    } // next dt
  } //end constant fits
  //************************************************************************************************

  /* pick the best constant fit to choose a start Form Factor for one_exp fits */

  vector<std::pair<int,int>> tmin_cnst;

  param_value F_cnst(control.F_start, control.F_err); F_cnst.minos = control.minos;

  if( fits.size() > 0 )
  {
    map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
    minuit_fit_result                           best_cnst = (ordered.begin()->second)->get_result();
  
    {
      vector<bool> active = (ordered.begin()->second)->get_active_data();
      /* in the current case, the indexing of active_data corresponds to the count */

      int indx = 0;
      std::vector<int>::iterator it_snk = control.tsnk.begin();
      std::vector<int>::iterator it_src = control.tsrc.begin();

      for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it)
      {
        while(it_snk != control.tsnk.end() && it_src != control.tsrc.end()){

          int count = indx;
          std::pair<int,int> tmp = make_pair(*it, 0);

          while(!active[count]){ tmp.second++; count++;}

          tmin_cnst.push_back (tmp);
          indx += *it_snk - *it_src + 1;
          ++it_snk; ++it_src;
        }
      }
    }
  
    F_cnst = (best_cnst.par_values.find("F"))->second;
  }


  else{
    ct = 0;
    for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it)
    {
      std::pair<int,int> tmp = make_pair(*it, control.tmin[ct]);
      tmin_cnst.push_back (tmp);
      ct++;
    }
  }

  //************************************************************************************************
  /* perform one and two exp fits */
  if(!control.only_cnst)
  {
    //Function* cnst_one_exp = new ThreePointtimesliceCorrNExp( true, 1);  
    // src_exp fits
    {
      map<string, param_value> start_params;
      {
        F_cnst.error *= 5.0; /* boost the error on the mass */   //Constraints?
        start_params.insert( make_pair("F", F_cnst) );
        
        param_value dmi(2.0,2.0);
        dmi.low_limited = true; dmi.low_limit = 0.0;   // Constraints
        dmi.minos = control.minos;
        start_params.insert( make_pair("dmi", dmi ) );

        param_value Fi(0.04,0.4);
        Fi.minos = control.minos;
        start_params.insert( make_pair("Fi", Fi ) );
        
      }

    
      vector<bool> previous( accepted_data.size(), false );
      ct = 0;
      for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it)
      {
        for(int tlow = control.tmax[ct] - control.Nt_min; tlow >= control.tmin[ct]; tlow--)
        {
          if(tlow == *it){continue;}
          
          vector<bool> active_data;
          {
            Abscissa* x_tlow = new PairIntAbscissa(make_pair(*it, tlow - 1));
            active_data = data_above_x( data, x_tlow );
            delete x_tlow;
          }
          active_data = active_data && accepted_data;
          
          if( (count_active(active_data) >= control.Nt_min) && (active_data != previous) )
          {
            stringstream name; name << "cnst_src_exp_dt" << *it << "_tmin" << tlow;
      
            AvgFit* this_fit = new AvgFit( data, active_data, cnst_src_exp, start_params, minuit_controls, control.correlated, name.str() );
            fits.push_back( this_fit );
            fits_src_exp.push_back( this_fit );
            previous = active_data;
      
            /*cout << "==================================================================================" << endl;
            cout << "*** 1 EXP , tmin = " << tlow  << " ***" << endl;
            cout << data.print_data(active_data) << endl;
            cout << this_fit->report() << endl;
            cout << "==================================================================================" << endl; */
          }
        } // next tlow
        ct++;
      } // next dt
    } //end of src_exp fits

    //snk_exp fits
    {
      map<string, param_value> start_params;
      {
        //F_cnst.error *= 5.0; /* boost the error on the mass */    // Constraints?
        start_params.insert( make_pair("F", F_cnst) );
        
        param_value dmf(2.0,2.0);
        dmf.low_limited = true; dmf.low_limit = 0.0;   // Constraints
        dmf.minos = control.minos;
        start_params.insert( make_pair("dmf", dmf ) );

        param_value Ff(0.04,0.4);
        Ff.minos = control.minos;
        start_params.insert( make_pair("Ff", Ff ) );
        
      }
    
      vector<bool> previous( accepted_data.size(), false );
      int count = 0;
      for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it) 
      {
        for(int tlow = control.tmax[ct] - control.Nt_min; tlow >= control.tmin[ct]; tlow--)
        {
          if(tlow == *it){continue;}
          
          vector<bool> active_data;
          {
            Abscissa* x_tlow = new PairIntAbscissa(make_pair(*it, tlow - 1));
            active_data = data_above_x( data, x_tlow );
            delete x_tlow;
          }
          active_data = active_data && accepted_data;
          
          if( (count_active(active_data) >= control.Nt_min) && (active_data != previous) )
          {
            stringstream name; name << "cnst_snk_exp_dt" << *it << "_tmin" << tlow ;
      
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_snk_exp, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          fits_snk_exp.push_back( this_fit );
          previous = active_data;
      
            /*cout << "==================================================================================" << endl;
            cout << "*** 1 EXP , tmin = " << tlow  << " ***" << endl;
            cout << data.print_data(active_data) << endl;
            cout << this_fit->report() << endl;
            cout << "==================================================================================" << endl; */
          }
          
        } // next tlow
        ct++;
      } //next dt
    } //end of snk_exp fits

    //************************************************************************************************
  
    /* pick the best cnst_src_exp fit and cnst_snk_exp to choose a start mass for cnst_two_exp fits */

    vector<std::pair<int,int>>  tmin_cnst_snk_exp;
    //vector<std::pair<int,int>>  tmin_cnst_src_exp = 0;

    param_value F_cnst_one_exp(control.F_start, control.F_err);
    F_cnst_one_exp.minos = control.minos;

    param_value dmi_cnst_src_exp(2.0,2.0); 
    dmi_cnst_src_exp.low_limited = true; dmi_cnst_src_exp.low_limit = 0.0;   // Constraints
    dmi_cnst_src_exp.minos = control.minos;

    param_value Fi_cnst_src_exp(0.04, 0.7);
    Fi_cnst_src_exp.minos = control.minos;

    param_value dmf_cnst_snk_exp(2.0,2.0);
    dmf_cnst_snk_exp.low_limited = true; dmf_cnst_snk_exp.low_limit = 0.0;   // Constraints
    dmf_cnst_snk_exp.minos = control.minos;

    param_value Ff_cnst_snk_exp(0.04, 0.7);
    Ff_cnst_snk_exp.minos = control.minos;    


  
    if( fits.size() > 0 )
    {
      map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
      minuit_fit_result                           best_cnst_one_exp = (ordered.begin()->second)->get_result();

      map<double, AvgFit*, std::greater<double> > ordered_src = make_ordered_list( fits_src_exp, fit_qual );
      minuit_fit_result                           best_cnst_src_exp = (ordered.begin()->second)->get_result();

      map<double, AvgFit*, std::greater<double> > ordered_snk = make_ordered_list( fits_snk_exp, fit_qual );
      minuit_fit_result                           best_cnst_snk_exp = (ordered.begin()->second)->get_result();



      {
      vector<bool> active = (ordered.begin()->second)->get_active_data();
      /* in the current case, the indexing of active_data corresponds to the count */

      int indx = 0;
      std::vector<int>::iterator it_snk = control.tsnk.begin();
      std::vector<int>::iterator it_src = control.tsrc.begin();

      for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it)
      {
        while(it_snk != control.tsnk.end() && it_src != control.tsrc.end()){

          int count = indx;
          std::pair<int,int> tmp = make_pair(*it, 0);

          while(!active[count]){ tmp.second++; count++;}

          tmin_cnst.push_back (tmp);
          indx += *it_snk - *it_src + 1;
          ++it_snk; ++it_src;
        }
      }
    }
    
      F_cnst_one_exp = (best_cnst_one_exp.par_values.find("F"))->second;
      Fi_cnst_src_exp = (best_cnst_src_exp.par_values.find("Fi"))->second;
      dmi_cnst_src_exp = (best_cnst_src_exp.par_values.find("dmi"))->second;
      Ff_cnst_snk_exp = (best_cnst_snk_exp.par_values.find("Ff"))->second;
      dmf_cnst_snk_exp = (best_cnst_snk_exp.par_values.find("dmf"))->second;      
    }

    else{
      ct = 0;
      for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it)
      {
        std::pair<int,int> tmp = make_pair(*it, control.tmin[ct]);
        tmin_cnst_snk_exp.push_back (tmp);
        ct++;
      }
    }

  
    //************************************************************************************************
    /* perform two exp fits */
    if(!control.only_one_exp)
    {
      {
        map<string, param_value> start_params;
        {
          F_cnst_one_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair("F", F_cnst_one_exp) );

          Fi_cnst_src_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair("Fi", Fi_cnst_src_exp) );

          Ff_cnst_snk_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair("Ff", Ff_cnst_snk_exp) );

          dmi_cnst_src_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair("dmi", dmi_cnst_src_exp) );

          dmf_cnst_snk_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair("dmf", dmf_cnst_snk_exp) );

        }
        
        vector<bool> previous( accepted_data.size(), false );  
        ct = 0; 
        for(std::vector<int>::iterator it = control.dt.begin() ; it != control.dt.end(); ++it)
        {
          for(int tlow = control.tmax[ct] - control.Nt_min; tlow >= control.tmin[ct]; tlow--)
          {
            if(tlow == *it){continue;}
          
            vector<bool> active_data;
            {
              Abscissa* x_tlow = new PairIntAbscissa(make_pair(*it, tlow - 1));
              active_data = data_above_x( data, x_tlow );
              delete x_tlow;
            }
            active_data = active_data && accepted_data;
          
            if( (count_active(active_data) >= control.Nt_min) && (active_data != previous) )
            {	
              stringstream name; name << "cnst_two_exp_dt" << *it  << "_tmin" << tlow ;
      
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
          ct++;
        }// next dt
      
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
    param_fit_variation F_param = fit_selector.get_param( "F" );
    out.F                       = F_param.ensem;
    out.F_fit_variation         = F_param.report();

    /* make the best fit data description in case you want to look for outliers */
    vector<bool> ensem_fit_active_data = fit_selector.get_ensem_fit()->get_active_data();
    vector<Abscissa*> x;
    vector<double>    y, y_err;
    
    Function* best_fit_fn = (fit_selector.get_ensem_fit())->get_function();  
    ThreePointtimesliceCorrFunction* ft = static_cast<ThreePointtimesliceCorrFunction*>(best_fit_fn);
      
    map<string, ensem_param_value> ensem_pars = (fit_selector.get_ensem_fit())->get_ensem_param_values();
    
    
    int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
    vector<pair<double,double>> active_timeslices;
    for(int i = 0; i < x.size(); i++)
    {
      pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x();
      active_timeslices.push_back( make_pair(double(t_i.first),double(t_i.second)) );
    }

    std::vector<double> delt(control.dt.begin(), control.dt.end());

    plot_three_point_timeslice_function_data t_y_yerr = plot_three_point_timeslice_ensem_function(ft, ensem_pars, active_timeslices, delt );

    for(int i = 0; i < t_y_yerr.central.size(); i++)
    {
      for(int j = 0; j < t_y_yerr.dt.size(); j++ )
      {
        int k =0;
        while(k < t_y_yerr.n_points[j]){
          out.best_data_description.insert( make_pair( 
                (t_y_yerr.central[i].first),
                make_pair( t_y_yerr.central[i].second, t_y_yerr.upper[i].second - t_y_yerr.central[i].second )
                )
            );
            k++;}

      }
      
      /* make the plot data text if requested */
      if(control.plots)
      {
        stringstream plot;
        vector<AvgFit*> accepted_fits = fit_selector.get_accepted_fits();

        for(int j = 0; j < control.dt.size(); j++){
          plot << "## tmin= " << control.tmin[j] << " tmax= " << control.tmax[j] << " nfits = " << accepted_fits.size() << endl;
        }
        

        /* change this to a nice format write */
        {
          char buff[10]; int nn = sprintf(buff, "%7.5f +/- %7.5f", toDouble(mean(out.F)), toDouble(sqrt(variance(out.F))) );
          plot << "## F= " << buff << endl;
        }

        chisq_dof chisq_desc = (fit_selector.get_ensem_fit())->get_chisq();
        plot << "## chisq/ndof= " << chisq_desc.one_line_report() << endl;
        plot << endl << endl;
        
        { /* write active data */
          plot << "# active data" << endl;
          int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
    
          for(int i = 0; i < x.size(); i++)
          {
            pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x(); plot <<  t_i.first << t_i.second  << " "; 
            plot << y[i] * F_param.ensem_mean << " ";
            plot << y_err[i] * F_param.ensem_mean << endl;
          }
            plot << endl << endl;
        }
        
        
        { /* write inactive data */
          plot << "# inactive data" << endl;
          int count = data.get_active_data(!ensem_fit_active_data, x, y, y_err);
    
          for(int i = 0; i < x.size(); i++)
          {
            pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x(); plot <<  t_i.first << t_i.second  << " ";  
            plot << y[i] * F_param.ensem_mean << " ";
            plot << y_err[i] * F_param.ensem_mean << endl;
          }
          plot << endl << endl;
        }

        vector<pair<double,double>> plot_t; 
        
        for(int j = 0; j < control.dt.size(); j++){
          int n_points = ((control.tsnk[j] - control.tsrc[j]))*10 + 1; double dt = 0.1;
          for(int i = 0; i <= n_points; i++){ plot_t.push_back( make_pair(double(control.dt[j]),double(control.tsrc[j]) + double(i)*dt) ); }
        }
        
        { /* write ensem fit function */
          plot << "# ensem fit" << endl;
          plot_three_point_timeslice_function_data plot_fn = plot_three_point_timeslice_ensem_function(ft, ensem_pars, plot_t, delt );

          for(int i = 0; i < plot_fn.central.size(); i++)
          {
            for(int j = 0; j < plot_fn.dt.size(); j++ )
              {
              int k =0;
              while(k < plot_fn.n_points[j]){
                pair<double,double> t = plot_fn.central[i].first;
                plot << t.second << "  "
                << plot_fn.lower[i].second   * F_param.ensem_mean << "  "
                << plot_fn.central[i].second * F_param.ensem_mean << "  "
                << plot_fn.upper[i].second   * F_param.ensem_mean << endl;
              }
            }
            plot << endl << endl;
          }
        }
        
        { /* write avg fit variation functions */ 
          for(vector<AvgFit*>::iterator fit = accepted_fits.begin(); fit != accepted_fits.end(); fit++)
          {
            plot << "# avg fit: " << (*fit)->get_fit_name() << endl;
      
            Function* fit_fn = (*fit)->get_function();  
            ThreePointtimesliceCorrFunction* ft = static_cast<ThreePointtimesliceCorrFunction*>(fit_fn);
      
            map<string, param_value> pars = ( (*fit)->get_result()).par_values;
      
            vector< pair<pair<double,double>,double> >  plot_fn = plot_three_point_timeslice_function(ft, pars, plot_t);
      
            for(int i = 0; i < plot_fn.size(); i++){
              pair<double,double> t = plot_fn[i].first;
              plot << t.second << "  " << plot_fn[i].second * F_param.ensem_mean << endl;
            }
      
            plot << endl << endl;
          }//next fit variation
        }
        
        out.plot_data = plot.str();
      }// end if plots

    }//end if success
    
    /* clean up pointers */
    for(vector<AvgFit*>::iterator it = fits.begin(); it!= fits.end(); it++){ delete *it; }
    delete cnst; 
    delete cnst_src_exp; delete cnst_snk_exp; delete cnst_two_exp;
    

    return out; 

 }; 
}





//************************************************************************
// PLOTTING UTILITIES
//************************************************************************

plot_three_point_timeslice_function_data plot_three_point_timeslice_ensem_function( ThreePointtimesliceCorrFunction* fn,
							    const map<string, ensem_param_value>& ensem_pars,
							    const vector<pair<double,double>>& t_values, const vector<double>& delt )
{
  if( ensem_pars.size() == 0 ){ cerr << "plot_three_point_timeslice_ensem_function:: no parameters provided, exiting " << endl; exit(1); }
  int ncfgs = ((ensem_pars.begin())->second).ensem.size();
  
  plot_three_point_timeslice_function_data out;
  out.dt = delt;

  for(auto it = delt.begin() ; it != delt.end(); ++it)
  {
    int count = 0;
    int tmin = 0;
    int tmax = 0;

    for(auto it1 = t_values.begin() ; it1 != t_values.end(); ++it1)
    {
      if(*it == it1->first){ count++; if(it1->second < tmin){tmin = it1->second;} if(it1->second > tmax) {tmax = it1->second;}}
      else{continue;}
    }
      out.n_points.push_back(count);
      out.tmin.push_back(tmin);
      out.tmax.push_back(tmin);
      out.dt.push_back(*it);

  } /* assuming t_values is ordered ! */
  
  vector< map<string, double> > scaled_down_pars(ncfgs);
  for(map<string, ensem_param_value>::const_iterator par = ensem_pars.begin(); par != ensem_pars.end(); par++)
    {
    EnsemReal down = rescaleEnsemDown( (par->second).ensem );
    for(int cfg = 0; cfg < ncfgs; cfg++){
      (scaled_down_pars[cfg]).insert( make_pair( par->first, toDouble( peekEnsem(down, cfg) ) ) );
    }
  }
  
  for(int i = 0; i < t_values.size(); i++){
    pair<double,double> t = make_pair(t_values[i].first,t_values[i].second);
    
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


vector< pair<pair<double,double>,double> > plot_three_point_timeslice_function( ThreePointtimesliceCorrFunction* fn,
						       const map<string, param_value>& pars,
						       const vector<std::pair<double,double>>& t_values){
  vector< pair<pair<double,double>,double> > out;

  map<string, double> central_pars;
  for(map<string, param_value>::const_iterator par = pars.begin(); par != pars.end(); par++){
    central_pars.insert( make_pair( par->first, (par->second).value ) );
  }
  
  for(int i = 0; i < t_values.size(); i++){
    pair<double,double> t = make_pair(t_values[i].first,t_values[i].second);
    double y = (*fn)( t, central_pars );
    out.push_back( make_pair( t , y ) );
  }//next t
  
  return out;
}

