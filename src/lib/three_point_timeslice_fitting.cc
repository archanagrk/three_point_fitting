#include "three_point_timeslice_fitting.h"

//************************************************************************
// DRIVER OVER DIFFERENT SOURCE-SINK SEPERATIONS
//************************************************************************


fit_three_point_output fit( const vector<Data>& data,                        
				    const vector<vector<bool>> accepted_data,       
				    vector<fit_three_point_control> control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff, double num         
				    )

{

  fit_three_point_output out;
  vector<single_dt_avg_fit_output> single_dt_fits;

  for(int i = 0; i < num; i++){ 
    single_dt_fits.push_back( single_dt_fit(data[i], accepted_data[i], control[i], fit_qual, chisq_ndof_cutoff) );
    }
  out = fit_three_point_corr(data[num], accepted_data[num], control[num], fit_qual, single_dt_fits, chisq_ndof_cutoff);
 

  return out;

};


//************************************************************************
// FIX THE RANGES FOR EACH DT
//************************************************************************


single_dt_avg_fit_output single_dt_fit( const Data& data,                        
				    const vector<bool> accepted_data,       
				    fit_three_point_control& control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff             
				    )
{
  /* set minuit controls using defaults */
  MinuitControl minuit_controls; minuit_controls.minos = control.minos; 
  minuit_controls.limit_max_fcn_calls = control.limit_max_fcn_calls;
  minuit_controls.max_fcn_calls = control.max_fcn_calls; 

  /* the output object */
  single_dt_avg_fit_output out;

  
  /* hold the various fits */
  vector<AvgFit*> fits;
  vector<AvgFit*> fits_src_exp;
  vector<AvgFit*> fits_snk_exp;

  //************************************************************************************************

  //************************************************************************************************
      /*
      THIS IS EXP FITTING WITH CONSTANT FOR FORM FACTORS
      */
  //************************************************************************************************



  Function* cnst         = new ThreePointtimesliceCorrNExp( {0}, {0}, control.dt );  
  Function* cnst_src_exp = new ThreePointtimesliceCorrNExp( {1}, {0}, control.dt );  
  Function* cnst_snk_exp = new ThreePointtimesliceCorrNExp( {0}, {1}, control.dt ); 
  Function* cnst_two_exp = new ThreePointtimesliceCorrNExp( {1}, {1}, control.dt );


  /* perform constant fits */
  {

    map<string, param_value> start_params;
    vector<PairIntAbscissa*> x_tlow;
    vector<PairIntAbscissa*> x_thigh;

    control.F_start = data.get_all_y_mean()[ int(control.dt[0]/2) ];
    control.F_err   = data.get_all_y_err()[ int(control.dt[0]/2) ]*5;
  
    param_value F(control.F_start, control.F_err); //the constant that would be the ff. Constraints?
    F.minos = control.minos;
    F.fixed = false;
    start_params.insert( make_pair("F", F ) ); 


    std::get<1>(control.tmin_max[0]) = int(std::get<0>(control.tmin_max[0])/2 - control.Nt_min[0]/2);
    std::get<2>(control.tmin_max[0]) = int(std::get<0>(control.tmin_max[0])/2 + control.Nt_min[0]/2);


  
    vector<bool> previous( accepted_data.size(), false ); 

    std::vector<pair<pair<int,int>,pair<int,int>>> t_pairs  = get_single_dt_t_ranges(control,1,1,1);

    cout << t_pairs.size() << endl;

    //Looping over all t ranges
    cout << "-------------------------------------------------" << endl;   
    for(auto t_pair = t_pairs.begin(); t_pair != t_pairs.end(); t_pair++){

      PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pair->first);
      PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pair->second);


      vector<bool> active_data = !(data.make_active_data()); //all false
        
      cout  << "cnst__t_low  =  "<< *x_t_low_dt << " | ";
      cout << "cnst__t_high  =  " << *x_t_high_dt << endl;


      Abscissa* xl = new PairIntAbscissa(make_pair(control.dt[0], t_pair->first.second - 1));
      Abscissa* xh = new PairIntAbscissa(make_pair(control.dt[0], t_pair->second.second + 1));
      
      active_data = active_data || data_in_x_range( data, make_pair(xl,xh) );
      active_data = active_data && accepted_data;
      delete xl; delete xh;

      if( (count_active(active_data) >= control.Nt_min[0]) && (active_data != previous) )
      {
        
        stringstream name; name << "c";

        name << "_Dt" << x_t_low_dt->get_x().first << "_" << x_t_low_dt->get_x().second << "-" << x_t_high_dt->get_x().second;

        AvgFit* this_fit = new AvgFit( data, active_data, cnst, start_params, minuit_controls, control.correlated, name.str() );
        if(this_fit->is_valid())
        {
          cout << "        " << "SUCCESSFUL CONSTANT AVG FIT" << "   " << endl;
          fits.push_back( this_fit );
          previous = active_data;
        }
        else{delete this_fit;}
  
        /*cout << "==================================================================================" << endl;
        cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
        cout << data.print_data(active_data) << endl;
        cout << this_fit->report() << endl;
        cout << "==================================================================================" << endl; */
      }
    cout << "-------------------------------------------------" << endl;   
    delete x_t_low_dt; delete x_t_high_dt; 
    } //next t_low
    t_pairs.clear();
  } //end constant fits
  //************************************************************************************************
  

  /* pick the best constant fit to choose a start Form Factor for one_exp fits */

  int t_min_cnst; int t_max_cnst;


  param_value F_cnst(control.F_start, control.F_err); F_cnst.minos = control.minos;

  if( fits.size() > 0 )
  {
    map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
    minuit_fit_result                           best_cnst = (ordered.begin()->second)->get_result();
  
    {
      vector<bool> active = (ordered.begin()->second)->get_active_data();


      /* in the current case, the indexing of active_data corresponds to the t_slice */
        t_min_cnst = 0; t_max_cnst = control.dt[0];


        while( !active[t_min_cnst] ){ t_min_cnst++; if(t_min_cnst > t_max_cnst ){t_min_cnst = 0; break;};}
        while( !active[t_max_cnst] ){ t_max_cnst--; if(t_min_cnst > t_max_cnst ){t_max_cnst = 0; break;};}



    }
  
    F_cnst = (best_cnst.par_values.find("F"))->second;
  }


  else{t_min_cnst = std::get<1>(control.tmin_max[0]); t_max_cnst = std::get<2>(control.tmin_max[0]);}

  //************************************************************************************************
  

  //************************************************************************************************
  /* perform one and two exp fits */
  if(!control.only_cnst)
  {
    // src_exp fits
    {
      map<string, param_value> start_params;

      F_cnst.fixed = false;
      vector<PairIntAbscissa*> x_tlow;
      vector<PairIntAbscissa*> x_thigh;

      {
        F_cnst.error *= 5.0; /* boost the error on the mass */   //Constraints?
        start_params.insert( make_pair("F", F_cnst) );
        
        param_value Ei(2.0,2.0);
        Ei.low_limited = true; Ei.low_limit = control.Ei_min[0]; Ei.fixed = false;  // Constraints
        Ei.high_limited = true; Ei.high_limit = control.Ei_max[0];
        Ei.minos = control.minos;

        stringstream fi; fi << "Fi" << "_Dt" << control.dt[0];
        stringstream ei; ei << "Ei" << "_Dt" << control.dt[0];

        start_params.insert( make_pair(ei.str(), Ei ) );
        start_params.insert( make_pair(fi.str(), F_cnst) );

      }
    
      vector<bool> previous( accepted_data.size(), false );

      std::get<1>(control.tmin_max[0]) = t_min_cnst;
      std::get<2>(control.tmin_max[0]) = t_max_cnst;



      std::vector<pair<pair<int,int>,pair<int,int>>> t_pairs  = get_single_dt_t_ranges(control,0,1,0);

      //looping over all t lows and highs
      cout << "-------------------------------------------------" << endl;  
      for(auto t_pair = t_pairs.begin() ; t_pair != t_pairs.end(); t_pair++){ 

        PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pair->first);
        PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pair->second);

        cout << "src_exp__low = " << *x_t_low_dt << " | ";
        cout << "src_exp__high  = " << *x_t_high_dt << endl;
    
        vector<bool> active_data = !(data.make_active_data()); //all false

        Abscissa* xl = new PairIntAbscissa(make_pair(control.dt[0], t_pair->first.second - 1));
        Abscissa* xh = new PairIntAbscissa(make_pair(control.dt[0], t_pair->second.second + 1));
        
        active_data = active_data || data_in_x_range( data, make_pair(xl,xh) );
        active_data = active_data && accepted_data;
        delete xl; delete xh;

        if( (count_active(active_data) >= control.Nt_min[0]) && (active_data != previous) )
        {
           
          stringstream name; name << "csrc";

          name << "_Dt" << x_t_low_dt->get_x().first << "_" << x_t_low_dt->get_x().second << "-" << x_t_high_dt->get_x().second;

    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_src_exp, start_params, minuit_controls, control.correlated, name.str() );

          if(this_fit->is_valid())
          {
            cout << "   " << "SUCCESSFUL SOURCE EXPONENTIAL AVG FIT" << "   " << endl;
            fits.push_back( this_fit );
            fits_src_exp.push_back( this_fit );
            previous = active_data;
          }
          else{delete this_fit;}

    
          /*cout << "==================================================================================" << endl;
          cout << "*** ONE EXP , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
      cout << "-------------------------------------------------" << endl;
      delete x_t_low_dt; delete x_t_high_dt;   
      } // next t_pairs
      t_pairs.clear();
    } //end of src_exp fits


    //************************************************************************************************

    //************************************************************************************************


    // snk_exp_fits
    {
      map<string, param_value> start_params;
      vector<PairIntAbscissa*> x_tlow;
      vector<PairIntAbscissa*> x_thigh;

      {
        F_cnst.error *= 5.0; /* boost the error on the mass */   //Constraints?
        start_params.insert( make_pair("F", F_cnst) );
        
        param_value Ef(2.0,2.0);
        Ef.low_limited = true; Ef.low_limit = control.Ef_min[0];   // Constraints
        Ef.high_limited = true; Ef.high_limit = control.Ef_max[0];   // Constraints
        Ef.minos = control.minos;
        Ef.fixed = false;

        stringstream ff; ff << "Ff" << "_Dt" << control.dt[0];
        stringstream ef; ef << "Ef" << "_Dt" << control.dt[0];

        start_params.insert( make_pair(ef.str(), Ef ) );
        start_params.insert( make_pair(ff.str(), F_cnst) );
        
      }

    
      vector<bool> previous( accepted_data.size(), false );

      std::vector<pair<pair<int,int>,pair<int,int>>> t_pairs  = get_single_dt_t_ranges(control,0,0,1);

      //looping over all t_mins and t_maxs
      cout << "-------------------------------------------------" << endl;   
      for(auto t_pair = t_pairs.begin() ; t_pair != t_pairs.end() ; t_pair++){

        PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pair->first);
        PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pair->second);

        cout << "snk_exp__low =" << *x_t_low_dt << "  | ";
        cout << "snk_exp__high  =" << *x_t_high_dt << endl;

    
        vector<bool> active_data = !(data.make_active_data()); //all false

        Abscissa* xl = new PairIntAbscissa(make_pair(control.dt[0], t_pair->first.second - 1));
        Abscissa* xh = new PairIntAbscissa(make_pair(control.dt[0], t_pair->second.second + 1));
        
        active_data = active_data || data_in_x_range( data, make_pair(xl,xh) );
        active_data = active_data && accepted_data;
        delete xl; delete xh;

        if( (count_active(active_data) >= control.dt.size() * control.Nt_min[0]) && (active_data != previous) )
        {

          
          stringstream name; name << "csnk";

          name << "_Dt" << x_t_low_dt->get_x().first << "_" << x_t_low_dt->get_x().second << "-" << x_t_high_dt->get_x().second;

    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_snk_exp, start_params, minuit_controls, control.correlated, name.str() );

          if(this_fit->is_valid())
          {
            cout << "   " << "SUCCESSFUL SINK EXPONENTIAL AVG FIT" << "   " << endl;
            fits.push_back( this_fit );
            fits_snk_exp.push_back( this_fit );
            previous = active_data;
          }
          else{delete this_fit;}

    
          /*cout << "==================================================================================" << endl;
          cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
      cout << "-------------------------------------------------" << endl; 
      delete x_t_low_dt; delete x_t_high_dt;  
      } // next t_pair
      t_pairs.clear();
    } //end of snk_exp fits


    //************************************************************************************************
  
    /* pick the best cnst_src_exp fit and cnst_snk_exp to choose a start mass for cnst_two_exp fits */

    param_value F_cnst_one_exp(control.F_start, control.F_err);
    F_cnst_one_exp.minos = control.minos;
    F_cnst_one_exp.fixed = false;
    

    param_value Ei_cnst_src_exp(2.0,2.0); 
    Ei_cnst_src_exp.low_limited = true; Ei_cnst_src_exp.low_limit = control.Ei_min[0];   // Constraints
    Ei_cnst_src_exp.high_limited = true; Ei_cnst_src_exp.high_limit = control.Ei_max[0];   // Constraints    
    Ei_cnst_src_exp.minos = control.minos;
    Ei_cnst_src_exp.fixed = false;

    param_value Fi_cnst_src_exp(control.F_start, control.F_err);
    Fi_cnst_src_exp.minos = control.minos;
    Fi_cnst_src_exp.fixed = false;

    param_value Ef_cnst_snk_exp(2.0,2.0);
    Ef_cnst_snk_exp.low_limited = true; Ef_cnst_snk_exp.low_limit = control.Ef_min[0];   // Constraints
    Ef_cnst_snk_exp.high_limited = true; Ef_cnst_snk_exp.high_limit = control.Ef_max[0];   // Constraints   
    Ef_cnst_snk_exp.minos = control.minos;
    Ef_cnst_snk_exp.fixed = false;

    param_value Ff_cnst_snk_exp(control.F_start, control.F_err);
    Ff_cnst_snk_exp.minos = control.minos; 
    Ff_cnst_snk_exp.fixed = false;

    int t_min_cnst_one_exp = 0; int t_max_cnst_one_exp = control.dt[0];


    if( fits.size() > 0 )
    {
      {

        map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
        minuit_fit_result                           best_cnst_one_exp = (ordered.begin()->second)->get_result();

      
        {
          vector<bool> active = (ordered.begin()->second)->get_active_data();

          /* in the current case, the indexing of active_data corresponds to the timeslice */

          while( !active[t_min_cnst_one_exp] ){ t_min_cnst_one_exp++; if(t_min_cnst_one_exp > t_max_cnst_one_exp ){t_min_cnst_one_exp = 0; break;}; }
          while( !active[t_max_cnst_one_exp] ){ t_max_cnst_one_exp--; if(t_min_cnst_one_exp > t_max_cnst_one_exp ){t_max_cnst_one_exp = 0; break;}; }

        }
      
        F_cnst_one_exp = (best_cnst_one_exp.par_values.find("F"))->second;

      }

      if( fits_src_exp.size() > 0 )
      {

        map<double, AvgFit*, std::greater<double> > ordered_src = make_ordered_list( fits_src_exp, fit_qual );
        minuit_fit_result                           best_cnst_src_exp = (ordered_src.begin()->second)->get_result();
  
      
        stringstream fi; fi << "Fi" << "_Dt" << control.dt[0];
        stringstream ei; ei << "Ei" << "_Dt" << control.dt[0];

        Fi_cnst_src_exp = (best_cnst_src_exp.par_values.find(fi.str()))->second;
        Ei_cnst_src_exp = (best_cnst_src_exp.par_values.find(ei.str()))->second;

      }

      if( fits_snk_exp.size() > 0 )
      {

         map<double, AvgFit*, std::greater<double> > ordered_snk = make_ordered_list( fits_snk_exp, fit_qual );
         minuit_fit_result                           best_cnst_snk_exp = (ordered_snk.begin()->second)->get_result();

         stringstream ff; ff << "Ff" << "_Dt" << control.dt[0];
         stringstream ef; ef << "Ef" << "_Dt" << control.dt[0];

         Ff_cnst_snk_exp = (best_cnst_snk_exp.par_values.find(ff.str()))->second;
         Ef_cnst_snk_exp = (best_cnst_snk_exp.par_values.find(ef.str()))->second;  

      }  
    } 


    else
    {
      int max = control.dt.size();
      for(int i = 0 ; i != max; i++){

      t_min_cnst_one_exp = std::get<1>(control.tmin_max[i]) ;
      t_max_cnst_one_exp = control.dt[i], std::get<2>(control.tmin_max[i]);
      }
    } 

  
    //************************************************************************************************
    /* perform two exp fits */
    if(!control.only_one_exp)
    {

      {

        map<string, param_value> start_params;
        vector<PairIntAbscissa*> x_tlow;
        vector<PairIntAbscissa*> x_thigh;

        {
          F_cnst_one_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair("F", F_cnst_one_exp) );

          stringstream ff; ff << "Ff" << "_Dt" << control.dt[0];
          stringstream ef; ef << "Ef" << "_Dt" << control.dt[0];

          stringstream fi; fi << "Fi" << "_Dt" << control.dt[0];
          stringstream ei; ei << "Ei" << "_Dt" << control.dt[0];

          Fi_cnst_src_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair(fi.str(), Fi_cnst_src_exp) );

          Ff_cnst_snk_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair(ff.str(), Ff_cnst_snk_exp) );

          Ei_cnst_src_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair(ei.str(), Ei_cnst_src_exp) );

          Ef_cnst_snk_exp.error *= 5.0; /* boost the error on the from factor */
          start_params.insert( make_pair(ef.str(), Ef_cnst_snk_exp) );

        }


      vector<bool> previous( accepted_data.size(), false );

      std::get<1>(control.tmin_max[0]) =  t_min_cnst_one_exp;
      std::get<2>(control.tmin_max[0]) =  t_max_cnst_one_exp;


      std::vector<pair<pair<int,int>,pair<int,int>>> t_pairs  = get_single_dt_t_ranges(control,0,1,1);

        // looping over all tmins and tmaxs
        cout << "-------------------------------------------------" << endl; 
        for(auto t_pair = t_pairs.begin() ; t_pair != t_pairs.end() ; t_pair++){

          PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pair->first);
          PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pair->second);

          x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);

      
          vector<bool> active_data = !(data.make_active_data()); //all false


          cout  << "two_exp_t__low ="  << *x_t_low_dt << " |  " ;
          cout << "two_exp_t__high =" << *x_t_high_dt << endl;


          Abscissa* xl = new PairIntAbscissa(make_pair(control.dt[0], t_pair->first.second - 1));
          Abscissa* xh = new PairIntAbscissa(make_pair(control.dt[0], t_pair->second.second + 1));
        
          active_data = active_data || data_in_x_range( data, make_pair(xl,xh) );
          active_data = active_data && accepted_data;
          delete xl; delete xh;


          if( (count_active(active_data) >= control.Nt_min[0]) && (active_data != previous) )
          {
              
            stringstream name; name << "c2";
            name << "_Dt" << x_t_low_dt->get_x().first << "_" << x_t_low_dt->get_x().second << "-" << x_t_high_dt->get_x().second;

      
            AvgFit* this_fit = new AvgFit( data, active_data, cnst_two_exp, start_params, minuit_controls, control.correlated, name.str() );
            if(this_fit->is_valid())
            {
              cout << "   " << "SUCCESSFUL SOURCE AND SINK EXP AVG FIT" << "   " << endl;
              fits.push_back( this_fit );
              previous = active_data;
            }
            else{delete this_fit;}

      
            /*cout << "==================================================================================" << endl;
            cout << "*** 2 EXP , tmin = " << tlow  << " ***" << endl;
            cout << data.print_data(active_data) << endl;
            cout << this_fit->report() << endl;
            cout << "==================================================================================" << endl;*/
          }  
          cout << "-------------------------------------------------" << endl;  
          delete x_t_low_dt; delete x_t_high_dt;     
        } // next dt
        t_pairs.clear();
      }//end two exp fits

    }//end of !only exp loop

    
  }//end of !only cnst loop

  //************************************************************************************************

  //************************************************************************************************
  /* done with avg fitting, select the ranges ... */

  if(fits.size() > 0){
    map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
    minuit_fit_result                           best_cnst = (ordered.begin()->second)->get_result();

    if(ordered.size()){out.success = true;}
    auto end = ordered.begin();
    if(ordered.size() > 20){std::advance(end, 20);}
    else{end = ordered.end();}

  
    {
      for(auto it=ordered.begin(); it!=end; ++it){

        if(it->second->get_chisq().chisq_per_ndof() < chisq_ndof_cutoff)
        {
          vector<bool> active = (it->second)->get_active_data();
          
          /* in the current case, the indexing of active_data corresponds to the timeslice*/
          int t_min_best = 0; int t_max_best = control.dt[0];


          while( !active[t_min_best] ){ t_min_best++; if(t_min_best > t_max_best ){t_min_best = 0; break;};  }
          while( !active[t_max_best] ){ t_max_best--; if(t_min_best > t_max_best ){t_max_best = 0; break;}; }

          std::get<1>(control.tmin_max[0]) = t_min_best;  std::get<2>(control.tmin_max[0]) = t_max_best; 
          pair<pair<int,int>,pair<int,int>> range = make_pair(make_pair(std::get<0>(control.tmin_max[0]),t_min_best),make_pair(std::get<0>(control.tmin_max[0]),t_max_best) );

          out.avg_fits.insert(make_pair(it->second, range));
          active.clear();
        }

      }

    }

    /* clean up pointers after best n */
    for (auto it = end ; it != ordered.end(); ++it){delete it->second;} 
    fits.clear();    fits_src_exp.clear();    fits_snk_exp.clear();  ordered.clear(); 

  }

  else{
    out.success = false;
  }

  
  out.dt = control.dt[0];


  //************************************************************************************************

    
  //************************************************************************************************




  cout << "*******************************" << endl;;
  cout << "   " << "Fixed range for DT = " << control.dt[0] << endl;
  cout << "*******************************" << endl;;


  delete cnst; delete cnst_src_exp; delete cnst_two_exp; delete cnst_snk_exp;

  return out;

};


//************************************************************************
// FIT ALL SOURCE-SINK SEPERATION
//************************************************************************


fit_three_point_output fit_three_point_corr( const Data& data,                        
				    const vector<bool> accepted_data,       
				    fit_three_point_control& control,
				    FitQuality* fit_qual, vector<single_dt_avg_fit_output> single_dt_fits,                 
				    double chisq_ndof_cutoff            
				    )
{
  /* set minuit controls using defaults */
  MinuitControl minuit_controls; minuit_controls.minos = control.minos; 
  minuit_controls.limit_max_fcn_calls = control.limit_max_fcn_calls; 
  minuit_controls.max_fcn_calls = control.max_fcn_calls; 
  
  /* hold the fits */
  vector<AvgFit*> fits;


  /* the output object */
  fit_three_point_output out;

  //************************************************************************************************

  //************************************************************************************************
      /*
      TAKE THE BEST FIT INPUTS IN EACH DT AND TRY TO DO A COMBINED FIT
      */
  //************************************************************************************************

  /* perform the fits */

  std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >> t_pairs;
  std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> > tmp;

  for(auto single_dt_fit = single_dt_fits.begin(); single_dt_fit != single_dt_fits.end(); single_dt_fit++){
    if(single_dt_fit->success){

      map<AvgFit*, pair< pair<int,int>,pair<int,int> > >::iterator it;

      for(it = single_dt_fit->avg_fits.begin(); it != single_dt_fit->avg_fits.end(); it++){
        tmp.push_back(make_pair(it->first, it->second));
      }
      t_pairs.push_back(tmp);
      tmp.clear();
    
    }
  }


  std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >> dt_fits;
  std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>>> t_pairsTemp;
  cart_product(dt_fits, t_pairsTemp, t_pairs.begin(), t_pairs.end());
  t_pairsTemp.clear(); t_pairs.clear();
  
  cout << "-------------------------------------------------" << endl; 
  
  Function* fit_fn;
  
  int max_fits_count = dt_fits.size();
  for(int fits_count = 0; fits_count != max_fits_count; fits_count++)
  {

    map<string, param_value>  par_vals, tmp_val;
    vector<bool> active_data = !(data.make_active_data()); //all false
    vector<bool> no_active( accepted_data.size(), false ); 
    stringstream name; 
    vector<int> src;
    vector<int> snk;
    double F_tmp = 0;

    vector<PairIntAbscissa*> x_tlow;
    vector<PairIntAbscissa*> x_thigh;


    int max_dt_count = dt_fits[fits_count].size();
    for(int dt_count = 0; dt_count != max_dt_count; dt_count++)
    {

      PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(dt_fits[fits_count][dt_count].second.first);
      PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(dt_fits[fits_count][dt_count].second.second);

      x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);

      Abscissa* xl = new PairIntAbscissa(make_pair(control.dt[dt_count], dt_fits[fits_count][dt_count].second.first.second - 1 ));
      Abscissa* xh = new PairIntAbscissa(make_pair(control.dt[dt_count], dt_fits[fits_count][dt_count].second.second.second + 1));
        
      active_data = active_data || data_in_x_range( data, make_pair(xl,xh) );
      delete xl; delete xh;
      
      tmp_val = dt_fits[fits_count][dt_count].first->get_result().par_values;
      par_vals.insert(tmp_val.begin(), tmp_val.end() );
      F_tmp += (tmp_val.find("F"))->second.value/max_dt_count;
      tmp_val.clear();

      name << dt_fits[fits_count][dt_count].first->get_fit_name() ;

      if(dt_count != max_dt_count-1){name << "x";}

      if((dt_fits[fits_count][dt_count].first->get_fit_name().find("csrc_") != std::string::npos) || (dt_fits[fits_count][dt_count].first->get_fit_name().find("c2_") != std::string::npos) ){src.push_back(1);}
      else{src.push_back(0);}

      if((dt_fits[fits_count][dt_count].first->get_fit_name().find("csnk_") != std::string::npos) || (dt_fits[fits_count][dt_count].first->get_fit_name().find("c2_") != std::string::npos) ){snk.push_back(1);}
      else{snk.push_back(0);}

      delete x_t_low_dt; delete x_t_high_dt;

    } //looping over all dt 

    par_vals.find("F")->second.value = F_tmp;

    active_data = active_data && accepted_data;

    int Nt = 0;
    accumulate(control.Nt_min.begin(), control.Nt_min.end(), Nt);

    if( (count_active(active_data) >= Nt) && (active_data != no_active) )
    {
      
      fit_fn = new ThreePointtimesliceCorrNExp( src, snk, control.dt );  


      AvgFit* this_fit = new AvgFit( data, active_data, fit_fn, par_vals, minuit_controls, control.correlated, name.str() );
      if(this_fit->is_valid()){
        cout << "                         " << "SUCCESSFUL AVG FIT" << "   " << endl;
        cout << "----------------------------------------------------------------------------------" << endl; 
        cout << this_fit->report() << endl;
        fits.push_back( this_fit );
      }
      else{delete this_fit;}


      /*cout << "==================================================================================" << endl;
      cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
      cout << data.print_data(active_data) << endl;
      cout << this_fit->report() << endl;
      cout << "==================================================================================" << endl; */
    }

    x_tlow.clear(); x_thigh.clear();
    src.clear();  snk.clear(); active_data.clear(); no_active.clear();


  }



  //************************************************************************************************

    
  //************************************************************************************************
  /* done with avg fitting, select a best fit, do an ensem fit ... */

  cout << endl << "# fits = " << fits.size() << endl;
  if(!fits.size()){
    cerr << "No suitale fits for ensemfitter. Exiting..." << endl;
    out.fit_summary = "No fits found";
    EnsemReal F; F.resize(1);
    out.F = F;
    return out;
  }

  map<double, AvgFit*, std::greater<double> > ordered_fits = make_ordered_list( fits, fit_qual );

  vector<AvgFit*> best_fits;
  auto end_best_fits = ordered_fits.begin();
  std::advance(end_best_fits, 20);

  for(auto fit_n = ordered_fits.begin(); fit_n != end_best_fits; fit_n++){
    best_fits.push_back(fit_n->second);
  }
  
  FitSelector fit_selector( best_fits, fit_qual, chisq_ndof_cutoff );


  if(!fit_selector.ensem_success()){

    auto cent_best_fits = ordered_fits.begin();
    std::advance(cent_best_fits, 100);

   for(auto fit_n = end_best_fits; fit_n != cent_best_fits; fit_n++){
    best_fits.push_back(fit_n->second);
    }   
    FitSelector fit_selector( best_fits, fit_qual, chisq_ndof_cutoff );
  }
  
  if(!fit_selector.ensem_success()){FitSelector fit_selector( fits, fit_qual, chisq_ndof_cutoff );}

  //************************************************************************************************

  cout << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  
  cout << "success=" << fit_selector.ensem_success();

  cout << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

  
  //************************************************************************************************
  /* build the output object */

  out.success     = fit_selector.ensem_success();

  if( !out.success  || control.long_log ){
    out.fit_long_log = fit_selector.get_long_log();
    cout << fit_selector.get_summary();
  }

  if( !out.success ){
    cerr << "no successful ensemfits found" << endl; 
    out.fit_summary = "No fits found";
    EnsemReal F; F.resize(1);
    out.F = F;
    return out;
  }


  out.fit_summary = fit_selector.get_summary();


  if(out.success)
  {

    EnsemFit* efit = fit_selector.get_ensem_fit();
    cout << endl << fit_selector.get_summary() << endl;
    cout << endl << efit->report() << endl;

    param_fit_variation F_param = fit_selector.get_param( "F");
    out.F                       = F_param.ensem;
    out.F_fit_variation         = F_param.report();

    /* make the best fit data description in case you want to look for outliers */
    vector<bool> ensem_fit_active_data = fit_selector.get_ensem_fit()->get_active_data();
    vector<Abscissa*> x;
    vector<double>    y, y_err;

    Function* best_fit_fn = (fit_selector.get_ensem_fit())->get_function();  
    ThreePointtimesliceCorrFunction* ft = static_cast<ThreePointtimesliceCorrFunction*>(best_fit_fn);
      
    map<string, ensem_param_value> ensem_pars = (fit_selector.get_ensem_fit())->get_ensem_param_values();
    
    {
      //int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
      data.get_active_data(ensem_fit_active_data, x, y, y_err);
      vector<pair<double,double>> active_timeslices;
      for(auto x_i = x.begin(); x_i != x.end(); x_i++){
        pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(*x_i) )->get_x();
        active_timeslices.push_back( make_pair(double(t_i.first),double(t_i.second)) );
      }

      vector<double> delt;
      for(auto d_i = control.dt.begin(); d_i != control.dt.end(); d_i++){delt.push_back(double(*d_i)); }

      plot_three_point_timeslice_function_data t_y_yerr = plot_three_point_timeslice_ensem_function(ft, ensem_pars, active_timeslices, delt );

      int max = t_y_yerr.central.size();
      for(int i = 0; i != max; i++)
      {

        out.best_data_description.insert( make_pair( make_pair( int(t_y_yerr.central[i].first.first), int(t_y_yerr.central[i].first.second) ),
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
      int prev_dt;

      plot << "## tmin= " << control.tsrc[control.dt.size()-1] << " tmax= ";

      plot << control.dt[0];
      auto begin = control.dt.begin(); std::advance(begin,1);

      for(auto it = begin; it != control.dt.end(); it++){
        plot << "," << *it;
      } 
    
      plot << " nfits= " << accepted_fits.size() << endl;
      

      /* change this to a nice format write */
      {
        char buff[10]; int nn = sprintf(buff, "%7.5f +/- %7.5f", toDouble(mean(out.F)), toDouble(sqrt(variance(out.F))) );
        plot << "## F= " << buff << endl;
      }

      chisq_dof chisq_desc = (fit_selector.get_ensem_fit())->get_chisq();
      plot << "## chisq/ndof= " << chisq_desc.one_line_report() << endl;
      plot << endl << endl;

      x.clear(); y.clear(); y_err.clear();


      { /* write active data */

        plot << "# active data" << endl;
        //int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
        data.get_active_data(ensem_fit_active_data, x, y, y_err);

        prev_dt = control.dt[0];
        int max = x.size();
        for(int i = 0; i != max; i++)
        {
          pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x(); 
          if(t_i.first != prev_dt){ plot << endl << endl << "# active data" << endl;}

          plot <<  t_i.first << " " << t_i.second  << " "; 
          plot << y[i] << " ";
          plot << y_err[i] << endl;

          prev_dt = t_i.first;
        }
          plot << endl << endl;

      }
      
      x.clear(); y.clear(); y_err.clear();

      { /* write inactive data */
        plot << "# inactive data" << endl;
        //int count = data.get_active_data(!ensem_fit_active_data, x, y, y_err);
        
        data.get_active_data(!ensem_fit_active_data, x, y, y_err);
        prev_dt = control.dt[0];
        int max = x.size();  
        for(int i = 0; i != max; i++)
        {
          pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x();
          if(t_i.first != prev_dt){ plot << endl << endl << "# inactive data" << endl;}

          plot <<  t_i.first << " "<< t_i.second  << " ";  
          plot << y[i] << " ";
          plot << y_err[i] << endl;

          prev_dt = t_i.first;
        }
        plot << endl << endl;
      }

      vector<pair<double,double>> plot_t; 

      int max = control.dt.size();
      for(int j = 0; j != max ; j++){
        int n_points = ((control.tsnk[j] - control.tsrc[j]))*10 + 1; double dt = 0.1;
        for(int i = 0; i <= n_points; i++){ plot_t.push_back( make_pair(double(control.dt[j]),double(control.tsrc[j]) + double(i)*dt) ); }
      }
      
      { /* write ensem fit function */
        plot << "# ensem fit" << endl;

        vector<double> delt;
        for(auto it = control.dt.begin(); it != control.dt.end(); it++){delt.push_back(double(*it)); }

        plot_three_point_timeslice_function_data plot_fn = plot_three_point_timeslice_ensem_function(ft, ensem_pars, plot_t, delt );

        prev_dt = control.dt[0];

        int max = plot_fn.central.size();
        for(int i = 0; i != max;  i++)
        {
          pair<double,double> t = plot_fn.central[i].first;
          if(t.first != prev_dt){ plot << endl << endl << "# ensem fit" << endl;}

          plot << t.first << "  " << t.second << "  "
          << plot_fn.lower[i].second   << "  "
          << plot_fn.central[i].second  << "  "
          << plot_fn.upper[i].second   << endl;

          prev_dt = t.first;
        }

        plot << endl << endl;

      }
      
      { /* write avg fit variation functions */ 
        for(vector<AvgFit*>::iterator fit = accepted_fits.begin(); fit != accepted_fits.end(); fit++)
        {
          plot << "# avg fit: " << (*fit)->get_fit_name() << endl;
    
          Function* fit_fn = (*fit)->get_function();  
          ThreePointtimesliceCorrFunction* ft = static_cast<ThreePointtimesliceCorrFunction*>(fit_fn);
    
          map<string, param_value> pars = ( (*fit)->get_result()).par_values;
    
          vector< pair<pair<double,double>,double> >  plot_fn = plot_three_point_timeslice_function(ft, pars, plot_t);
    
          for(auto it = plot_fn.begin(); it != plot_fn.end(); it++){
            pair<double,double> t = it->first;
            if(t.first != prev_dt){plot << endl << endl << "# avg fit: " << (*fit)->get_fit_name() << endl;}

            plot << t.first << "  " << t.second << "  " << it->second << endl;

            prev_dt = t.first;
          }
    
          plot << endl << endl;
        }//next fit variation
      }
      
      out.plot_data = plot.str();
    }// end if plots

  }//end if success

  cout << "*********************************" << endl;;
  
  /* clean up pointers */

  delete fit_fn;

  for (auto it = fits.begin() ; it != fits.end(); ++it){delete (*it);} 
  fits.clear(); 

  return out; 

};


//************************************************************************
// PLOTTING UTILITIES
//************************************************************************

// PLOTTING FOR THE ENSEM

plot_three_point_timeslice_function_data plot_three_point_timeslice_ensem_function( ThreePointtimesliceCorrFunction* fn,
							    const map<string, ensem_param_value>& ensem_pars,
							    const vector<pair<double,double>>& t_values, const vector<double>& delt )
{
  if( ensem_pars.size() == 0 ){ cerr << "plot_three_point_timeslice_ensem_function:: no parameters provided, exiting " << endl; exit(1); }
  int ncfgs = ((ensem_pars.begin())->second).ensem.size();
  
  plot_three_point_timeslice_function_data out;
  out.dt = delt;

  for(auto delt_i = delt.begin(); delt_i != delt.end(); ++delt_i)
  {
    int count = 0;
    int tmin = 0;
    int tmax = 0;

    for(auto it = t_values.begin() ; it != t_values.end(); ++it)
    {
      if(*delt_i == it->first){ count++; if(it->second < tmin){tmin = it->second;} if(it->second > tmax) {tmax = it->second;}}
      else{continue;}
    }
      out.n_points.push_back(count);
      out.tmin.push_back(tmin);
      out.tmax.push_back(tmax);
      out.dt.push_back(*delt_i);
  } /* assuming t_values is ordered ! */
  
  vector< map<string, double> > scaled_down_pars(ncfgs);
  for(map<string, ensem_param_value>::const_iterator par = ensem_pars.begin(); par != ensem_pars.end(); par++)
    {
    EnsemReal down = rescaleEnsemDown( (par->second).ensem );
    for(int cfg = 0; cfg < ncfgs; cfg++){
      (scaled_down_pars[cfg]).insert( make_pair( par->first, toDouble( peekEnsem(down, cfg) ) ) );
    }
  }
  
  for(auto it = t_values.begin(); it != t_values.end(); it++){

    pair<double,double> t = make_pair(it->first,it->second);
    
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

// PLOTTING FOR THE AVG FITS

vector< pair<pair<double,double>,double> > plot_three_point_timeslice_function( ThreePointtimesliceCorrFunction* fn,
						       const map<string, param_value>& pars,
						       const vector<std::pair<double,double>>& t_values)
{
  vector< pair<pair<double,double>,double> > out;

  map<string, double> central_pars;
  for(map<string, param_value>::const_iterator par = pars.begin(); par != pars.end(); par++){
    central_pars.insert( make_pair( par->first, (par->second).value ) );
  }
  
  for(auto it = t_values.begin(); it != t_values.end(); it++){
    pair<double,double> t = make_pair(it->first,it->second);
    double y = (*fn)( t, central_pars );
    out.push_back( make_pair( t , y ) );
  }//next t
  
  return out;
}


//************************************************************************
// FUNCTION TO LOOP OVER ALL THE T_MINS AND T_MAXS
//************************************************************************

// FUNCTION TO LOOP OVER ALL THE T_MINS AND T_MAXS FOR ALL SOURCE SINK SEPERATIONS

std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> get_all_t_ranges(fit_three_point_control control, bool slide, bool src, bool snk,
                  std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> range, int count_dt)
                  
  {

    /* slide = if you want to slide the range across the source sink seperations 
      src    = if you want to decrease the lower limit or keep it constant usually for src exp
      snk    = if you want to increase the upper limit or keep it constant usually for snk exp
    */

    std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> trange_bins_i;
    std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> trange_bins;

    pair<int,int> tmin, tmax;

    if(slide && src && snk){

      //loop over lower limits
      for(int tlow =  std::get<1>(control.tmin_max[count_dt]); tlow >= control.tsrc[count_dt]; tlow--)
      { 
        // loop over upper limits
        for(int thigh =  std::get<2>(control.tmin_max[count_dt]); thigh <= control.tsnk[count_dt]; thigh++)
        { 
          // loop for sliding
          for(int t_slide_low = control.tsrc[count_dt]; t_slide_low <= control.tsnk[count_dt] - (thigh - tlow ); t_slide_low++){

            trange_bins_i = range;


            int t_slide_high = t_slide_low + thigh - tlow;

            //cout << "....t_slide_low" << t_slide_low << "....t_slide_high" << t_slide_high << endl;
            {

              tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),t_slide_low);
              tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),t_slide_high);

              int t_bins_i_size = trange_bins_i.size();
              if(t_bins_i_size)
              {
                for(int i = 0; i < t_bins_i_size; i++){
                  trange_bins_i[i].push_back(make_pair( tmin, tmax ));
              }

              trange_bins.insert( trange_bins.end(), trange_bins_i.begin(), trange_bins_i.end() );
              trange_bins_i.clear();
            }

            else{
              std::vector<pair<pair<int,int>,pair<int,int>>> tmp; tmp.push_back(make_pair( tmin, tmax ));
              trange_bins.push_back(tmp);
              tmp.clear();
            }

            }
          } //next tslide
        } // next thigh
      } // next tlow


    }

    else if(!slide && !src && snk){

      int tlow = std::get<1>(control.tmin_max[count_dt]);

      // loop over upper limits
      for(int thigh =  std::get<2>(control.tmin_max[count_dt]); thigh <= control.tsnk[count_dt]; thigh++)
      { //cout << "...t_high" << thigh << endl ;

        trange_bins_i = range;

        {

          tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),tlow);
          tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),thigh);

          int t_bins_i_size = trange_bins_i.size();
          if(t_bins_i_size)
          {
            for(int i = 0; i < t_bins_i_size; i++){
              trange_bins_i[i].push_back(make_pair( tmin, tmax ));
            }

            trange_bins.insert( trange_bins.end(), trange_bins_i.begin(), trange_bins_i.end() );
            trange_bins_i.clear();
          }

          else{
            std::vector<pair<pair<int,int>,pair<int,int>>> tmp; tmp.push_back(make_pair( tmin, tmax ));
            trange_bins.push_back(tmp);
            tmp.clear();      
          }
        }
      } // next thigh


    }


    else if(!slide && src && !snk){


      int thigh =  std::get<2>(control.tmin_max[count_dt]);

      // loop over lower limits    
      for(int tlow =  std::get<1>(control.tmin_max[count_dt]); tlow >= control.tsrc[count_dt]; tlow--)
      { 

        trange_bins_i = range;

        {

          tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),tlow);
          tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),thigh);

          int t_bins_i_size = trange_bins_i.size();
          if(t_bins_i_size)
          {
            for(int i = 0; i < t_bins_i_size; i++){
              trange_bins_i[i].push_back(make_pair( tmin, tmax ));
            }

            trange_bins.insert( trange_bins.end(), trange_bins_i.begin(), trange_bins_i.end() );
            trange_bins_i.clear();
          }

          else
          {
            std::vector<pair<pair<int,int>,pair<int,int>>> tmp; tmp.push_back(make_pair( tmin, tmax ));
            trange_bins.push_back(tmp);
            tmp.clear();
          }

        }

      } // next tlow



    }

    else if(!slide && src && snk){

      // loop over lower limits          
      for(int tlow =  std::get<1>(control.tmin_max[count_dt]); tlow >= control.tsrc[count_dt]; tlow--)
      { 

        // loop over upper limits
        for(int thigh =  std::get<2>(control.tmin_max[count_dt]); thigh <= control.tsnk[count_dt]; thigh++)
        { 

          trange_bins_i = range;

          {

            tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),tlow);
            tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),thigh);

            int t_bins_i_size = trange_bins_i.size();
            if(t_bins_i_size)
            {
              for(int i = 0; i < t_bins_i_size; i++){
                trange_bins_i[i].push_back(make_pair( tmin, tmax ));
              }

              trange_bins.insert( trange_bins.end(), trange_bins_i.begin(), trange_bins_i.end() );
              trange_bins_i.clear();
            }

            else{
              std::vector<pair<pair<int,int>,pair<int,int>>> tmp; tmp.push_back(make_pair( tmin, tmax ));
              trange_bins.push_back(tmp);
              tmp.clear();
            }

          }

        } // next thigh
      } // next tlow


    }

    
    return trange_bins;

  };


// FUNCTION TO LOOP OVER ALL THE T_MINS AND T_MAXS FOR SINGLE SOURCE SINK SEPERATION

std::vector<pair<pair<int,int>,pair<int,int>>> get_single_dt_t_ranges(fit_three_point_control control, bool slide, bool src, bool snk)
  {

    std::vector<pair<pair<int,int>,pair<int,int>>> trange_bins;

    pair<int,int> tmin, tmax;


    if(slide && src && snk){

      // loop over lower limits
      for(int tlow =  std::get<1>(control.tmin_max[0]); tlow >= control.tsrc[0]; tlow--)
      { 

        // loop over upper limits
        for(int thigh =  std::get<2>(control.tmin_max[0]); thigh <= control.tsnk[0]; thigh++)
        { 

          // loop over sliding
          for(int t_slide_low = control.tsrc[0]; t_slide_low <= control.tsnk[0] - (thigh - tlow ); t_slide_low++){


            int t_slide_high = t_slide_low + thigh - tlow;


            {

              tmin = make_pair(std::get<0>(control.tmin_max[0]),t_slide_low);
              tmax = make_pair(std::get<0>(control.tmin_max[0]),t_slide_high);


              trange_bins.push_back(make_pair( tmin, tmax ));

            }
          } //next tslide
        } // next thigh
      } // next tlow


    }

    else if(!slide && !src && snk){

      int tlow = std::get<1>(control.tmin_max[0]);

      // loop over upper limits
      for(int thigh =  std::get<2>(control.tmin_max[0]); thigh <= control.tsnk[0]; thigh++)
      {


        {

          tmin = make_pair(std::get<0>(control.tmin_max[0]),tlow);
          tmax = make_pair(std::get<0>(control.tmin_max[0]),thigh);

          trange_bins.push_back(make_pair( tmin, tmax ));
      

        }
      } // next thigh


    }


    else if(!slide && src && !snk){


      int thigh =  std::get<2>(control.tmin_max[0]);

      // loop over lower limits    
      for(int tlow =  std::get<1>(control.tmin_max[0]); tlow >= control.tsrc[0]; tlow--)
      { 

        {

          tmin = make_pair(std::get<0>(control.tmin_max[0]),tlow);
          tmax = make_pair(std::get<0>(control.tmin_max[0]),thigh);

          trange_bins.push_back(make_pair( tmin, tmax ));

        }

      } // next tlow



    }

    else if(!slide && src && snk){


      // loop over lower limits         
      for(int tlow =  std::get<1>(control.tmin_max[0]); tlow >= control.tsrc[0]; tlow--)
      { 

        // loop over upper limits
        for(int thigh =  std::get<2>(control.tmin_max[0]); thigh <= control.tsnk[0]; thigh++)
        { 

          {

            tmin = make_pair(std::get<0>(control.tmin_max[0]),tlow);
            tmax = make_pair(std::get<0>(control.tmin_max[0]),thigh);

            trange_bins.push_back(make_pair( tmin, tmax ));

          }

        } // next thigh
      } // next tlow


    }

    
    return trange_bins;

  };


// CARTESIAN PRODUCT 

void cart_product(
    std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >>& rvvi,  // final result
    std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>>>&  rvi,   // current result 
    std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >>::const_iterator me, // current input
    std::vector< std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>> >>::const_iterator end) // final input
  {
      if(me == end) {
          // terminal condition of the recursion. We no longer have
          // any input vectors to manipulate. Add the current result (rvi)
          // to the total set of results (rvvvi).
          rvvi.push_back(rvi);
          return;
      }

      // need an easy name for my vector-of-ints
      const std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>>>& mevi = *me;
      for(std::vector<pair<AvgFit*, pair<pair<int,int>,pair<int,int>>>>::const_iterator it = mevi.begin();
          it != mevi.end();
          it++) {
          // final rvi will look like "a, b, c, ME, d, e, f"
          // At the moment, rvi already has "a, b, c"
          rvi.push_back(*it);  // add ME
          cart_product(rvvi, rvi, me+1, end); //add "d, e, f"
          rvi.pop_back(); // clean ME off for next round
      }
  };


