#include "three_point_timeslice_fitting.h"

//************************************************************************
// DRIVER OVER DIFFERENT SOURCE-SINK SEPERATIONS
//************************************************************************

fit_three_point_output fit( const vector<Data>& data,                        
				    const vector<vector<bool>> accepted_data,       
				    vector<fit_three_point_control> control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff,
            int num               
				    )

{

  fit_three_point_output out;
  std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > range;

  for(int i = 0; i < num-1; i++){
    

    range = get_range(data[i], accepted_data[i], control[i], fit_qual, chisq_ndof_cutoff, range, i);
    for(int j = 0; j <= i; j++){control[i+1].tmin_max[j] = control[i].tmin_max[j];}

  }

  out = fit_three_point_corr(data[num-1], accepted_data[num-1], control[num-1], fit_qual, chisq_ndof_cutoff, range, num-1);

  return out;

};

//************************************************************************
// FIX THE RANGES FOR EACH DT
//************************************************************************

std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > get_range( const Data& data,                        
				    const vector<bool> accepted_data,       
				    fit_three_point_control& control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff, std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > range, 
            int count_dt               
				    )
{
  /* set minuit controls using defaults */
  MinuitControl minuit_controls; minuit_controls.minos = control.minos;

  /* the output object */
  std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > out;
  std::vector<pair<pair<int,int>,pair<int,int>>> out_i;

  
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

  Function* cnst = new ThreePointtimesliceCorrNExp( 0, 0, control.dt );  
  Function* cnst_src_exp = new ThreePointtimesliceCorrNExp( 1, 0, control.dt );  
  Function* cnst_snk_exp = new ThreePointtimesliceCorrNExp( 0, 1, control.dt ); 
  Function* cnst_two_exp = new ThreePointtimesliceCorrNExp( 1, 1, control.dt );


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


    std::get<1>(control.tmin_max[count_dt]) = int(std::get<0>(control.tmin_max[count_dt])/2 - control.Nt_min/2);
    std::get<2>(control.tmin_max[count_dt]) = int(std::get<0>(control.tmin_max[count_dt])/2 + control.Nt_min/2);


  
    vector<bool> previous( accepted_data.size(), false ); 

    std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,1,1,1, range, count_dt);

    cout << t_pairs.size() << endl;


    for(int i = 0 ; i < t_pairs.size() ; i++){
    cout << "----------------------------------------------------------------------------------" << endl;

      for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
      {

        PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
        PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

        x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);


      } //looping over all dt 

      vector<bool> active_data = !(data.make_active_data()); //all false
        
      for(int j = 0 ; j < control.tmin_max.size() ; j++){

        cout  << "cnst__t_low  =  "<< *x_tlow[j] << " | ";
        cout << "cnst__t_high  =  " << *x_thigh[j] << endl;
        active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
      }

      active_data = active_data && accepted_data;


      if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
      {
        cout << "***  SUCCESSFUL CONST AVG FIT  ***" << endl;
        cout << "----------------------------------------------------------------------------------" << endl;
        stringstream name; name << "c";
        for(int k = 0 ; k < control.tmin_max.size() ; k++){
          name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
        }
  
        AvgFit* this_fit = new AvgFit( data, active_data, cnst, start_params, minuit_controls, control.correlated, name.str() );
        fits.push_back( this_fit );
        //fits_src_exp.push_back( this_fit );
        //fits_snk_exp.push_back( this_fit );
        previous = active_data;
  
        /*cout << "==================================================================================" << endl;
        cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
        cout << data.print_data(active_data) << endl;
        cout << this_fit->report() << endl;
        cout << "==================================================================================" << endl; */
      }

      x_tlow.clear(); x_thigh.clear();
    } //next t_low
    t_pairs.clear();
  } //end constant fits
  //************************************************************************************************
  

  /* pick the best constant fit to choose a start Form Factor for one_exp fits */

  vector<std::pair<int,int>> tmin_cnst;
  vector<std::pair<int,int>> tmax_cnst;

  param_value F_cnst(control.F_start, control.F_err); F_cnst.minos = control.minos;

  if( fits.size() > 0 )
  {
    map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
    minuit_fit_result                           best_cnst = (ordered.begin()->second)->get_result();
  
    {
      vector<bool> active = (ordered.begin()->second)->get_active_data();
      int dt_count = 0;
      for(int i = 0 ; i < control.dt.size() ; i++){
      /* in the current case, the indexing of active_data corresponds to the count */
        int t_min_cnst = 0; int t_max_cnst = control.dt[i];


        while( !active[dt_count + t_min_cnst] ){ t_min_cnst++; if(t_min_cnst > t_max_cnst ){t_min_cnst = 1; break;};}
        while( !active[dt_count + t_max_cnst] ){ t_max_cnst--; if(t_min_cnst > t_max_cnst ){t_max_cnst = -1; break;};}
        t_min_cnst--;t_max_cnst++;
        tmin_cnst.push_back(make_pair(control.dt[i], t_min_cnst)); tmax_cnst.push_back(make_pair(control.dt[i], t_max_cnst));

        dt_count += control.dt[i]+1;
      }

    }
  
    F_cnst = (best_cnst.par_values.find("F"))->second;
  }


  else{
    for(int i = 0 ; i < control.dt.size() ; i++){

      tmin_cnst.push_back(make_pair(control.dt[i], std::get<1>(control.tmin_max[i]) ));
      tmax_cnst.push_back(make_pair(control.dt[i], std::get<2>(control.tmin_max[i]) ));
    }
  }

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
        Ei.low_limited = true; Ei.low_limit = 0.0; Ei.fixed = false;  // Constraints
        Ei.minos = control.minos;

        for(int t = 0; t < control.dt.size(); t++){

          stringstream fi; fi << "Fi" << "_Dt" << control.dt[t];
          stringstream ei; ei << "Ei" << "_Dt" << control.dt[t];

          start_params.insert( make_pair(ei.str(), Ei ) );
          start_params.insert( make_pair(fi.str(), F_cnst) );
        }
      }

    
      vector<bool> previous( accepted_data.size(), false );


      for(int i = 0; i < control.tmin_max.size(); i++){ //if the ordering is right it works
        std::get<1>(control.tmin_max[i]) = tmin_cnst[i].second;
        std::get<2>(control.tmin_max[i]) = tmax_cnst[i].second;
      }


      std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,0,1,0, range, count_dt);


      for(int i = 0 ; i < t_pairs.size() ; i++){
      cout << "----------------------------------------------------------------------------------" << endl;
        for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
        {

          PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
          PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

          x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);

          cout << "src_exp__low = " << *x_t_low_dt << " | ";
          cout << "src_exp__high  = " << *x_t_high_dt << endl;

        } //looping over all dt 
    
        vector<bool> active_data = !(data.make_active_data()); //all false
          
        for(int j = 0 ; j < control.tmin_max.size() ; j++){
          active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
        }

        active_data = active_data && accepted_data;

        if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
        {
          cout << "***  SUCCESSFUL SOURCE EXPONENTIAL AVG FIT ***" << endl;
          cout << "----------------------------------------------------------------------------------" << endl;
          stringstream name; name << "csrc";
          for(int k = 0 ; k < control.tmin_max.size() ; k++){
            name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
          }
    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_src_exp, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          fits_src_exp.push_back( this_fit );
          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** ONE EXP , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
        x_tlow.clear(); x_thigh.clear();
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
        Ef.low_limited = true; Ef.low_limit = 0.0;   // Constraints
        Ef.minos = control.minos;
        Ef.fixed = false;

        for(int t = 0; t < control.dt.size(); t++){
          stringstream ff; ff << "Ff" << "_Dt" << control.dt[t];
          stringstream ef; ef << "Ef" << "_Dt" << control.dt[t];

          start_params.insert( make_pair(ef.str(), Ef ) );
          start_params.insert( make_pair(ff.str(), F_cnst) );
        }
        
      }

    
      vector<bool> previous( accepted_data.size(), false );


      // for(int i = 0; i < control.tmin_max.size(); i++){
      //   std::get<1>(control.tmin_max[i]) = tmin_cnst[i].second;
      //   std::get<2>(control.tmin_max[i]) = tmax_cnst[i].second;
      // }

      std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,0,0,1, range, count_dt);


      for(int i = 0 ; i < t_pairs.size() ; i++){
      cout << "----------------------------------------------------------------------------------" << endl;

        for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
        {

          PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
          PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

          x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);

          cout << "snk_exp__low =" << *x_t_low_dt << "  | ";
          cout << "snk_exp__high  =" << *x_t_high_dt << endl;

        } //looping over all dt 
    
        vector<bool> active_data = !(data.make_active_data()); //all false
          
        for(int j = 0 ; j < control.tmin_max.size() ; j++){
          active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
        }

        active_data = active_data && accepted_data;

        if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
        {
          cout << "***  SUCCESSFUL SINK EXPONENTIAL AVG FIT ***" << endl;
          cout << "----------------------------------------------------------------------------------" << endl;
          stringstream name; name << "csnk";
          for(int k = 0 ; k < control.tmin_max.size() ; k++){
            name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
          }
    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_snk_exp, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          fits_snk_exp.push_back( this_fit );
          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
        x_tlow.clear(); x_thigh.clear();
      } // next t_pair
      t_pairs.clear();
    } //end of snk_exp fits


    //************************************************************************************************
  
    /* pick the best cnst_src_exp fit and cnst_snk_exp to choose a start mass for cnst_two_exp fits */

    vector<std::pair<int,int>>  tmin_cnst_one_exp;
    vector<std::pair<int,int>>  tmax_cnst_one_exp; 

    param_value F_cnst_one_exp(control.F_start, control.F_err);
    F_cnst_one_exp.minos = control.minos;
    F_cnst_one_exp.fixed = false;
    
    vector<param_value>  Ei_cnst_src_exp,  Ef_cnst_snk_exp, Fi_cnst_src_exp, Ff_cnst_snk_exp;

    param_value Emi(2.0,2.0); 
    Emi.low_limited = true; Emi.low_limit = 0.0;   // Constraints
    Emi.minos = control.minos;
    Emi.fixed = false;

    param_value Fi(control.F_start, control.F_err);
    Fi.minos = control.minos;
    Fi.fixed = false;

    param_value Emf(2.0,2.0);
    Emf.low_limited = true; Emf.low_limit = 0.0;   // Constraints
    Emf.minos = control.minos;
    Emf.fixed = false;

    param_value Ff(control.F_start, control.F_err);
    Ff.minos = control.minos; 
    Ff.fixed = false;

    for(int t = 0; t < control.dt.size(); t++){
      Ei_cnst_src_exp.push_back(Emi); Ef_cnst_snk_exp.push_back(Emf); 
      Fi_cnst_src_exp.push_back(Fi); Ff_cnst_snk_exp.push_back(Ff);
    }


    if( fits.size() > 0 )
    {
      {

        map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
        minuit_fit_result                           best_cnst_one_exp = (ordered.begin()->second)->get_result();

      
        {
          vector<bool> active = (ordered.begin()->second)->get_active_data();
          int dt_count = 0;
          for(int i = 0 ; i < control.dt.size() ; i++){
          /* in the current case, the indexing of active_data corresponds to the count */
            int t_min_cnst_one_exp = 0; int t_max_cnst_one_exp = control.dt[i];

            while( !active[dt_count + t_min_cnst_one_exp] ){ t_min_cnst_one_exp++; if(t_min_cnst_one_exp > t_max_cnst_one_exp ){t_min_cnst_one_exp = 1; break;}; }
            while( !active[dt_count + t_max_cnst_one_exp] ){ t_max_cnst_one_exp--; if(t_min_cnst_one_exp > t_max_cnst_one_exp ){t_max_cnst_one_exp = -1; break;}; }
            t_min_cnst_one_exp--;t_max_cnst_one_exp++;
            tmin_cnst_one_exp.push_back(make_pair(control.dt[i], t_min_cnst_one_exp)); tmax_cnst_one_exp.push_back(make_pair(control.dt[i], t_max_cnst_one_exp));

            dt_count += control.dt[i]+1;
          }

        }
      
        F_cnst_one_exp = (best_cnst_one_exp.par_values.find("F"))->second;

      }

      if( fits_src_exp.size() > 0 )
      {

        map<double, AvgFit*, std::greater<double> > ordered_src = make_ordered_list( fits_src_exp, fit_qual );
        minuit_fit_result                           best_cnst_src_exp = (ordered_src.begin()->second)->get_result();
  
      
        for(int t = 0; t < control.dt.size(); t++){
          stringstream fi; fi << "Fi" << "_Dt" << control.dt[t];
          stringstream ei; ei << "Ei" << "_Dt" << control.dt[t];

          Fi_cnst_src_exp[t] = (best_cnst_src_exp.par_values.find(fi.str()))->second;
          Ei_cnst_src_exp[t] = (best_cnst_src_exp.par_values.find(ei.str()))->second;
        }


      }

      if( fits_snk_exp.size() > 0 )
      {

         map<double, AvgFit*, std::greater<double> > ordered_snk = make_ordered_list( fits_snk_exp, fit_qual );
         minuit_fit_result                           best_cnst_snk_exp = (ordered_snk.begin()->second)->get_result();

        for(int t = 0; t < control.dt.size(); t++){
          stringstream ff; ff << "Ff" << "_Dt" << control.dt[t];
          stringstream ef; ef << "Ef" << "_Dt" << control.dt[t];

          Ff_cnst_snk_exp[t] = (best_cnst_snk_exp.par_values.find(ff.str()))->second;
          Ef_cnst_snk_exp[t] = (best_cnst_snk_exp.par_values.find(ef.str()))->second;  
        }

      }  
    } 


    else
    {
      for(int i = 0 ; i < control.dt.size() ; i++){

      tmin_cnst.push_back(make_pair(control.dt[i], std::get<1>(control.tmin_max[i]) ));
      tmax_cnst.push_back(make_pair(control.dt[i], std::get<2>(control.tmin_max[i]) ));
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

          for(int t = 0; t < control.dt.size(); t++){
            stringstream ff; ff << "Ff" << "_Dt" << control.dt[t];
            stringstream ef; ef << "Ef" << "_Dt" << control.dt[t];

            stringstream fi; fi << "Fi" << "_Dt" << control.dt[t];
            stringstream ei; ei << "Ei" << "_Dt" << control.dt[t];

            Fi_cnst_src_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(fi.str(), Fi_cnst_src_exp[t]) );

            Ff_cnst_snk_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(ff.str(), Ff_cnst_snk_exp[t]) );

            Ei_cnst_src_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(ei.str(), Ei_cnst_src_exp[t]) );

            Ef_cnst_snk_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(ef.str(), Ef_cnst_snk_exp[t]) );
          }

        }


     vector<bool> previous( accepted_data.size(), false );


      for(int i = 0; i < control.tmin_max.size(); i++){
        std::get<1>(control.tmin_max[i]) =  tmin_cnst_one_exp[i].second;
        std::get<2>(control.tmin_max[i]) =  tmax_cnst_one_exp[i].second;
      }


      std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,0,1,1, range, count_dt);


        for(int i = 0 ; i < t_pairs.size() ; i++){
        cout << "----------------------------------------------------------------------------------" << endl;
          for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
          {

            PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
            PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

            x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);


          } //looping over all dt 
      
          vector<bool> active_data = !(data.make_active_data()); //all false

          for(int j = 0 ; j < control.tmin_max.size() ; j++){

            cout  << "two_exp_t__low ="  << *x_tlow[j] << " |  " ;
            cout << "two_exp_t__high =" << *x_thigh[j] << endl;
            active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
          }

          active_data = active_data && accepted_data;

          if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
          {
          cout << "***  SUCCESSFUL SOURCE AND SINK EXP AVG FIT  ***" << endl;
          cout << "----------------------------------------------------------------------------------" << endl;
          stringstream name; name << "c2";
          for(int k = 0 ; k < control.tmin_max.size() ; k++){
            name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
          }
    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_two_exp, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** 2 EXP , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl;*/
          }  
          x_tlow.clear(); x_thigh.clear();      
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

    auto end = ordered.begin();
    std::advance(end, 1);
  
    {
      for(auto it=ordered.begin(); it!=end; ++it){

        vector<bool> active = (it->second)->get_active_data();
        int dt_count = 0;
        
        for(int i = 0 ; i < control.dt.size() ; i++){
        /* in the current case, the indexing of active_data corresponds to the count */
          int t_min_best = 0; int t_max_best = control.dt[i];


          while( !active[dt_count + t_min_best] ){ t_min_best++; if(t_min_best > t_max_best ){t_min_best = 1; break;};  }
          while( !active[dt_count + t_max_best] ){ t_max_best--; if(t_min_best > t_max_best ){t_max_best = -1; break;}; }
          t_min_best--;t_max_best++;

          std::get<1>(control.tmin_max[i]) = t_min_best;  std::get<2>(control.tmin_max[i]) = t_max_best; 

          dt_count += control.dt[i]+1;
          out_i.push_back(make_pair(make_pair(std::get<0>(control.tmin_max[i]),t_min_best),make_pair(std::get<0>(control.tmin_max[i]),t_max_best) ));
        }

        out.push_back(out_i);

        out_i.clear(); active.clear();

      }

    }

  }

  else{
    std::get<1>(control.tmin_max[count_dt]) = 0;  std::get<2>(control.tmin_max[count_dt]) = 0; 
    out = range;
    for(int i = 0; i < range.size();i++){out[i].push_back(make_pair(make_pair(std::get<0>(control.tmin_max[count_dt]),0),make_pair(std::get<0>(control.tmin_max[count_dt]),0) ));}
  }


  //************************************************************************************************

    
  //************************************************************************************************




  cout << "#################################################################################" << endl;
  cout << "Fixed range for DT = " << control.dt[count_dt] << endl;
  cout << "#################################################################################" << endl;
  /* clean up pointers */
  delete cnst; delete cnst_src_exp; delete cnst_two_exp; delete cnst_snk_exp;


  return out;

};

//************************************************************************
// FIT EACH SOURCE-SINK SEPERATION
//************************************************************************

fit_three_point_output fit_three_point_corr( const Data& data,                        
				    const vector<bool> accepted_data,       
				    fit_three_point_control& control,
				    FitQuality* fit_qual,                   
				    double chisq_ndof_cutoff, std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > range, 
            int count_dt               
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

  Function* cnst = new ThreePointtimesliceCorrNExp( 0, 0, control.dt );  
  Function* cnst_src_exp = new ThreePointtimesliceCorrNExp( 1, 0, control.dt );  
  Function* cnst_snk_exp = new ThreePointtimesliceCorrNExp( 0, 1, control.dt ); 
  Function* cnst_two_exp = new ThreePointtimesliceCorrNExp( 1, 1, control.dt );


  /* perform constant fits */
  {

    map<string, param_value> start_params;
    vector<PairIntAbscissa*> x_tlow;
    vector<PairIntAbscissa*> x_thigh;

    control.F_start = data.get_all_y_mean()[ int(control.dt[0]/2) ];
    control.F_err   = data.get_all_y_err()[ int(control.dt[0]/2) ]*2;
  
    param_value F(control.F_start, control.F_err); //the constant that would be the ff. Constraints?
    F.minos = control.minos;
    F.fixed = false;
    start_params.insert( make_pair("F", F ) ); 


    std::get<1>(control.tmin_max[count_dt]) = int(std::get<0>(control.tmin_max[count_dt])/2 - control.Nt_min/2);
    std::get<2>(control.tmin_max[count_dt]) = int(std::get<0>(control.tmin_max[count_dt])/2 + control.Nt_min/2);


  
    vector<bool> previous( accepted_data.size(), false ); 

    std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,1,1,1, range, count_dt);



    for(int i = 0 ; i < t_pairs.size() ; i++){
    cout << "----------------------------------------------------------------------------------" << endl;
      for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
      {

        PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
        PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

        x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);


      } //looping over all dt 

      vector<bool> active_data = !(data.make_active_data()); //all false
        
      for(int j = 0 ; j < control.tmin_max.size() ; j++){

        cout <<  "cnst_t_low = " << *x_tlow[j] << " | " ;
        cout << "cnst_t_high = " <<  *x_thigh[j] << endl;
        active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
      }

      active_data = active_data && accepted_data;


      if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
      {
        cout << "*** SUCCESSFUL CONST AVG FIT ***" << endl;
        cout << "----------------------------------------------------------------------------------" << endl;
        stringstream name; name << "c";
        for(int k = 0 ; k < control.tmin_max.size() ; k++){
          name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
        }
  
        AvgFit* this_fit = new AvgFit( data, active_data, cnst, start_params, minuit_controls, control.correlated, name.str() );
        fits.push_back( this_fit );
        //fits_src_exp.push_back( this_fit );
        //fits_snk_exp.push_back( this_fit );
        previous = active_data;
  
        /*cout << "==================================================================================" << endl;
        cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
        cout << data.print_data(active_data) << endl;
        cout << this_fit->report() << endl;
        cout << "==================================================================================" << endl; */
      }

      x_tlow.clear(); x_thigh.clear();
    } //next t_low
    t_pairs.clear();
  } //end constant fits
  //************************************************************************************************
  

  /* pick the best constant fit to choose a start Form Factor for one_exp fits */

  vector<std::pair<int,int>> tmin_cnst;
  vector<std::pair<int,int>> tmax_cnst;

  param_value F_cnst(control.F_start, control.F_err); F_cnst.minos = control.minos;

  if( fits.size() > 0 )
  {
    map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
    minuit_fit_result                           best_cnst = (ordered.begin()->second)->get_result();
  
    {
      vector<bool> active = (ordered.begin()->second)->get_active_data();
      int dt_count = 0;
      for(int i = 0 ; i < control.dt.size() ; i++){
      /* in the current case, the indexing of active_data corresponds to the count */
        int t_min_cnst = 0; int t_max_cnst = control.dt[i];


        while( !active[dt_count + t_min_cnst] ){ t_min_cnst++; if(t_min_cnst > t_max_cnst ){t_min_cnst = 1; break;}; }
        while( !active[dt_count + t_max_cnst] ){ t_max_cnst--; if(t_min_cnst > t_max_cnst ){t_max_cnst = -1; break;}; }
        t_min_cnst--;t_max_cnst++;
        tmin_cnst.push_back(make_pair(control.dt[i], t_min_cnst)); tmax_cnst.push_back(make_pair(control.dt[i], t_max_cnst));

        dt_count += control.dt[i]+1;
      }

    }
  
    F_cnst = (best_cnst.par_values.find("F"))->second;
  }


  else{
    for(int i = 0 ; i < control.dt.size() ; i++){

      tmin_cnst.push_back(make_pair(control.dt[i], std::get<1>(control.tmin_max[i]) ));
      tmax_cnst.push_back(make_pair(control.dt[i], std::get<2>(control.tmin_max[i]) ));
    }
  }

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
        Ei.low_limited = true; Ei.low_limit = 0.0; Ei.fixed = false;  // Constraints
        Ei.minos = control.minos;

        for(int t = 0; t < control.dt.size(); t++){

          stringstream fi; fi << "Fi" << "_Dt" << control.dt[t];
          stringstream ei; ei << "Ei" << "_Dt" << control.dt[t];

          start_params.insert( make_pair(ei.str(), Ei ) );
          start_params.insert( make_pair(fi.str(), F_cnst) );
        }
      }


    
      vector<bool> previous( accepted_data.size(), false );

      for(int i = 0; i < control.tmin_max.size(); i++){ //if the ordering is right it works
        std::get<1>(control.tmin_max[i]) = tmin_cnst[i].second;
        std::get<2>(control.tmin_max[i]) = tmax_cnst[i].second;
      }

      std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,0,1,0, range, count_dt);


      for(int i = 0 ; i < t_pairs.size() ; i++){
      cout << "----------------------------------------------------------------------------------" << endl;
        for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
        {

          PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
          PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

          x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);

          cout << "src_exp__low =" << *x_t_low_dt << "  | ";
          cout << "src_exp__high  =" << *x_t_high_dt << endl;

        } //looping over all dt 
    
        vector<bool> active_data = !(data.make_active_data()); //all false
          
        for(int j = 0 ; j < control.tmin_max.size() ; j++){
          active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
        }

        active_data = active_data && accepted_data;

        if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
        {
          cout << "*** SUCCESSFUL SOURCE EXPONENTIAL AVG FIT ***" << endl;
          cout << "----------------------------------------------------------------------------------" << endl;
          stringstream name; name << "csrc";
          for(int k = 0 ; k < control.tmin_max.size() ; k++){
            name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
          }
    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_src_exp, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          fits_src_exp.push_back( this_fit );
          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** ONE EXP , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
        x_tlow.clear(); x_thigh.clear();
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
        Ef.low_limited = true; Ef.low_limit = 0.0;   // Constraints
        Ef.minos = control.minos;
        Ef.fixed = false;

        for(int t = 0; t < control.dt.size(); t++){
          stringstream ff; ff << "Ff" << "_Dt" << control.dt[t];
          stringstream ef; ef << "Ef" << "_Dt" << control.dt[t];

          start_params.insert( make_pair(ef.str(), Ef ) );
          start_params.insert( make_pair(ff.str(), F_cnst) );
        }
        
      }

    
      vector<bool> previous( accepted_data.size(), false );


      // for(int i = 0; i < control.tmin_max.size(); i++){
      //   std::get<1>(control.tmin_max[i]) = tmin_cnst[i].second;
      //   std::get<2>(control.tmin_max[i]) = tmax_cnst[i].second;
      // }

      std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,0,0,1, range, count_dt);



      for(int i = 0 ; i < t_pairs.size() ; i++){
      cout << "----------------------------------------------------------------------------------" << endl;

        for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
        {

          PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
          PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

          x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);

          cout << "snk_exp__low   =" << *x_t_low_dt << "  | ";
          cout << "snk_exp__high  =" << *x_t_high_dt << endl;

        } //looping over all dt 
    
        vector<bool> active_data = !(data.make_active_data()); //all false
          
        for(int j = 0 ; j < control.tmin_max.size() ; j++){
          active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
        }

        active_data = active_data && accepted_data;

        if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
        {
          cout << "*** SUCCESSFUL SINK EXPONENTIAL AVG FIT ***" << endl;
          cout << "----------------------------------------------------------------------------------" << endl;
          stringstream name; name << "csnk";
          for(int k = 0 ; k < control.tmin_max.size() ; k++){
            name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
          }
    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_snk_exp, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          fits_snk_exp.push_back( this_fit );
          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** CONST , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl; */
        }
        x_tlow.clear(); x_thigh.clear();
      } // next t_pair
      t_pairs.clear();
    } //end of snk_exp fits


    //************************************************************************************************
  
    /* pick the best cnst_src_exp fit and cnst_snk_exp to choose a start mass for cnst_two_exp fits */

    vector<std::pair<int,int>>  tmin_cnst_one_exp;
    vector<std::pair<int,int>>  tmax_cnst_one_exp; 

    param_value F_cnst_one_exp(control.F_start, control.F_err);
    F_cnst_one_exp.minos = control.minos;
    F_cnst_one_exp.fixed = false;
    
    vector<param_value>  Ei_cnst_src_exp,  Ef_cnst_snk_exp, Fi_cnst_src_exp, Ff_cnst_snk_exp;

    param_value Emi(2.0,2.0); 
    Emi.low_limited = true; Emi.low_limit = 0.0;   // Constraints
    Emi.minos = control.minos;
    Emi.fixed = false;

    param_value Fi(control.F_start, control.F_err);
    Fi.minos = control.minos;
    Fi.fixed = false;

    param_value Emf(2.0,2.0);
    Emf.low_limited = true; Emf.low_limit = 0.0;   // Constraints
    Emf.minos = control.minos;
    Emf.fixed = false;

    param_value Ff(control.F_start, control.F_err);
    Ff.minos = control.minos; 
    Ff.fixed = false;

    for(int t = 0; t < control.dt.size(); t++){
      Ei_cnst_src_exp.push_back(Emi); Ef_cnst_snk_exp.push_back(Emf); 
      Fi_cnst_src_exp.push_back(Fi); Ff_cnst_snk_exp.push_back(Ff);
    }


    if( fits.size() > 0 )
    {
      {

        map<double, AvgFit*, std::greater<double> > ordered = make_ordered_list( fits, fit_qual );
        minuit_fit_result                           best_cnst_one_exp = (ordered.begin()->second)->get_result();

      
        {
          vector<bool> active = (ordered.begin()->second)->get_active_data();
          int dt_count = 0;
          for(int i = 0 ; i < control.dt.size() ; i++){
          /* in the current case, the indexing of active_data corresponds to the count */
            int t_min_cnst_one_exp = 0; int t_max_cnst_one_exp = control.dt[i];

            while( !active[dt_count + t_min_cnst_one_exp] ){ t_min_cnst_one_exp++; if(t_min_cnst_one_exp > t_max_cnst_one_exp ){t_min_cnst_one_exp = 1; break;}; }
            while( !active[dt_count + t_max_cnst_one_exp] ){ t_max_cnst_one_exp--; if(t_min_cnst_one_exp > t_max_cnst_one_exp ){t_max_cnst_one_exp = 0; break;}; }
            t_min_cnst_one_exp--;t_max_cnst_one_exp++;
            tmin_cnst_one_exp.push_back(make_pair(control.dt[i], t_min_cnst_one_exp)); tmax_cnst_one_exp.push_back(make_pair(control.dt[i], t_max_cnst_one_exp));

            dt_count += control.dt[i]+1;
          }

        }
      
        F_cnst_one_exp = (best_cnst_one_exp.par_values.find("F"))->second;

      }

      if( fits_src_exp.size() > 0 )
      {

        map<double, AvgFit*, std::greater<double> > ordered_src = make_ordered_list( fits_src_exp, fit_qual );
        minuit_fit_result                           best_cnst_src_exp = (ordered_src.begin()->second)->get_result();
  
      
        for(int t = 0; t < control.dt.size(); t++){
          stringstream fi; fi << "Fi" << "_Dt" << control.dt[t];
          stringstream ei; ei << "Ei" << "_Dt" << control.dt[t];

          Fi_cnst_src_exp[t] = (best_cnst_src_exp.par_values.find(fi.str()))->second;
          Ei_cnst_src_exp[t] = (best_cnst_src_exp.par_values.find(ei.str()))->second;
        }


      }

      if( fits_snk_exp.size() > 0 )
      {

         map<double, AvgFit*, std::greater<double> > ordered_snk = make_ordered_list( fits_snk_exp, fit_qual );
         minuit_fit_result                           best_cnst_snk_exp = (ordered_snk.begin()->second)->get_result();

        for(int t = 0; t < control.dt.size(); t++){
          stringstream ff; ff << "Ff" << "_Dt" << control.dt[t];
          stringstream ef; ef << "Ef" << "_Dt" << control.dt[t];

          Ff_cnst_snk_exp[t] = (best_cnst_snk_exp.par_values.find(ff.str()))->second;
          Ef_cnst_snk_exp[t] = (best_cnst_snk_exp.par_values.find(ef.str()))->second;  
        }

      } 
    } 


    else
    {
      for(int i = 0 ; i < control.dt.size() ; i++){

      tmin_cnst.push_back(make_pair(control.dt[i], std::get<1>(control.tmin_max[i]) ));
      tmax_cnst.push_back(make_pair(control.dt[i], std::get<2>(control.tmin_max[i]) ));
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

          for(int t = 0; t < control.dt.size(); t++){
            stringstream ff; ff << "Ff" << "_Dt" << control.dt[t];
            stringstream ef; ef << "Ef" << "_Dt" << control.dt[t];

            stringstream fi; fi << "Fi" << "_Dt" << control.dt[t];
            stringstream ei; ei << "Ei" << "_Dt" << control.dt[t];

            Fi_cnst_src_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(fi.str(), Fi_cnst_src_exp[t]) );

            Ff_cnst_snk_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(ff.str(), Ff_cnst_snk_exp[t]) );

            Ei_cnst_src_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(ei.str(), Ei_cnst_src_exp[t]) );

            Ef_cnst_snk_exp[t].error *= 5.0; /* boost the error on the from factor */
            start_params.insert( make_pair(ef.str(), Ef_cnst_snk_exp[t]) );
          }

        }


     vector<bool> previous( accepted_data.size(), false );

      for(int i = 0; i < control.tmin_max.size(); i++){
        std::get<1>(control.tmin_max[i]) =  tmin_cnst_one_exp[i].second;
        std::get<2>(control.tmin_max[i]) =  tmax_cnst_one_exp[i].second;
      }

      std::vector< std::vector<pair<pair<int,int>,pair<int,int>>> > t_pairs  = get_all_t_ranges(control,0,1,1, range, count_dt);


        for(int i = 0 ; i < t_pairs.size() ; i++){
        cout << "----------------------------------------------------------------------------------" << endl;
          for(int dt_count = 0 ; dt_count < t_pairs[i].size() ; dt_count++)
          {

            PairIntAbscissa* x_t_low_dt = new PairIntAbscissa(t_pairs[i][dt_count].first);
            PairIntAbscissa* x_t_high_dt = new PairIntAbscissa(t_pairs[i][dt_count].second);

            x_tlow.push_back(x_t_low_dt); x_thigh.push_back(x_t_high_dt);


          } //looping over all dt 
      
          vector<bool> active_data = !(data.make_active_data()); //all false

          for(int j = 0 ; j < control.tmin_max.size() ; j++){

            cout << "two_exp_t__low ="<< *x_tlow[j] << "  | ";
            cout << "two_exp_t__high =" << *x_thigh[j] << endl;
            active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
          }

          active_data = active_data && accepted_data;

          if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
          {
          cout << "*** SUCCESSFUL SOURCE AND SINK EXP AVG FIT ***" << endl;
          cout << "----------------------------------------------------------------------------------" << endl;
          stringstream name; name << "c2";
          for(int k = 0 ; k < control.tmin_max.size() ; k++){
            name << "_Dt" << x_tlow[k]->get_x().first << "_" << x_tlow[k]->get_x().second << "-" << x_thigh[k]->get_x().second;
          }
    
          AvgFit* this_fit = new AvgFit( data, active_data, cnst_two_exp, start_params, minuit_controls, control.correlated, name.str() );
          fits.push_back( this_fit );
          previous = active_data;
    
          /*cout << "==================================================================================" << endl;
          cout << "*** 2 EXP , tmin = " << tlow  << " ***" << endl;
          cout << data.print_data(active_data) << endl;
          cout << this_fit->report() << endl;
          cout << "==================================================================================" << endl;*/
          }  
          x_tlow.clear(); x_thigh.clear();      
        } // next dt
        t_pairs.clear();
      }//end two exp fits

    }//end of !only exp loop

    
  }//end of !only cnst loop

  //************************************************************************************************

  //************************************************************************************************

    
  //************************************************************************************************
  /* done with avg fitting, select a best fit, do an ensem fit ... */

  cout << endl << "# fits = " << fits.size() << endl;

  map<double, AvgFit*, std::greater<double> > ordered_fits = make_ordered_list( fits, fit_qual );

  FitSelector fit_selector( fits, fit_qual, chisq_ndof_cutoff );

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
    cerr << "no successful ensemfits found, exiting.........." << endl; exit(1);
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
      int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);
      vector<pair<double,double>> active_timeslices;
      for(int i = 0; i < x.size(); i++){
        pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x();
        active_timeslices.push_back( make_pair(double(t_i.first),double(t_i.second)) );
      }

      vector<double> delt;
      for(int j = 0; j < control.dt.size(); j++){delt.push_back(double(control.dt[j])); }

      plot_three_point_timeslice_function_data t_y_yerr = plot_three_point_timeslice_ensem_function(ft, ensem_pars, active_timeslices, delt );
      
      for(int i = 0; i < t_y_yerr.central.size(); i++)
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

      for(int s = 1; s < control.dt.size(); s++){
        plot << "," << control.dt[s];
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
      
      { /* write active data */

        plot << "# active data" << endl;
        int count = data.get_active_data(ensem_fit_active_data, x, y, y_err);

        prev_dt = control.dt[0];
  
        for(int i = 0; i < x.size(); i++)
        {
          pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x(); 
          if(t_i.first != prev_dt){ plot << endl << endl << "# active data" << endl;}

          plot <<  t_i.first << " " << t_i.second  << " "; 
          plot << y[i] * F_param.ensem_mean << " ";
          plot << y_err[i] * F_param.ensem_mean << endl;

          prev_dt = t_i.first;
        }
          plot << endl << endl;

      }
      
      
      { /* write inactive data */
        plot << "# inactive data" << endl;
        int count = data.get_active_data(!ensem_fit_active_data, x, y, y_err);

        prev_dt = control.dt[0];
  
        for(int i = 0; i < x.size(); i++)
        {
          pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x();
          if(t_i.first != prev_dt){ plot << endl << endl << "# inactive data" << endl;}

          plot <<  t_i.first << " "<< t_i.second  << " ";  
          plot << y[i] * F_param.ensem_mean << " ";
          plot << y_err[i] * F_param.ensem_mean << endl;

          prev_dt = t_i.first;
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

        vector<double> delt;
        for(int j = 0; j < control.dt.size(); j++){delt.push_back(double(control.dt[j])); }

        plot_three_point_timeslice_function_data plot_fn = plot_three_point_timeslice_ensem_function(ft, ensem_pars, plot_t, delt );

        prev_dt = control.dt[0];

        for(int i = 0; i < plot_fn.central.size(); i++)
        {
          pair<double,double> t = plot_fn.central[i].first;
          if(t.first != prev_dt){ plot << endl << endl << "# ensem fit" << endl;}

          plot << t.first << "  " << t.second << "  "
          << plot_fn.lower[i].second   * F_param.ensem_mean << "  "
          << plot_fn.central[i].second * F_param.ensem_mean << "  "
          << plot_fn.upper[i].second   * F_param.ensem_mean << endl;

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
    
          for(int i = 0; i < plot_fn.size(); i++){
            pair<double,double> t = plot_fn[i].first;
            if(t.first != prev_dt){plot << endl << endl << "# avg fit: " << (*fit)->get_fit_name() << endl;}

            plot << t.first << "  " << t.second << "  " << plot_fn[i].second * F_param.ensem_mean << endl;

            prev_dt = t.first;
          }
    
          plot << endl << endl;
        }//next fit variation
      }
      
      out.plot_data = plot.str();
    }// end if plots

  }//end if success

  cout << "#################################################################################" << endl;
  
  /* clean up pointers */
  delete cnst; delete cnst_src_exp; delete cnst_two_exp; delete cnst_snk_exp;
  return out; 

};


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

  for(int i = 0 ; i < delt.size(); ++i)
  {
    int count = 0;
    int tmin = 0;
    int tmax = 0;

    for(auto it = t_values.begin() ; it != t_values.end(); ++it)
    {
      if(delt[i] == it->first){ count++; if(it->second < tmin){tmin = it->second;} if(it->second > tmax) {tmax = it->second;}}
      else{continue;}
    }
      out.n_points.push_back(count);
      out.tmin.push_back(tmin);
      out.tmax.push_back(tmax);
      out.dt.push_back(delt[i]);
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

//************************************************************************
// FUNCTION TO LOOP OVER ALL THE T_MINS AND T_MAXS
//************************************************************************


std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> get_all_t_ranges(fit_three_point_control control, bool slide, bool src, bool snk,
                  std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> range, int count_dt){


  std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> trange_bins_i;
  std::vector<std::vector<pair<pair<int,int>,pair<int,int>>>> trange_bins;

  pair<int,int> tmin, tmax;
    //cout << "in" << endl;

  if(slide && src && snk){


    for(int tlow =  std::get<1>(control.tmin_max[count_dt]); tlow >= control.tsrc[count_dt]; tlow--)
    { //cout << "...t_low" << tlow <<  endl;

      for(int thigh =  std::get<2>(control.tmin_max[count_dt]); thigh <= control.tsnk[count_dt]; thigh++)
      { //cout << "...t_high" << thigh << endl ;

        //if(tlow == control.tsrc[i] || tlow == control.tsnk[i] ){continue;}
        //if(thigh == control.tsrc[i] || thigh == control.tsnk[i] ){continue;}

        for(int t_slide_low = control.tsrc[count_dt]; t_slide_low <= control.tsnk[count_dt] - (thigh - tlow ); t_slide_low++){

          trange_bins_i = range;


          int t_slide_high = t_slide_low + thigh - tlow;

          //cout << "....t_slide_low" << t_slide_low << "....t_slide_high" << t_slide_high << endl;
          {

            tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),t_slide_low);
            tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),t_slide_high);

          if(trange_bins_i.size())
          {
            for(int i = 0; i < trange_bins_i.size(); i++){
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


    for(int thigh =  std::get<2>(control.tmin_max[count_dt]); thigh <= control.tsnk[count_dt]; thigh++)
    { //cout << "...t_high" << thigh << endl ;

      trange_bins_i = range;

      {

        tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),tlow);
        tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),thigh);

        if(trange_bins_i.size())
        {
          for(int i = 0; i < trange_bins_i.size(); i++){
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
  
    for(int tlow =  std::get<1>(control.tmin_max[count_dt]); tlow >= control.tsrc[count_dt]; tlow--)
    { //cout << "...t_low" << tlow <<  endl;

      trange_bins_i = range;

      {

        tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),tlow);
        tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),thigh);

        if(trange_bins_i.size())
        {
          for(int i = 0; i < trange_bins_i.size(); i++){
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

        
    for(int tlow =  std::get<1>(control.tmin_max[count_dt]); tlow >= control.tsrc[count_dt]; tlow--)
    { //cout << "...t_low" << tlow <<  endl;

      for(int thigh =  std::get<2>(control.tmin_max[count_dt]); thigh <= control.tsnk[count_dt]; thigh++)
      { //cout << "...t_high" << thigh << endl ;

        trange_bins_i = range;

        {

          tmin = make_pair(std::get<0>(control.tmin_max[count_dt]),tlow);
          tmax = make_pair(std::get<0>(control.tmin_max[count_dt]),thigh);

          if(trange_bins_i.size())
          {
            for(int i = 0; i < trange_bins_i.size(); i++){
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



