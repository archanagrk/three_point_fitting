#include "three_point_timeslice_fitting.h"

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


  /* perform constant fits */
  {

    map<string, param_value> start_params;
    vector<Abscissa*> x_tlow;
    vector<Abscissa*> x_thigh;

    control.F_start = data.get_all_y_mean()[ int(control.dt[0]/2) ];
    control.F_err   = data.get_all_y_err()[ int(control.dt[0]/2) ]*2;
  
    param_value F(control.F_start, control.F_err); //the constant that would be the ff. Constraints?
    F.minos = control.minos;
    F.fixed = false;
    start_params.insert( make_pair("F", F ) ); 

  
    vector<bool> previous( accepted_data.size(), false );

  
    for(int i = 0 ; i < control.tmin_max.size() ; i++)
    {
      std::get<1>(control.tmin_max[i]) = int(std::get<0>(control.tmin_max[i])/2 - control.Nt_min/2);
      std::get<2>(control.tmin_max[i]) = int(std::get<0>(control.tmin_max[i])/2 + control.Nt_min/2);

      Abscissa* x_t_low_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), int(std::get<0>(control.tmin_max[i])/2 - control.Nt_min/2)) );
      Abscissa* x_t_high_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), int(std::get<0>(control.tmin_max[i])/2 + control.Nt_min/2)) );

      x_tlow.push_back(x_t_low_i); x_thigh.push_back(x_t_high_i);
      //delete x_t_low_i; delete x_t_high_i;

    }

    for(int i = 0 ; i < control.tmin_max.size() ; i++)
    {
      //cout << i << "...dt" << std::get<0>(control.tmin_max[i]);
      for(int tlow =  int(std::get<0>(control.tmin_max[i])/2 - control.Nt_min/2); tlow >= control.tsrc[i]; tlow--)
      {
        //cout << "...t_low" << tlow;
        for(int thigh =  int(std::get<0>(control.tmin_max[i])/2 + control.Nt_min/2); thigh <= control.tsnk[i]; thigh++)
        {
          //cout << "...t_high" << thigh ;
          if(tlow == control.tsrc[i] || tlow == control.tsnk[i] ){continue;}
          if(thigh == control.tsrc[i] || thigh == control.tsnk[i] ){continue;}

          for(int t_slide_low = control.tsrc[i]; t_slide_low < control.tsnk[i] - (thigh - tlow ); t_slide_low++){

            int t_slide_high = t_slide_low + thigh - tlow;

            cout << "....t_slide_low" << t_slide_low << "....t_slide_high" << t_slide_high << endl;
        
            vector<bool> active_data = !(data.make_active_data()); //all false
            {

              Abscissa* x_t_low_i = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_low));
              Abscissa* x_t_high_i = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_high));

              x_tlow[i] = x_t_low_i;
              x_thigh[i] = x_t_high_i;


              std::get<1>(control.tmin_max[i]) = t_slide_low;
              std::get<2>(control.tmin_max[i]) = t_slide_high;
              
              for(int j = 0 ; j < control.tmin_max.size() ; j++){
                active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
              }
              //delete x_t_low; delete x_t_high;
            }

            for(int m =0; m < control.tmin_max.size(); m++){

              cout << "low" << *x_tlow[m] << endl;
              cout << "high" << *x_thigh[m] << endl;
            }

            active_data = active_data && accepted_data;

            cout << "here" << count_active(active_data) << endl;

            if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
            {
              cout << "...SUCCESSFUL CONST AVG FIT!!!!..." << endl;
              stringstream name; name << "cnst";
              for(int k = 0 ; k < control.tmin_max.size() ; k++){
                name << "_dt" << std::get<0>(control.tmin_max[k]) << "_tmin" << std::get<1>(control.tmin_max[k]) << "_tmax" << std::get<2>(control.tmin_max[k]);
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
          } //next tslide 
        } // next tlow
      } //next thigh
    } // next dt
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
        int t_min_cnst = 0; int t_max_cnst = 0;

        while( !active[dt_count + t_min_cnst] ){ t_min_cnst++; }
        while( !active[dt_count + control.dt[i] - t_max_cnst] ){ t_max_cnst++; }
        tmin_cnst.push_back(make_pair(control.dt[i], t_min_cnst)); tmax_cnst.push_back(make_pair(control.dt[i],  control.dt[i] - t_max_cnst));

        dt_count += control.dt[i]+1;
      }

    }
  
    F_cnst = (best_cnst.par_values.find("F"))->second;
  }


  else{
    for(int i = 0 ; i < control.dt.size() ; i++){

      tmin_cnst.push_back(make_pair(control.dt[i],control.tsrc[i]));
      tmax_cnst.push_back(make_pair(control.dt[i],control.tsnk[i]));
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
      vector<Abscissa*> x_tlow;
      vector<Abscissa*> x_thigh;

      {
        F_cnst.error *= 5.0; /* boost the error on the mass */   //Constraints?
        start_params.insert( make_pair("F", F_cnst) );
        
        param_value dmi(2.0,2.0);
        dmi.low_limited = true; dmi.low_limit = 0.0; dmi.fixed = false;  // Constraints
        dmi.minos = control.minos;
        start_params.insert( make_pair("dmi", dmi ) );

        start_params.insert( make_pair("Fi", F_cnst) );

        
      }

    
      vector<bool> previous( accepted_data.size(), false );

      for(int i = 0 ; i < control.tmin_max.size() ; i++)
      {
        std::get<1>(control.tmin_max[i]) = tmin_cnst[i].second;
        std::get<2>(control.tmin_max[i]) = tmax_cnst[i].second;

        Abscissa* x_t_low_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), tmin_cnst[i].second ) );
        Abscissa* x_t_high_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), tmax_cnst[i].second ) );

        x_tlow.push_back(x_t_low_i); x_thigh.push_back(x_t_high_i);
        //delete x_t_low; delete x_t_high;

      }


      for(int i = 0 ; i < control.tmin_max.size() ; i++)
      {
        cout << i << "...dt" << std::get<0>(control.tmin_max[i]);
        for(int tlow =  tmin_cnst[i].second; tlow >= control.tsrc[i]; tlow--)
        {
          cout << "...t_low" << tlow;
          for(int thigh =  tmax_cnst[i].second; thigh <= control.tsnk[i]; thigh++)
          {
            cout << "...t_high" << thigh ;
            if(tlow == control.tsrc[i] || tlow == control.tsnk[i] ){continue;}
            if(thigh == control.tsrc[i] || thigh == control.tsnk[i] ){continue;}

            for(int t_slide_low = control.tsrc[i]; t_slide_low < control.tsnk[i] - (thigh - tlow ); t_slide_low++){

              int t_slide_high = t_slide_low + thigh - tlow;

              cout << "....t_slide_low" << t_slide_low << "....t_slide_high" << t_slide_high;
          
              vector<bool> active_data = !(data.make_active_data()); //all false
              {
                Abscissa* x_t_low = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_low));
                Abscissa* x_t_high = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_high));
                
                *x_tlow[i] = *x_t_low;
                *x_thigh[i] = *x_t_high;

                std::get<1>(control.tmin_max[i]) = t_slide_low;
                std::get<2>(control.tmin_max[i]) = t_slide_high;
                
                for(int j = 0 ; j < control.tmin_max.size() ; j++){
                  active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
                }
                delete x_t_low; delete x_t_high;
              }

              active_data = active_data && accepted_data;

              if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
              {
                cout << "...SUCCESSFUL SOURCE EXPONENTIAL AVG FIT!!!!..." << endl;
                stringstream name; name << "csrc";
                for(int k = 0 ; k < control.tmin_max.size() ; k++){
                  name << "_dt" << std::get<0>(control.tmin_max[k]) << "_tmin" << std::get<1>(control.tmin_max[k]) << "_tmax" << std::get<2>(control.tmin_max[k]);
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
            } //next tslide 
          } // next tlow
        } //next thigh
      } // next dt
    } //end of src_exp fits


    //************************************************************************************************

    //************************************************************************************************


    // snk_exp_fits
    {
      map<string, param_value> start_params;
      vector<Abscissa*> x_tlow;
      vector<Abscissa*> x_thigh;

      {
        F_cnst.error *= 5.0; /* boost the error on the mass */   //Constraints?
        start_params.insert( make_pair("F", F_cnst) );
        
        param_value dmf(2.0,2.0);
        dmf.low_limited = true; dmf.low_limit = 0.0;   // Constraints
        dmf.minos = control.minos;
        dmf.fixed = false;
        start_params.insert( make_pair("dmf", dmf ) );

        start_params.insert( make_pair("Ff", F_cnst) );
        
      }

    
      vector<bool> previous( accepted_data.size(), false );

      for(int i = 0 ; i < control.tmin_max.size() ; i++)
      {
        std::get<1>(control.tmin_max[i]) = tmin_cnst[i].second;
        std::get<2>(control.tmin_max[i]) = tmax_cnst[i].second;

        Abscissa* x_t_low_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), tmin_cnst[i].second ) );
        Abscissa* x_t_high_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), tmax_cnst[i].second ) );

        x_tlow.push_back(x_t_low_i); x_thigh.push_back(x_t_high_i);
        //delete x_t_low; delete x_t_high;

      }


      for(int i = 0 ; i < control.tmin_max.size() ; i++)
      {
        cout << i << "...dt" << std::get<0>(control.tmin_max[i]);
        for(int tlow =  tmin_cnst[i].second; tlow >= control.tsrc[i]; tlow--)
        {
          cout << "...t_low" << tlow;
          for(int thigh =  tmax_cnst[i].second; thigh <= control.tsnk[i]; thigh++)
          {
            cout << "...t_high" << thigh ;
            if(tlow == control.tsrc[i] || tlow == control.tsnk[i] ){continue;}
            if(thigh == control.tsrc[i] || thigh == control.tsnk[i] ){continue;}

            for(int t_slide_low = control.tsrc[i]; t_slide_low < control.tsnk[i] - (thigh - tlow ); t_slide_low++){

              int t_slide_high = t_slide_low + thigh - tlow;

              cout << "....t_slide_low" << t_slide_low << "....t_slide_high" << t_slide_high;
          
              vector<bool> active_data = !(data.make_active_data()); //all false
              {
                Abscissa* x_t_low = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_low));
                Abscissa* x_t_high = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_high));
                *x_tlow[i] = *x_t_low;
                *x_thigh[i] = *x_t_high;

                std::get<1>(control.tmin_max[i]) = t_slide_low;
                std::get<2>(control.tmin_max[i]) = t_slide_high;
                
                for(int j = 0 ; j < control.tmin_max.size() ; j++){
                  active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
                }
                delete x_t_low; delete x_t_high;
              }

              active_data = active_data && accepted_data;

              if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
              {
                cout << "...SUCCESSFUL SINK EXPONENTIAL AVG FIT!!!!..." << endl;
                stringstream name; name << "csnk";
                for(int k = 0 ; k < control.tmin_max.size() ; k++){
                  name << "_dt" << std::get<0>(control.tmin_max[k]) << "_tmin" << std::get<1>(control.tmin_max[k]) << "_tmax" << std::get<2>(control.tmin_max[k]);
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
            } //next tslide 
          } // next tlow
        } //next thigh
      } // next dt
    } //end of snk_exp fits


    //************************************************************************************************
  
    /* pick the best cnst_src_exp fit and cnst_snk_exp to choose a start mass for cnst_two_exp fits */

    vector<std::pair<int,int>>  tmin_cnst_one_exp;
    vector<std::pair<int,int>>  tmax_cnst_one_exp; 

    param_value F_cnst_one_exp(control.F_start, control.F_err);
    F_cnst_one_exp.minos = control.minos;
    F_cnst_one_exp.fixed = false;

    param_value dmi_cnst_src_exp(2.0,2.0); 
    dmi_cnst_src_exp.low_limited = true; dmi_cnst_src_exp.low_limit = 0.0;   // Constraints
    dmi_cnst_src_exp.minos = control.minos;
    dmi_cnst_src_exp.fixed = false;

    param_value Fi_cnst_src_exp(control.F_start, control.F_err);
    Fi_cnst_src_exp.minos = control.minos;
    Fi_cnst_src_exp.fixed = false;

    param_value dmf_cnst_snk_exp(2.0,2.0);
    dmf_cnst_snk_exp.low_limited = true; dmf_cnst_snk_exp.low_limit = 0.0;   // Constraints
    dmf_cnst_snk_exp.minos = control.minos;
    dmf_cnst_snk_exp.fixed = false;

    param_value Ff_cnst_snk_exp(control.F_start, control.F_err);
    Ff_cnst_snk_exp.minos = control.minos; 
    Ff_cnst_snk_exp.fixed = false;



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
            int t_min_cnst_one_exp = 0; int t_max_cnst_one_exp = 0;

            while( !active[dt_count + t_min_cnst_one_exp] ){ t_min_cnst_one_exp++; }
            while( !active[dt_count + control.dt[i] - t_max_cnst_one_exp] ){ t_max_cnst_one_exp++; }
            tmin_cnst_one_exp.push_back(make_pair(control.dt[i], t_min_cnst_one_exp)); tmax_cnst_one_exp.push_back(make_pair(control.dt[i],  control.dt[i] - t_max_cnst_one_exp));

            dt_count += control.dt[i]+1;
          }

        }
      
        F_cnst_one_exp = (best_cnst_one_exp.par_values.find("F"))->second;

      }

      if( fits_src_exp.size() > 0 )
      {

        map<double, AvgFit*, std::greater<double> > ordered_src = make_ordered_list( fits_src_exp, fit_qual );
        minuit_fit_result                           best_cnst_src_exp = (ordered_src.begin()->second)->get_result();
  
      
        Fi_cnst_src_exp = (best_cnst_src_exp.par_values.find("Fi"))->second;
        dmi_cnst_src_exp = (best_cnst_src_exp.par_values.find("dmi"))->second;


      }

      if( fits_snk_exp.size() > 0 )
      {

         map<double, AvgFit*, std::greater<double> > ordered_snk = make_ordered_list( fits_snk_exp, fit_qual );
         minuit_fit_result                           best_cnst_snk_exp = (ordered_snk.begin()->second)->get_result();

    
         Ff_cnst_snk_exp = (best_cnst_snk_exp.par_values.find("Ff"))->second;
         //dmf_cnst_snk_exp = (best_cnst_snk_exp.par_values.find("dmf"))->second;  

      }  
    } 


    else
    {
      for(int i = 0 ; i < control.dt.size() ; i++){

        tmin_cnst_one_exp.push_back(make_pair(control.dt[i],control.tsrc[i]));
        tmax_cnst_one_exp.push_back(make_pair(control.dt[i],control.tsnk[i]));
      }
    } 

  
    //************************************************************************************************
    /* perform two exp fits */
    if(!control.only_one_exp)
    {

      {

        map<string, param_value> start_params;
        vector<Abscissa*> x_tlow;
        vector<Abscissa*> x_thigh;

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

        for(int i = 0 ; i < control.tmin_max.size() ; i++)
        {
          std::get<1>(control.tmin_max[i]) = tmin_cnst_one_exp[i].second;
          std::get<2>(control.tmin_max[i]) = tmax_cnst_one_exp[i].second;

          Abscissa* x_t_low_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), tmin_cnst_one_exp[i].second ) );
          Abscissa* x_t_high_i = new PairIntAbscissa(make_pair( std::get<0>(control.tmin_max[i]), tmax_cnst_one_exp[i].second ) );

          x_tlow.push_back(x_t_low_i); x_thigh.push_back(x_t_high_i);
          //delete x_t_low; delete x_t_high;

        }


        for(int i = 0 ; i < control.tmin_max.size() ; i++)
        {
          cout << i << "...dt" << std::get<0>(control.tmin_max[i]);
          for(int tlow =  tmin_cnst_one_exp[i].second; tlow >= control.tsrc[i]; tlow--)
          {
            cout << "...t_low" << tlow;
            for(int thigh =  tmax_cnst_one_exp[i].second; thigh <= control.tsnk[i]; thigh++)
            {
              cout << "...t_high" << thigh ;
              if(tlow == control.tsrc[i] || tlow == control.tsnk[i] ){continue;}
              if(thigh == control.tsrc[i] || thigh == control.tsnk[i] ){continue;}

              for(int t_slide_low = control.tsrc[i]; t_slide_low < control.tsnk[i] - (thigh - tlow ); t_slide_low++){

                int t_slide_high = t_slide_low + thigh - tlow;

                cout << "....t_slide_low" << t_slide_low << "....t_slide_high" << t_slide_high;
            
                vector<bool> active_data = !(data.make_active_data()); //all false
                {
                  Abscissa* x_t_low = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_low));
                  Abscissa* x_t_high = new PairIntAbscissa(make_pair(std::get<0>(control.tmin_max[i]), t_slide_high));
                  *x_tlow[i] = *x_t_low;
                  *x_thigh[i] = *x_t_high;

                  std::get<1>(control.tmin_max[i]) = t_slide_low;
                  std::get<2>(control.tmin_max[i]) = t_slide_high;
                  
                  for(int j = 0 ; j < control.tmin_max.size() ; j++){
                    active_data = active_data || data_in_x_range( data, make_pair(x_tlow[j],x_thigh[j]) );
                  }
                  delete x_t_low; delete x_t_high;
                }

                active_data = active_data && accepted_data;

                if( (count_active(active_data) >= control.dt.size() * control.Nt_min) && (active_data != previous) )
                {
                  cout << "...SUCCESSFUL SOURCE AND SINK EXP AVG FIT!!!!..." << endl;
                  stringstream name; name << "c2";
                  for(int k = 0 ; k < control.tmin_max.size() ; k++){
                    name << "_dt" << std::get<0>(control.tmin_max[k]) << "_tmin" << std::get<1>(control.tmin_max[k]) << "_tmax" << std::get<2>(control.tmin_max[k]);
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
              } //next tslide 
            } // next tlow
          } //next thigh
        } // next dt

      
      }//end two exp fits


    }//end of !only exp loop

    
  }//end of !only cnst loop

  //************************************************************************************************

  //************************************************************************************************

    
  //************************************************************************************************
  /* done with avg fitting, select a best fit, do an ensem fit ... */

  cout << endl << "fit" << fits.size() << endl;
  cout << endl << "carrying out ensem fits............." << endl;

  FitSelector fit_selector( fits, fit_qual, chisq_ndof_cutoff );

  //************************************************************************************************

  cout << endl << "#################################################################################" << endl;
  
  cout << endl << "success=" << fit_selector.ensem_success() << endl;
  //cout << fit_selector.get_summary();

    cout << endl << "#################################################################################" << endl;
  
  //************************************************************************************************
  /* build the output object */

  out.success     = fit_selector.ensem_success();

  if( !out.success  || control.long_log ){
    out.fit_long_log = fit_selector.get_long_log();
  }

  if( !out.success ){
    cerr << "no successful ensemfits found, exiting.........." << endl; exit(1);
  }


  out.fit_summary = fit_selector.get_summary();


  if(out.success)
  {
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

      /*for(int j = 0; j < control.dt.size(); j++){
        plot << "## tmin= " << control.tsrc[j] << " tmax= " << control.tsnk[j] << " nfits = " << accepted_fits.size() << endl;
      }*/
      

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
          pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x(); plot <<  t_i.first << " " << t_i.second  << " "; 
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
          pair<int,int> t_i = ( static_cast<PairIntAbscissa*>(x[i]) )->get_x(); plot <<  t_i.first << " "<< t_i.second  << " ";  
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

        vector<double> delt;
        for(int j = 0; j < control.dt.size(); j++){delt.push_back(double(control.dt[j])); }

        plot_three_point_timeslice_function_data plot_fn = plot_three_point_timeslice_ensem_function(ft, ensem_pars, plot_t, delt );

        for(int i = 0; i < plot_fn.central.size(); i++)
        {
          pair<double,double> t = plot_fn.central[i].first;
          plot << t.first << "  " << t.second << "  "
          << plot_fn.lower[i].second   * F_param.ensem_mean << "  "
          << plot_fn.central[i].second * F_param.ensem_mean << "  "
          << plot_fn.upper[i].second   * F_param.ensem_mean << endl;
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
            plot << t.first << "  " << t.second << "  " << plot_fn[i].second * F_param.ensem_mean << endl;
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
  cout << "ncfgs............" << ncfgs; 
  
  plot_three_point_timeslice_function_data out;
  out.dt = delt;

  for(int i = 0 ; i < delt.size(); ++i)
  {
    cout << i << "dt..." << delt[i];
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
    cout << "ensm";
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

