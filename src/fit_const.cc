/*
  read a prin corr ensem file and do a set of fits to it for various timeslice regions and one or two exponentials 
*/

#include "pair_int_abscissa_data_builder.h"
#include "fitting_lib/data_selectors.h"
#include "three_point_timeslice_fitting.h"


int main(int argc, char** argv){

  /* the factories are really required here, as we're only considering a very limited set 
     of functions and fit_qualities */
  //=================
  ThreePointFunctionEnv::registerAll();
  ThreePtFitQualityEnv::registerAll();
  //=================

  string error_msg = "fit_3pt_corr <Inifile> \n ";
  int n_args = 1;
  
  if( argc < n_args + 1 ){ cerr << error_msg; exit(1); }
  
  string xmlini;  {istringstream a(argv[1]); a >> xmlini;};  


  //============================
  //===== READ THE XML =====

  int num, dt, tmin, tmax, tsrc, tsnk, Nt_min;
  vector<string> filename;
  string file; 
  double noise, chisq_cut;
  string qual_type;

  vector<int> dt_v, tmin_v, tmax_v, tsrc_v, tsnk_v;
  vector<std::tuple<int,int,int>> tmin_max_v;

  XMLReader xml_in(xmlini);
  {
    read(xml_in,"/ThreeptIniParams/inputProps/NumdbFiles",num);


    for(int i = 1; i <= num; i++ ){ 
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/name", file);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/dt", dt);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/tsrc", tsrc);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/tsnk", tsnk);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/tmin", tmin);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/tmax", tmax);

      std::tuple<int,int,int> tmin_max = std::make_tuple(dt, tsrc,tsnk);

      filename.push_back(file);
      dt_v.push_back(dt); tmin_v.push_back(tmin); tmax_v.push_back(tmax); tsrc_v.push_back(tsrc); tsnk_v.push_back(tsnk);
      tmin_max_v.push_back(tmin_max);  
    }

    read(xml_in,"/ThreeptIniParams/FitProps/NtMin",Nt_min);
    read(xml_in,"/ThreeptIniParams/FitProps/fit_qual/type",qual_type);
    read(xml_in,"/ThreeptIniParams/FitProps/chisq_cutoff",chisq_cut);
    read(xml_in,"/ThreeptIniParams/FitProps/noise_cutoff",noise);


  }



  //============================
  //===== LOAD THE DATA =====

  vector<pair<int,int>> x_int; vector<EnsemReal> y_ensem; EnsemVectorReal tmp;


  for(int i = 0; i < num; i++ ){ 

    read(filename[i], tmp); 
    
    for(int t = 0; t < tmp.numElem(); t++){

      x_int.push_back(make_pair(dt_v[i],t));
      y_ensem.push_back( peekObs(tmp, t) );
      
    }
  }

  Data data; make_pair_int_abscissa_data(x_int, y_ensem, data );
  
  cout << "loaded data::" << endl << data.print_data() << endl;
  //============================

  //=================================
  //===== REMOVE UNWANTED POINTS ====
  vector<bool> remove_dt = !(data.make_active_data()); //all false
  {
    for(int i = 0; i < num; i++ ){ 
      Abscissa* x_dt = new PairIntAbscissa(make_pair(dt_v[i],dt_v[i]));
      remove_dt = remove_dt || !( data_at_x(data, x_dt ) );
      delete x_dt;
    }
  }
  
  vector<bool> remove_noisy = data_below_y_noise_ratio( data, noise);
  
  vector<bool> tmin_tmax = !(data.make_active_data()); //all false;
  {
    for(int i = 0; i < num; i++ ){ 
      Abscissa* x_tmin = new PairIntAbscissa(make_pair(dt_v[i], tmin_v[i] ));
      Abscissa* x_tmax = new PairIntAbscissa(make_pair(dt_v[i], tsnk_v[i] ));
      tmin_tmax = tmin_tmax || data_in_x_range( data, make_pair(x_tmin, x_tmax) );
      delete x_tmin; delete x_tmax;
    }
  }
  
  vector<bool> keep = ( remove_dt && remove_noisy && tmin_tmax );
  
  cout << "acceptable data::" << endl << data.print_data(keep) << endl;

  cout << "###############################" << endl;
  //===================================


  //============================
  //======= DO THE FITS ======== 
  fit_three_point_control control;

  control.tsrc = tsrc_v; control.tsnk = tsnk_v; control.Nt_min = Nt_min;
  control.plots = true ; control.dt = dt_v; control.tmin_max = tmin_max_v;

  //===================================

  if( count_active(keep) < control.dt.size() * Nt_min )
    { cerr << "fewer than " << Nt_min << " timeslices survive your restrictions, no fits will be acceptable, exiting ..." << endl; exit(1); }

  //===================================

  FitQuality* fit_qual;

  /* try a bit of xml magic to use the factory */
  /* for the fit_qual types used here, no xml content is required */
  {
    /* if you wanted to push into xml you could use something like this 
    XMLBufferWriter xml_write;
    push(xml_write, "Stuff");
    int pp = 10; write(xml_write, "fred", pp); // requires that a writer for the type of pp be known
    pop(xml_write);
    istringstream xml_s( xml_write.str() );
    XMLReader xml( xml_s );
    cout << xml_s.str() << endl;
    */
    
    XMLReader xml_in(xmlini);
    fit_qual = TheFitQualityFactory::Instance().createObject( qual_type, xml_in, "/ThreeptIniParams/FitProps/fit_qual/params");
  }
  
  fit_three_point_output output =  fit_three_point_corr( data, keep, control, fit_qual, chisq_cut);
  //============================


  
  
  //=============================
  //======= OUTPUT STUFF ========  

  /* write log to screen */
  cout << output.fit_summary << endl;
  
  /* write log to file */
  {
    stringstream s; s << xmlini << "_three_pt_exp_fit.log"; 
    ofstream out; out.open(s.str().c_str());
    out << output.fit_summary;
    out.close();
  }
  
  /* write mass ensem file */
  {  ostringstream outfile; outfile << xmlini << "_three_pt_exp_fit.jack";  write(outfile.str(), output.F );  }

  /* write mass variations to file */
  {
    stringstream s; s << xmlini << "_three_pt_exp_fit.syst"; 
    ofstream out; out.open(s.str().c_str());
    out << output.F_fit_variation;
    out.close();
  }

  // /* write plot data to file */
   {
     stringstream s; s << xmlini << "_three_pt_exp_fit.plot"; 
     ofstream out; out.open(s.str().c_str());
     out << output.plot_data;
     out.close();
   }
  
 

  /* for(map<int, pair<double,double> >::iterator it = output.best_data_description.begin(); it != output.best_data_description.end(); it++){
    cout << it->first << " " << (it->second).first << " " << (it->second).second << endl;
  }
  */

  /* can use this in 
vector<bool> reject_outliers( const Data& data,
                              const map< Abscissa*, pair<double,double> >& ref, 
                              double n_sigma );
  
to check for outliers 
  */
  
  delete fit_qual;
};
