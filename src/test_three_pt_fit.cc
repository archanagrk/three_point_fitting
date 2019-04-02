/*
  read a prin corr ensem file and do a set of fits to it for various timeslice regions and one or two exponentials 
*/

#include "pair_int_abscissa_data_builder.h"
#include "fitting_lib/data_selectors.h"
#include "three_point_timeslice_fitting.h"
#include "fitting_lib/tools.h" 


/* generate some fake data */
void make_three_points( const vector<int>& Dts, int ncfgs, Data& data){


  mapstringdouble par_values;
  par_values["F"] = 1.0;

  vector< Abscissa* > xs;
  vector<double> values, errs;

  int seed = 999;

  double noise_level = 0.1; // default to 0.1
  
  vector<int> dt;
  
  for( vector<int>::const_iterator Dt = Dts.begin(); Dt != Dts.end(); Dt++){

    dt.push_back(*Dt);
    
    par_values["Fi_Dt" + to_string(*Dt)] = 0.25;
    par_values["Ei_Dt" + to_string(*Dt)] = 0.2;

    par_values["Ff_Dt" + to_string(*Dt)] = -0.65;
    par_values["Ef_Dt" + to_string(*Dt)] = 0.14;

    for(int t = 0; t <= *Dt; t++){
      
      ThreePointtimesliceCorrNExp f(1,1, dt);
      xs.push_back( new PairIntAbscissa(make_pair(*Dt, t) ) );
      double y = f( *(xs.back()), par_values );
      
      values.push_back(y);
      errs.push_back( y * noise_level * double(*Dt) / double(Dts[0]) ); 
    }   
  }

  int n_data = values.size();

  /* just go with nearly uncorrelated for now */
  itpp::mat corr = 0.02*itpp::ones(n_data, n_data);
  for(int i = 0; i < n_data; i++){
    corr(i,i) = 1.0;
  }

  vector<EnsemReal> ys = generate_correlated_ensembles( values, errs, corr, ncfgs, seed, false);

  data.make(xs, ys);
};


int main(int argc, char** argv){

  /* the factories  */
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

  int num, dt, j, tmin, tmax, tsrc, tsnk, Nt_min;
  vector<string> filename;
  string file; 
  double noise, chisq_cut;
  string qual_type;

  vector<int> dt_v, tmin_v, tmax_v, tsrc_v, tsnk_v;
  vector<std::tuple<int,int,int>> tmin_max_v;

  XMLReader xml_in(xmlini);
  {
    read(xml_in,"/ThreeptIniParams/inputProps/NumdbFiles",num);

    /* Have to sort the elems if not already descending order in Dt */
    std::map<int,int> sort_dt;
    for(int i = 1; i <= num; i++ ){ 
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/dt", dt);
      sort_dt.insert(make_pair(dt,i));
    }


    for(std::map<int,int>::iterator it=sort_dt.begin(); it!=sort_dt.end(); ++it){ 
      j = it->second;
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/dt", dt);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tsrc", tsrc);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tsnk", tsnk);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tmin", tmin);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tmax", tmax);

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
  //===== LOAD THE FAKE DATA =====

  /* Generates vector of data with the size  = number of Dt. Each data[i] has data for a particular Dt and all the Dts less than it.
   First ensemble fits are done on the smallest Dt. The range is fixed and then moves to the next Dt and does a combined fit with the fixed range
   smaller Dt and the curret Dt to fix the range of the current Dt. Then finally data[num - 1] contains all the data the the output of this fit is used.
   reduces the complexity from n^m to n*m */

  vector<Data> data(num); 
  vector<int> Dts;

  for(int i = 0; i < num; i++ ){ 

    Dts.push_back(dt_v[i]);

    int ncfgs = 400;
    int seed  = 999;
    
    make_three_points( Dts, ncfgs, data[i] );
  }

  cout << "loaded data::" << endl << data[num-1].print_data() << endl;


  //============================

  //=================================
  //===== REMOVE UNWANTED POINTS ====


  vector<vector<bool> > keep(num);

  for(int i = 0; i < num; i++ ){

    vector<bool> remove_dt = !(data[i].make_active_data()); //all false
    {
      Abscissa* x_dt = new PairIntAbscissa(make_pair(dt_v[i],dt_v[i]));
      remove_dt = remove_dt || !( data_at_x(data[i], x_dt ) );
      delete x_dt;
    }
    
    vector<bool> remove_noisy = data_below_y_noise_ratio( data[i], noise);
    
    vector<bool> tmin_tmax = !(data[i].make_active_data()); //all false;

    {
      Abscissa* x_tmin = new PairIntAbscissa(make_pair(dt_v[i], tmin_v[i] ));
      Abscissa* x_tmax = new PairIntAbscissa(make_pair(dt_v[i], tsnk_v[i] ));
      tmin_tmax = tmin_tmax || data_in_x_range( data[i], make_pair(x_tmin, x_tmax) );
      delete x_tmin; delete x_tmax;
    }

    vector<bool> tmp = ( remove_dt && remove_noisy && tmin_tmax );

    if(i > 0){
      tmp.erase(tmp.begin(), tmp.begin() +  keep[i-1].size());
      keep[i].reserve( keep[i-1].size() + tmp.size() ); // preallocate memory
      keep[i].insert( keep[i].end(), keep[i-1].begin(), keep[i-1].end() );
      keep[i].insert( keep[i].end(), tmp.begin(), tmp.end() );
      tmp.clear();
    }

    else{keep[i] = tmp;}

  }
  
  cout << "acceptable data::" << endl << data[num-1].print_data(keep[num-1]) << endl;

  cout << "###############################" << endl;
  //===================================


  //============================
  //======= DO THE FITS ======== 
  vector<fit_three_point_control> control(num);

  for(int i = 0; i < num; i++){

    if(i > 0){  
      control[i].tsrc = control[i-1].tsrc;
      control[i].tsnk = control[i-1].tsnk;
      control[i].dt = control[i-1].dt;
      control[i].tmin_max = control[i-1].tmin_max;
    }

    control[i].tsrc.push_back(tsrc_v[i]);
    control[i].tsnk.push_back(tsnk_v[i]);
    control[i].dt.push_back(dt_v[i]); 
    control[i].tmin_max.push_back(tmin_max_v[i]);   


    control[i].Nt_min = Nt_min;
    control[i].plots = true ;
  }

  //===================================

  for(int i = 0; i < num; i++ ){

    if( count_active(keep[i]) < control[i].dt.size() * Nt_min )
      { cerr << "fewer than " << Nt_min << " timeslices survive your restrictions, no fits will be acceptable, exiting ..." << endl; exit(1); }
}

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
  
  fit_three_point_output output =  fit( data, keep, control, fit_qual, chisq_cut, num);
  //============================


  
  
  //=============================
  //======= OUTPUT STUFF ========  

  
  /* write log to file */
  {
    stringstream s; s << xmlini << "_three_pt_test.log"; 
    ofstream out; out.open(s.str().c_str());
    out << output.fit_summary;
    out.close();
  }
  
  /* write mass ensem file */
  {  ostringstream outfile; outfile << xmlini << "_three_pt_test.jack";  write(outfile.str(), output.F );  }

  /* write mass variations to file */
  {
    stringstream s; s << xmlini << "_three_pt_test.syst"; 
    ofstream out; out.open(s.str().c_str());
    out << output.F_fit_variation;
    out.close();
  }

  // /* write plot data to file */
   {
     stringstream s; s << xmlini << "_three_pt_test.plot"; 
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
