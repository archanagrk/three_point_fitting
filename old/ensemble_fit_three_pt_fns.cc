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


  int num, dt, j, tmin, tmax, tsrc, tsnk, Nt_min;
  vector<string> filename;
  string file; 
  double noise, chisq_cut;
  string qual_type;

  vector<int> dt_v, tmin_v, tmax_v, tsrc_v, tsnk_v, Nt_min_v;
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
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/name", file);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/dt", dt);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tsrc", tsrc);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tsnk", tsnk);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tmin", tmin);
      read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/tmax", tmax);
      read(xml_in,"/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/NtMin",Nt_min);

      std::tuple<int,int,int> tmin_max = std::make_tuple(dt, tsrc,tsnk);

      filename.push_back(file);
      dt_v.push_back(dt); tmin_v.push_back(tmin); tmax_v.push_back(tmax); tsrc_v.push_back(tsrc); tsnk_v.push_back(tsnk);
      tmin_max_v.push_back(tmin_max); Nt_min_v.push_back(Nt_min); 
    }

    read(xml_in,"/ThreeptIniParams/FitProps/fit_qual/type",qual_type);
    read(xml_in,"/ThreeptIniParams/FitProps/chisq_cutoff",chisq_cut);
    read(xml_in,"/ThreeptIniParams/FitProps/noise_cutoff",noise);


  }


  //============================
  //===== LOAD THE DATA =====
  
  /* Generates vector of data with the size  = number of Dt. Each data[i] has data for a particular Dt and all the Dts less than it.
   First ensemble fits are done on the smallest Dt. The range is fixed and then moves to the next Dt and does a combined fit with the fixed range
   smaller Dt and the curret Dt to fix the range of the current Dt. Then finally data[num - 1] contains all the data the the output of this fit is used.
   reduces the complexity from n^m to n*m */


  vector<pair<int,int>> x; vector<EnsemReal> y_ensem; EnsemVectorReal tmp;
  vector<pair<int,int>> x_all; vector<EnsemReal> y_ensem_all; 
  vector<Data> data(num+1); 


  for(int i = 0; i < num; i++ ){ 

    read(filename[i], tmp); 
    
    for(int t = 0; t < tmp.numElem(); t++){

      x.push_back(make_pair(dt_v[i],t));
      y_ensem.push_back( peekObs(tmp, t) );

      x_all.push_back(make_pair(dt_v[i],t));
      y_ensem_all.push_back( peekObs(tmp, t) );
      
    }
    
    make_pair_int_abscissa_data(x, y_ensem, data[i] );
    x.clear(); y_ensem.clear();
    
  }

  make_pair_int_abscissa_data(x_all, y_ensem_all, data[num] );
  cout << "loaded data::" << endl << data[num].print_data() << endl;
  //============================

  //=================================
  //===== REMOVE UNWANTED POINTS ====


  vector<vector<bool> > keep(num+1);

  vector<bool> remove_dt_all = !(data[num].make_active_data()); //all false
  vector<bool> remove_noisy_all = data_below_y_noise_ratio( data[num], noise); //removes all noisy data points
  vector<bool> tmin_tmax_all = !(data[num].make_active_data()); //all false;


  for(int i = 0; i < num; i++ ){

    vector<bool> remove_dt = !(data[i].make_active_data()); //all false
    {
      Abscissa* x_dt = new PairIntAbscissa(make_pair(dt_v[i],dt_v[i]));
      remove_dt = remove_dt || !( data_at_x(data[i], x_dt ) );
      remove_dt_all = remove_dt_all || !( data_at_x(data[num], x_dt ) );
      delete x_dt;
    }
    
    vector<bool> remove_noisy = data_below_y_noise_ratio( data[i], noise);
    
    vector<bool> tmin_tmax = !(data[i].make_active_data()); //all false;

    {
      Abscissa* x_tmin = new PairIntAbscissa(make_pair(dt_v[i], tmin_v[i] - 1 ));
      Abscissa* x_tmax = new PairIntAbscissa(make_pair(dt_v[i], tsnk_v[i] + 1 ));
      tmin_tmax = tmin_tmax || data_in_x_range( data[i], make_pair(x_tmin, x_tmax) );
      tmin_tmax_all = tmin_tmax_all || data_in_x_range( data[num], make_pair(x_tmin, x_tmax) );
      delete x_tmin; delete x_tmax;
    }

    keep[i] = ( remove_dt && remove_noisy && tmin_tmax );

  }

  keep[num] = ( remove_dt_all && remove_noisy_all && tmin_tmax_all );
  
  cout << "acceptable data::" << endl << data[num].print_data(keep[num]) << endl;

  cout << "###############################" << endl;
  //===================================


  //============================
  //======= DO THE FITS ======== 
  vector<fit_three_point_control> control(num+1);
  
  for(int i = 0; i < num; i++){

    control[i].tsrc.push_back(tsrc_v[i]);
    control[i].tsnk.push_back(tsnk_v[i]);
    control[i].dt.push_back(dt_v[i]); 
    control[i].tmin_max.push_back(tmin_max_v[i]); 

    control[num].tsrc.push_back(tsrc_v[i]);
    control[num].tsnk.push_back(tsnk_v[i]);
    control[num].dt.push_back(dt_v[i]); 
    control[num].tmin_max.push_back(tmin_max_v[i]);   


    control[i].Nt_min.push_back(Nt_min_v[i]);
    control[num].Nt_min.push_back(Nt_min_v[i]);
    control[i].plots = true ;
  }
  
  control[num].plots = true ;

  //===================================

  for(int i = 0; i < num; i++ ){

    if( count_active(keep[i]) < Nt_min_v[i] )
      { cerr << "fewer than " << Nt_min_v[i] << " timeslices survive your restrictions, no fits will be acceptable, exiting ..." << endl; exit(1); }
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
    stringstream s; s << xmlini << "_three_pt_fit.log"; 
    ofstream out; out.open(s.str().c_str());
    out << output.fit_summary;
    out.close();
  }
  
  /* write mass ensem file */
  {  ostringstream outfile; outfile << xmlini << "_three_pt_fit.jack";  write(outfile.str(), output.F );  }

  /* write mass variations to file */
  {
    stringstream s; s << xmlini << "_three_pt_fit.syst"; 
    ofstream out; out.open(s.str().c_str());
    out << output.F_fit_variation;
    out.close();
  }

  // /* write plot data to file */
   {
     stringstream s; s << xmlini << "_three_pt_fit.plot"; 
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
