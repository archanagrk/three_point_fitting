/*
  read a prin corr ensem file and do a set of fits to it for various timeslice regions and one or two exponentials 
*/

#include "fitting_lib/data_builders.h"
#include "fitting_lib/data_selectors.h"
#include "lib/z_timeslice_fitting.h"


int main(int argc, char** argv){

  /* the factories are really required here, as we're only considering a very limited set 
     of functions and fit_qualities */
  //=================
  FitQualityEnv::registerAll();
  ZFunctionEnv::registerAll();
  //=================

  string error_msg = "fit_zt <filename> <t0> <tmin> <tmax> <noise_cutoff> <fit_qual> <min_tslices> <chisq_cutoff>\n ";
  int n_args = 8;
  
  if( argc < n_args + 1 ){ cerr << error_msg; exit(1); }
  
  string filename;  {istringstream a(argv[1]); a >> filename;};  
  int t0;           {istringstream a(argv[2]); a >> t0;};
  int tmin;         {istringstream a(argv[3]); a >> tmin;}; 
  int tmax;         {istringstream a(argv[4]); a >> tmax;}; 
  double noise;     {istringstream a(argv[5]); a >> noise;};
  string qual_type; {istringstream a(argv[6]); a >> qual_type;};  
  int Nt_min;       {istringstream a(argv[7]); a >> Nt_min;}; 
  double chisq_cut; {istringstream a(argv[8]); a >> chisq_cut;};

  //============================
  //===== LOAD THE DATA =====
  EnsemVectorReal tmp; read(filename, tmp); 
  
  vector<int> x_int; vector<EnsemReal> y_ensem; 
  
  for(int t = 0; t < tmp.numElem(); t++){
    x_int.push_back(t);
    //cout << t;
    y_ensem.push_back( peekObs(tmp, t) );
  }
  
  Data data; make_int_abscissa_data(x_int, y_ensem, data );
  
  cout << "loaded data::" << endl << data.print_data() << endl;
  //============================

  //=================================
  //===== REMOVE UNWANTED POINTS ====
  vector<bool> remove_t0;
  {
    Abscissa* x_t0 = new IntAbscissa(t0);
    remove_t0 = !( data_at_x(data, x_t0 ) );
    delete x_t0;
  }
  
  vector<bool> remove_noisy = data_below_y_noise_ratio( data, noise);
  
  vector<bool> tmin_tmax;
  {
    Abscissa* x_tmin = new IntAbscissa(tmin - 1);
    Abscissa* x_tmax = new IntAbscissa(tmax + 1);
    tmin_tmax = data_in_x_range( data, make_pair(x_tmin, x_tmax) );
    delete x_tmin; delete x_tmax;
  }
  
  vector<bool> keep = ( remove_t0 && remove_noisy && tmin_tmax );
  
  cout << "acceptable data::" << endl << data.print_data(keep) << endl;

  if( count_active(keep) < Nt_min )
    { cerr << "fewer than " << Nt_min << " timeslices survive your restrictions, no fits will be acceptable, exiting ..." << endl; exit(1); }
  //===================================


  //============================
  //======= DO THE FITS ======== 
  fit_z_timeslice_control control;
  control.tmin = tmin; control.tmax = tmax; control.Nt_min = Nt_min;
  control.plots = true;

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
    
    XMLReader xml;
    fit_qual = TheFitQualityFactory::Instance().createObject( qual_type, xml, "/Stuff");
  }

  fit_z_timeslice_output output =  fit_z_timeslice( data, t0, keep, control, fit_qual, chisq_cut);
  //============================



  
  //=============================
  //======= OUTPUT STUFF ========  

  /* write log to screen */
  cout << output.fit_summary << endl;
  
  /* write log to file */
  {
    stringstream s; s << filename << "_z_exp_fit.log"; 
    ofstream out; out.open(s.str().c_str());
    out << output.fit_summary;
    out.close();
  }
  
  /* write mass ensem file */
  {  ostringstream outfile; outfile << filename << "_z_exp_fit.jack";  write(outfile.str(), output.m );  } 

  /* write mass variations to file */
  {
    stringstream s; s << filename << "_z_exp_fit.syst"; 
    ofstream out; out.open(s.str().c_str());
    out << output.m_fit_variation;
    out.close();
  }

  /* write plot data to file */
  {
    stringstream s; s << filename << "_z_exp_fit.plot"; 
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
