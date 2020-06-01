/*
  read a prin corr ensem file and do a set of fits to it for various timeslice regions and one or two exponentials 
*/

#include "fitting_lib/data_builders.h"
#include "fitting_lib/data_selectors.h"

#include "lib/double_abscissa_data_builder.h"
#include "lib/formfac_qsq_fitting.h"
#include "io/read_formfac_qsq_data.h"



int main(int argc, char** argv){

  /* the factories are really required here, as we're only considering a very limited set 
     of functions and fit_qualities */
  //=================
  FitQualityEnv::registerAll();
  FFunctionEnv::registerAll();
  //=================

  string error_msg = "fit_ff_qsq <fit_ff_qsq_ini_xml>\n ";
  int n_args = 1;
  
  if( argc < n_args + 1 ){ cerr << error_msg; exit(1); }

  string xmlini;  {istringstream a(argv[1]); a >> xmlini;};   

  string qual_type, qsq_list, fit_form;
  double qsqmin, qsqmax, m_p, tcut, fbound;
  int Nq_min, n_max;
  double noise, chisq_cut;
  bool topt;

 XMLReader xml_in(xmlini);
  try
    {
      read(xml_in,"/FitFFQsqIniParams/QsqList",qsq_list);
      read(xml_in,"/FitFFQsqIniParams/qsqMin",qsqmin);
      read(xml_in,"/FitFFQsqIniParams/qsqMax",qsqmax);

      read(xml_in,"/FitFFQsqIniParams/mp",m_p);
      read(xml_in,"/FitFFQsqIniParams/FitParam/Nqmin",Nq_min);
      read(xml_in,"/FitFFQsqIniParams/FitParam/fitForm",fit_form);

      if(fit_form == "all" || fit_form == "ff_qsq_z_exp" ){
        string tmp;
        read(xml_in,"/FitFFQsqIniParams/FitParam/ZFitParam/topt", tmp);
        topt = ( tmp == "true"?true:false);
        read(xml_in,"/FitFFQsqIniParams/FitParam/ZFitParam/tcut",tcut);
        read(xml_in,"/FitFFQsqIniParams/FitParam/ZFitParam/fbound",fbound);
        read(xml_in,"/FitFFQsqIniParams/FitParam/ZFitParam/n_max",n_max);

      }

      read(xml_in,"/FitFFQsqIniParams/FitQual/fitCrit",qual_type);
      read(xml_in,"/FitFFQsqIniParams/FitQual/noise",noise);      
      read(xml_in,"/FitFFQsqIniParams/FitQual/cutoff",chisq_cut);      

    }
  catch( const string& error ){
    cerr << "Error reading input file : " << error << endl;
    }


  //============================
  //===== LOAD THE DATA =====

  map<double , EnsemReal> ff_qsq_list = readlist::make_ff_qsq(qsq_list);
  
  vector<double> x_double; vector<EnsemReal> y_ensem; 
  
  for(auto it = ff_qsq_list.begin(); it != ff_qsq_list.end(); it++){

    x_double.push_back(it->first);
    y_ensem.push_back(it->second);

  }
  
  Data data; make_double_abscissa_data(x_double, y_ensem, data );
  
  cout << "loaded data::" << endl << data.print_data() << endl;
  //============================

  //=================================
  //===== REMOVE UNWANTED POINTS ====

  vector<bool> remove_noisy = data_below_y_noise_ratio( data, noise);
  
  vector<bool> qsqmin_qsqmax;
  {
    Abscissa* x_qsqmin = new DoubleAbscissa(qsqmin - 0.00001);
    Abscissa* x_qsqmax = new DoubleAbscissa(qsqmax + 0.00001);
    qsqmin_qsqmax = data_in_x_range( data, make_pair(x_qsqmin, x_qsqmax) );
    delete x_qsqmin; delete x_qsqmax;
  }
  
  vector<bool> keep = ( !remove_noisy && qsqmin_qsqmax );
  
  cout << "acceptable data::" << endl << data.print_data(keep) << endl;

  if( count_active(keep) < Nq_min )
    { cerr << "fewer than " << Nq_min << " timeslices survive your restrictions, no fits will be acceptable, exiting ..." << endl; exit(1); }
  //===================================


  //============================
  //======= DO THE FITS ======== 
  fit_ff_qsq_control control;
  control.qsqmin = qsqmin; control.qsqmax = qsqmax; control.Nq_min = Nq_min;
  control.tcut = tcut; control.topt = topt; control.n_max = n_max; control.fbound = fbound;
  control.m_p = m_p; control.fit_form = fit_form;
  control.plots = true;

  FitQuality* fit_qual;

  /* xml input for the factory */
  {
    
    XMLReader xml;
    fit_qual = TheFitQualityFactory::Instance().createObject( qual_type, xml, "/Stuff");
  }

  fit_ff_qsq_output output =  fit_ff_qsq( data, keep, control, fit_qual, chisq_cut);
  //============================

  
  //=============================
  //======= OUTPUT STUFF ========  

  /* write log to screen */
  cout << output.fit_summary << endl;
  
  /* write log to file */
  {
    stringstream s; s << qsq_list << "_ff_qsq_fit.log"; 
    ofstream out; out.open(s.str().c_str());
    out << output.fit_summary;
    out.close();
  }
  
  /* write mass ensem file */
  {  ostringstream outfile; outfile << qsq_list << "_ff_qsq_fit.jack";  write(outfile.str(), output.msq );  } 

  /* write mass variations to file */
  {
    stringstream s; s << qsq_list << "_ff_qsq_fit.syst"; 
    ofstream out; out.open(s.str().c_str());
    out << output.msq_fit_variation;
    out.close();
  }

  /* write plot data to file */
  {
    stringstream s; s << qsq_list << "_ff_qsq_fit.plot"; 
    ofstream out; out.open(s.str().c_str());
    out << "## name= " << qsq_list << endl;
    out << output.plot_data;
    out.close();
  }
  
 
  
  delete fit_qual;
};
