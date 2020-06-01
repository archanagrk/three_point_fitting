/*
  read a corr ensem file and do a set of fits to it for various timeslice regions with constant, one and two exponentials and find the best fit
*/

// three_pt_fit
#include "lib/pair_int_abscissa_data_builder.h"
#include "lib/three_point_timeslice_fitting.h"
#include "io/write_data.h"


// fitting lib
#include "fitting_lib/data_selectors.h"

using namespace SEMBLE;

int main(int argc, char** argv)
{


  map<string, vector< pair< pair<ENSEM::EnsemReal, ENSEM::EnsemReal> , ENSEM::EnsemReal >>> fqsq; 
  map<string, vector<ENSEM::EnsemReal> > favg, fmean;
  /* the factories are really required here, as we're only considering a very limited set 
     of functions and fit_qualities */

  //=================
  ThreePointFunctionEnv::registerAll();
  ThreePtFitQualityEnv::registerAll();
  //=================

  string error_msg = "fit_three_point_corr <Inifile> \n ";
  int n_args = 1;
  
  if( argc < n_args + 1 ){ cerr << error_msg; exit(1); }
  
  string xmlini;  {istringstream a(argv[1]); a >> xmlini;};  
  
  //============================
  //===== READ THE XML =====

  XMLReader xml_in(xmlini);
  XMLinidata_t xml_data = read_xml_ini(xml_in);
  
  //============================
  //===== INTERFACE TO EDB =====

  EDBdata_t edb_data = read_corrs(xml_data);

  //================================
  //===== MAKE THE DATA OBJECT =====


  /* Generates vector of data with the size  = number of Dt + 1. Each data[num_corr][i] has data for a particular Dt and all the Dts less than it.
    First ensemble fits are done on the smallest Dt. The range is fixed and then moves to the next Dt and does a combined fit with the fixed range
    smaller Dt and the curret Dt to fix the range of the current Dt. Then finally data[num_corr][size] contains all the data
    the the output of this fit is used - reduces the complexity from n^m to n*m */

  /* threading over corrs */
  int nthr = omp_get_max_threads(); 
  cout << "*** distributing corrs over " << nthr << " threads ***" << endl;
  #pragma omp parallel for num_threads(nthr) //schedule(dynamic)

  
  // START THE LOOP OVER THE CORRS IN THE EDB
  for(size_t num_corr = 0; num_corr < edb_data.x.size(); num_corr += 1)
  {

    int num_dt = edb_data.x[num_corr].size();
    vector<Data> data_rl(num_dt+1); vector<Data> data_im(num_dt+1);
    vector<pair<int,int> >x_total; vector<ENSEM::EnsemReal> y_ensem_total_rl; vector<ENSEM::EnsemReal> y_ensem_total_im;

    for(int i = 0; i < num_dt; i++ )
    {
      make_pair_int_abscissa_data(edb_data.x[num_corr][i], edb_data.y_ensem_rl[num_corr][i], data_rl[i] );
      make_pair_int_abscissa_data(edb_data.x[num_corr][i], edb_data.y_ensem_im[num_corr][i], data_im[i] );

      if(i == 0){x_total = edb_data.x[num_corr][i]; 
        y_ensem_total_rl = edb_data.y_ensem_rl[num_corr][i];
        y_ensem_total_im = edb_data.y_ensem_im[num_corr][i];}
      else
      {
        x_total.insert(std::end(x_total), std::begin(edb_data.x[num_corr][i]), std::end(edb_data.x[num_corr][i]));
        y_ensem_total_rl.insert(std::end(y_ensem_total_rl), std::begin(edb_data.y_ensem_rl[num_corr][i]), std::end(edb_data.y_ensem_rl[num_corr][i]));
        y_ensem_total_im.insert(std::end(y_ensem_total_im), std::begin(edb_data.y_ensem_im[num_corr][i]), std::end(edb_data.y_ensem_im[num_corr][i]));        
      }

    }

    make_pair_int_abscissa_data(x_total, y_ensem_total_rl, data_rl[num_dt] );
    make_pair_int_abscissa_data(x_total, y_ensem_total_im, data_im[num_dt] );
    
    cout << "loaded real(data)::" << endl << data_rl[num_dt].print_data() << endl;
    cout << "loaded imag(data)::" << endl << data_im[num_dt].print_data() << endl;


    //=================================

    //=================================
    //===== REMOVE UNWANTED POINTS ====


    vector<vector<bool> > keep_rl(num_dt+1);  vector<vector<bool> > keep_im(num_dt+1);

    vector<bool> remove_dt_all_rl = !(data_rl[num_dt].make_active_data()); //all false
    vector<bool> remove_noisy_all_rl = data_below_y_error(data_rl[num_dt], xml_data.noise); //removes all noisy data points depending on the absolute error in y
    vector<bool> tmin_tmax_all_rl = !(data_rl[num_dt].make_active_data()); //all false;

    vector<bool> remove_dt_all_im = !(data_im[num_dt].make_active_data()); //all false
    vector<bool> remove_noisy_all_im = data_below_y_error(data_im[num_dt], xml_data.noise); //removes all noisy data points depending on the absolute error in y
    vector<bool> tmin_tmax_all_im = !(data_im[num_dt].make_active_data()); //all false;
  
    for(int i = 0; i < num_dt; i++ )
    {

      vector<bool> remove_dt_rl = !(data_rl[i].make_active_data()); //remove the contact terms
      vector<bool> remove_dt_im = !(data_im[i].make_active_data()); //remove the contact terms
      {
        Abscissa* x_dt = new PairIntAbscissa(make_pair(xml_data.dt_v[i],xml_data.dt_v[i]));
        Abscissa* x_0 = new PairIntAbscissa(make_pair(xml_data.dt_v[i],0));

        vector<bool> remove_ends_rl = ( data_at_x(data_rl[i], x_dt ) ) || ( data_at_x(data_rl[i], x_0 ) ); 
        vector<bool> remove_ends_im = ( data_at_x(data_im[i], x_dt ) ) || ( data_at_x(data_im[i], x_0 ) ); 

        remove_dt_rl = remove_dt_rl || !(remove_ends_rl); remove_dt_im = remove_dt_im || !(remove_ends_im);

        vector<bool> remove_ends_all_rl = ( data_at_x(data_rl[num_dt], x_dt ) ) || ( data_at_x(data_rl[num_dt], x_0 ) ); 
        remove_dt_all_rl = remove_dt_all_rl || !(remove_ends_all_rl);

        vector<bool> remove_ends_all_im = ( data_at_x(data_im[num_dt], x_dt ) ) || ( data_at_x(data_im[num_dt], x_0 ) ); 
        remove_dt_all_im = remove_dt_all_im || !(remove_ends_all_im);

        delete x_dt;delete x_0;
      }
      
      vector<bool> remove_noisy_rl = data_below_y_error( data_rl[i], xml_data.noise); //removes all noisy data points depending on the absolute error in y
      vector<bool> remove_noisy_im = data_below_y_error( data_im[i], xml_data.noise); //removes all noisy data points depending on the absolute error in y     

      vector<bool> tmin_tmax_rl = !(data_rl[i].make_active_data()); //remove all the points below the tmin and above tmax provided by the user;
      vector<bool> tmin_tmax_im = !(data_im[i].make_active_data()); //remove all the points below the tmin and above tmax provided by the user;

      {
        Abscissa* x_tmin = new PairIntAbscissa(make_pair(xml_data.dt_v[i], xml_data.tmin_v[i] - 1 ));
        Abscissa* x_tmax = new PairIntAbscissa(make_pair(xml_data.dt_v[i], xml_data.tsnk_v[i] + 1 ));

        tmin_tmax_rl = tmin_tmax_rl || data_in_x_range( data_rl[i], make_pair(x_tmin, x_tmax) );
        tmin_tmax_all_rl = tmin_tmax_all_rl || data_in_x_range( data_rl[num_dt], make_pair(x_tmin, x_tmax) );

        tmin_tmax_im = tmin_tmax_im || data_in_x_range( data_im[i], make_pair(x_tmin, x_tmax) );
        tmin_tmax_all_im = tmin_tmax_all_im || data_in_x_range( data_im[num_dt], make_pair(x_tmin, x_tmax) );

        delete x_tmin; delete x_tmax;
      }

      keep_rl[i] = ( remove_dt_rl && remove_noisy_rl && tmin_tmax_rl );
      remove_dt_rl.clear(); remove_noisy_rl.clear(); tmin_tmax_rl.clear();

      keep_im[i] = ( remove_dt_im && remove_noisy_im && tmin_tmax_im );

      remove_dt_im.clear(); remove_noisy_im.clear(); tmin_tmax_im.clear();
    }

    keep_rl[num_dt] = ( remove_dt_all_rl && remove_noisy_all_rl && tmin_tmax_all_rl );
    keep_im[num_dt] = ( remove_dt_all_im && remove_noisy_all_im && tmin_tmax_all_im );    

    cout << "acceptable real(data)::" << endl << data_rl[num_dt].print_data(keep_rl[num_dt]) << endl;
    cout << "acceptable imag(data)::" << endl << data_im[num_dt].print_data(keep_im[num_dt]) << endl;

    cout << "###############################" << endl;


    remove_dt_all_rl.clear(); remove_noisy_all_rl.clear(); tmin_tmax_all_rl.clear();
    remove_dt_all_im.clear(); remove_noisy_all_im.clear(); tmin_tmax_all_im.clear();

    //===================================


    //============================
    //======= CONTROL FILE  ======== 

    
    vector<fit_three_point_control> control(num_dt+1);
    
    for(int i = 0; i < num_dt; i++)
    {

      control[i].tsrc.push_back(xml_data.tsrc_v[i]);
      control[i].tsnk.push_back(xml_data.tsnk_v[i]);
      control[i].dt.push_back(xml_data.dt_v[i]); 
      control[i].tmin_max.push_back(xml_data.tmin_max_v[i]); 

      control[num_dt].tsrc.push_back(xml_data.tsrc_v[i]);
      control[num_dt].tsnk.push_back(xml_data.tsnk_v[i]);
      control[num_dt].dt.push_back(xml_data.dt_v[i]); 
      control[num_dt].tmin_max.push_back(xml_data.tmin_max_v[i]);   

      control[i].Nt_min.push_back(xml_data.Nt_min_v[i]);
      control[i].Ei_min.push_back(xml_data.Ei_min_v[i]);
      control[i].Ef_min.push_back(xml_data.Ef_min_v[i]);
      control[i].Ei_max.push_back(xml_data.Ei_max_v[i]);
      control[i].Ef_max.push_back(xml_data.Ef_max_v[i]);

      control[num_dt].Nt_min.push_back(xml_data.Nt_min_v[i]);
      control[num_dt].Ei_min.push_back(xml_data.Ei_min_v[i]);
      control[num_dt].Ef_min.push_back(xml_data.Ef_min_v[i]);
      control[num_dt].Ei_max.push_back(xml_data.Ei_max_v[i]);
      control[num_dt].Ef_max.push_back(xml_data.Ef_max_v[i]);

      control[i].plots = true ;
    }
    
    control[num_dt].plots = true ;

    //===================================

    for(int i = 0; i < num_dt; i++ )
    {
      if( count_active(keep_rl[i]) < xml_data.Nt_min_v[i] )
        { cerr << "fewer than " << xml_data.Nt_min_v[i] << " timeslices survive your restrictions, no fits on the real part will be acceptable, exiting ..." << endl; exit(1); }
    }

    for(int i = 0; i < num_dt; i++ )
    {
      if( count_active(keep_im[i]) < xml_data.Nt_min_v[i] )
        { cerr << "fewer than " << xml_data.Nt_min_v[i] << " timeslices survive your restrictions, no fits on the imag part will be acceptable, exiting ..." << endl; exit(1); }
    }

    //===================================

    FitQuality* fit_qual;

    /* xml input for the factory */
    {

      XMLReader xml_in(xmlini);
      fit_qual = TheFitQualityFactory::Instance().createObject( xml_data.qual_type, xml_in, "/ThreeptIniParams/FitProps/fit_qual/params");
    }

    //============================
    //======= DIR NAMING ======== 

    /* Creates a subdirectory with the irrep names and the files for the respective fit is stored in this subdirectory */

    std::string name = edb_data.dir.find(num_corr)->second;
    std::stringstream ss;
    ss << "Q2_";
    ss << std::fixed << std::setprecision(6) << toDouble( mean(edb_data.pref.find(num_corr)->second.qsq) ); 
    ss << "/" + name;
    std::string path = SEMBLEIO::getPath() += ss.str();
    SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    cout << path << endl;

    //============================
    //======= DO THE FITS ======== 

    fit_three_point_output output_rl =  fit( data_rl, keep_rl, control, fit_qual, xml_data.chisq_cut, num_dt);
    fit_three_point_output output_im =  fit( data_im, keep_im, control, fit_qual, xml_data.chisq_cut, num_dt);

    //=============================
    //======= OUTPUT STUFF ======== 
    prefactor pf = edb_data.pref.find(num_corr)->second;
    print_output(path, name, pf, output_rl, output_im, xml_data.Zfile, fqsq, favg, fmean);

    data_rl.clear(); keep_rl.clear(); control.clear(); data_im.clear(); keep_im.clear(); 
    delete fit_qual;

    cout << "finisihed fitting " << num_corr + 1 << " corrs out of " <<  edb_data.x.size() << " corrs" << endl;
    cout << "====================================================================================" << endl;

  } // END OF THE LOOP OVER THE CORRS IN THE EBD

  cout << endl << "*******************************************************************************************************************************************" << endl;
  cout << " Finished fiiting all the three point correlation functions in the given edbs. " << endl << endl;


  /* Write the F(Q^2) vs Q^2 plot */

  if(xml_data.divkfac) print_extra_data(xmlini, xml_data.Zfile, fqsq, favg);

};



