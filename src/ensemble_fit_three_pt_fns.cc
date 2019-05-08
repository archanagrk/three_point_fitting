/*
  read a prin corr ensem file and do a set of fits to it for various timeslice regions and one or two exponentials 
*/

#include "pair_int_abscissa_data_builder.h"
#include "three_point_timeslice_fitting.h"
#include "edb_reader.h"
#include "key_struct.h"

//libraries for making a directory
#include <sys/types.h>
#include <sys/stat.h>
#include"semble/semble_file_management.h"


#include "fitting_lib/data_selectors.h"

#include "hadron/ensem_filenames.h"
#include "semble/semble_meta.h"



int main(int argc, char** argv)
{

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
  double noise, chisq_cut, Ei, Ef;
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
    read(xml_in,"/ThreeptIniParams/FitProps/Ei",Ei);      
    read(xml_in,"/ThreeptIniParams/FitProps/Ef",Ef);


  }

  //============================
  //===== INTERFACE TO EDB =====


  // Read all the correlators from an edb file the edbs need to be the same irreps they are not verified. User has to make sure
  std::map<key_struct::KeyHadronSUNNPartNPtIrrep_t, int> key;
  std::map<int,string> dir;
  int key_count = 0;
  int pn;


  vector<vector<vector<pair<int,int>>> >x; vector<vector< vector<ENSEM::EnsemReal> >> y_ensem;
  vector<pair<int,int> >x_total; vector<ENSEM::EnsemReal> y_ensem_total;
  vector<pair<int,int>> x_t; vector<ENSEM::EnsemReal> y_t_ensem; 
  vector<vector<pair<int,int>>> x_vec; vector< std::vector<ENSEM::EnsemReal> > y_vec;  

  for(int i = 0; i < num; i++ )
  {

    std::cout << "read file = " << filename[i] << std::endl;

    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > corr = edb_reader::getCorrs(filename[i]);

    std::cout << "end of file =  " << filename[i] << std::endl;

    //============================
    //===== LOAD THE DATA =====

    /* Generates vector of data with the size  = number of Dt. Each data[pn][i] has data for a particular Dt and all the Dts less than it.
      First ensemble fits are done on the smallest Dt. The range is fixed and then moves to the next Dt and does a combined fit with the fixed range
      smaller Dt and the curret Dt to fix the range of the current Dt. Then finally data[pn][size - 1] contains all the data the the output of this fit is used.
      reduces the complexity from n^m to n*m */
      
    for(ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> >::const_iterator  in = corr.begin(); in != corr.end(); ++in) 
    {
      
      //name of directory
      std::string name = Hadron::ensemFileName(in->first); removeSubstrs(name, ".dat");

      int np = in->first.npoint.size();

      vector<XMLArray::Array<int>> key_mom(np); vector<int> key_row(np); vector<string> key_irrep(np);


      for(int k = 1; k <= np; k++ ) //Array1DO has 1-based indices
      {
        key_row[k-1]   =  in->first.npoint[k].irrep.irrep_mom.row;
        key_mom[k-1]   =  in->first.npoint[k].irrep.irrep_mom.mom;
        key_irrep[k-1] =  in->first.npoint[k].irrep.op.ops[1].name;


        // for directory naming
        int t_str = in->first.npoint[k].t_slice;      
        if(t_str < 0){removeSubstrs(name, "tm"+std::to_string(std::abs(t_str))+",");}
        else{removeSubstrs(name, "t"+std::to_string(t_str)+",");}
      
      }

      //name of directory
      dir.insert(make_pair(pn,name));
        
      key_struct::KeyHadronSUNNPartNPtIrrep_t tmp_k(key_row, key_mom, key_irrep);

      if(key.find(tmp_k) == key.end()){key.insert(make_pair(tmp_k,key_count)); pn = key_count; key_count++;}
      else{pn = key.find(tmp_k)->second;}



      int dt = std::abs(in->first.npoint[1].t_slice - in->first.npoint[3].t_slice ); 

      for(int t = 0; t < in->second[0].numElem(); t++){ //numElem gives the number of observables in the ensem so use any index for convinience use 0

        ENSEM::EnsemReal rl; rl.resize(in->second.size()); //resize the ensem to the size of the vector = number of bins 

        for(int bin = 0; bin < in->second.size(); bin++){
          rl.elem(bin) = std::exp(Ef*(dt - t))*std::exp(Ei*t)*SEMBLE::toScalar(real(peekObs(in->second[bin], t)));
        }

        x_t.push_back(make_pair(dt, t));
        y_t_ensem.push_back( rl );


      }

      if(pn >= x.size()){
        x_vec.push_back(x_t); y_vec.push_back(y_t_ensem);
        x.resize(pn+1, x_vec); y_ensem.resize(pn+1, y_vec);
      }

      else{
        x_vec = x[pn]; y_vec = y_ensem[pn];
        x_vec.push_back(x_t); y_vec.push_back(y_t_ensem);
        x[pn] = x_vec; y_ensem[pn] = y_vec;
      }

      
      key_row.clear(); key_mom.clear(); key_irrep.clear(); x_vec.clear(); y_vec.clear();
      x_t.clear();  y_t_ensem.clear();
    
    }

    corr.clear();

  }


  //============================
  //===== MAKE THE DATA OBJECT =====

  /* threading over irreps */
  // int nthr = omp_get_max_threads(); 
  // cout << "*** distributing irreps over " << nthr << " threads ***" << endl;
  // #pragma omp parallel for //num_threads(nthr) //schedule(dynamic)
  
  for(int n = 0; n < x.size();n++) //loop over the irreps in the edb
  {

    int size = x[n].size();
    vector<Data> data(size+1);


    for(int i = 0; i < size; i++ )
    {
      make_pair_int_abscissa_data(x[n][i], y_ensem[n][i], data[i] );
      if(i == 0){x_total = x[n][i]; y_ensem_total = y_ensem[n][i];}
      else
      {
        x_total.insert(std::end(x_total), std::begin(x[n][i]), std::end(x[n][i]));
        y_ensem_total.insert(std::end(y_ensem_total), std::begin(y_ensem[n][i]), std::end(y_ensem[n][i]));        
      }

    }

    make_pair_int_abscissa_data(x_total, y_ensem_total, data[size] );
    
    cout << "loaded data::" << endl << data[size].print_data() << endl;

    //============================

    //=================================
    //===== REMOVE UNWANTED POINTS ====


    vector<vector<bool> > keep(size+1);

    vector<bool> remove_dt_all = !(data[size].make_active_data()); //all false
    vector<bool> remove_noisy_all = data_below_y_noise_ratio( data[size], noise); //removes all noisy data points
    vector<bool> tmin_tmax_all = !(data[size].make_active_data()); //all false;


    for(int i = 0; i < size; i++ )
    {

      vector<bool> remove_dt = !(data[i].make_active_data()); //all false
      {
        Abscissa* x_dt = new PairIntAbscissa(make_pair(dt_v[i],dt_v[i]));
        remove_dt = remove_dt || !( data_at_x(data[i], x_dt ) );
        remove_dt_all = remove_dt_all || !( data_at_x(data[size], x_dt ) );
        delete x_dt;
      }
      
      vector<bool> remove_noisy = data_below_y_noise_ratio( data[i], noise);
      
      vector<bool> tmin_tmax = !(data[i].make_active_data()); //all false;

      {
        Abscissa* x_tmin = new PairIntAbscissa(make_pair(dt_v[i], tmin_v[i] - 1 ));
        Abscissa* x_tmax = new PairIntAbscissa(make_pair(dt_v[i], tsnk_v[i] + 1 ));
        tmin_tmax = tmin_tmax || data_in_x_range( data[i], make_pair(x_tmin, x_tmax) );
        tmin_tmax_all = tmin_tmax_all || data_in_x_range( data[size], make_pair(x_tmin, x_tmax) );
        delete x_tmin; delete x_tmax;
      }

      keep[i] = ( remove_dt && remove_noisy && tmin_tmax );
      remove_dt.clear(); remove_noisy.clear(); tmin_tmax.clear();

    }

    keep[size] = ( remove_dt_all && remove_noisy_all && tmin_tmax_all );
    
    cout << "acceptable data::" << endl << data[size].print_data(keep[size]) << endl;

    cout << "###############################" << endl;


    remove_dt_all.clear(); remove_noisy_all.clear(); tmin_tmax_all.clear();

    //===================================


    //============================
    //======= CONTROL FILE  ======== 

    
    vector<fit_three_point_control> control(size+1);
    
    for(int i = 0; i < size; i++)
    {

      control[i].tsrc.push_back(tsrc_v[i]);
      control[i].tsnk.push_back(tsnk_v[i]);
      control[i].dt.push_back(dt_v[i]); 
      control[i].tmin_max.push_back(tmin_max_v[i]); 

      control[size].tsrc.push_back(tsrc_v[i]);
      control[size].tsnk.push_back(tsnk_v[i]);
      control[size].dt.push_back(dt_v[i]); 
      control[size].tmin_max.push_back(tmin_max_v[i]);   


      control[i].Nt_min.push_back(Nt_min_v[i]);
      control[size].Nt_min.push_back(Nt_min_v[i]);
      control[i].plots = true ;
    }
    
    control[size].plots = true ;

    //===================================

    for(int i = 0; i < size; i++ )
    {
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

    //============================
    //======= DO THE FITS ======== 

    fit_three_point_output output =  fit( data, keep, control, fit_qual, chisq_cut, size);

    //============================

    
    //=============================
    //======= OUTPUT STUFF ========  

    std::string name = dir.find(n)->second;
    std::stringstream ss;
    ss << name;
    std::string path = SEMBLE::SEMBLEIO::getPath() += ss.str();
    SEMBLE::SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");
    
    /* write log to file */
    {
      stringstream s; s << path << xmlini << "_three_pt_fit.log"; 
      ofstream out; out.open(s.str().c_str());
      out << output.fit_summary;
      out.close();
    }
    
    /* write mass ensem file */
    {  ostringstream outfile; outfile << path <<  xmlini << "_three_pt_fit.jack";  write(outfile.str(), output.F );  }

    /* write mass variations to file */
    {
      stringstream s; s << path <<  xmlini << "_three_pt_fit.syst"; 
      ofstream out; out.open(s.str().c_str());
      out << output.F_fit_variation;
      out.close();
    }

    // /* write plot data to file */
    {
      stringstream s; s << path << xmlini << "_three_pt_fit.plot"; 
      ofstream out; out.open(s.str().c_str());
      out << output.plot_data;
      out.close();
    }
  


    data.clear(); keep.clear(); control.clear(); 
    delete fit_qual;
  } // end of loop over the irreps in the edb
};
