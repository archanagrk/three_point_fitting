/*
  read a corr ensem file and do a set of fits to it for various timeslice regions with constant, one and two exponentials and find the best fit
*/

// three_pt_fit
#include "lib/pair_int_abscissa_data_builder.h"
#include "lib/three_point_timeslice_fitting.h"
#include "lib/edb_reader.h"
#include "lib/key_struct.h"

// fitting lib
#include "fitting_lib/data_selectors.h"

// semble
#include"semble/semble_file_management.h"
#include "semble/semble_meta.h"

// adat
#include "hadron/ensem_filenames.h"

using namespace key_struct;




int main(int argc, char** argv)
{

  const double PI = (atan(double(1)) * double(4.0));

  map<string, vector< pair< pair<ENSEM::EnsemReal, ENSEM::EnsemReal> , ENSEM::EnsemComplex >>> fqsqmv; 
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
  double Ei_min, Ef_min, Ei_max, Ef_max;

  vector<string> filename;
  string file;

  ADAT::Array1dO<string> Ei_jack, Ef_jack ; 
  double noise, chisq_cut;
  string qual_type;

  vector<int> dt_v, tmin_v, tmax_v, tsrc_v, tsnk_v, Nt_min_v; 
  vector<double> Ei_min_v, Ef_min_v, Ei_max_v, Ef_max_v;
  vector<std::tuple<int,int,int>> tmin_max_v;

  string kfac_xml;
  bool divkfac;

  XMLReader xml_in(xmlini);
  try
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
        read(xml_in,"/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/EiMin",Ei_min);
        read(xml_in,"/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/EfMin",Ef_min);
        read(xml_in,"/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/EiMax",Ei_max);
        read(xml_in,"/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/EfMax",Ef_max);
        read(xml_in,"/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/NtMin",Nt_min);




        std::tuple<int,int,int> tmin_max = std::make_tuple(dt, tsrc,tsnk);

        filename.push_back(file);
        dt_v.push_back(dt); tmin_v.push_back(tmin); tmax_v.push_back(tmax); tsrc_v.push_back(tsrc); tsnk_v.push_back(tsnk);
        tmin_max_v.push_back(tmin_max); Nt_min_v.push_back(Nt_min); Ei_min_v.push_back(Ei_min); Ef_min_v.push_back(Ef_min);
        Ei_max_v.push_back(Ei_max); Ef_max_v.push_back(Ef_max);
      }

      read(xml_in,"/ThreeptIniParams/FitProps/fit_qual/type",qual_type);
      read(xml_in,"/ThreeptIniParams/FitProps/chisq_cutoff",chisq_cut);
      read(xml_in,"/ThreeptIniParams/FitProps/noise_cutoff",noise);
      read(xml_in,"/ThreeptIniParams/FitProps/divKfac",divkfac);
      read(xml_in,"/ThreeptIniParams/FitProps/kfac",kfac_xml);
      read(xml_in,"/ThreeptIniParams/FitProps/Ei",Ei_jack);      
      read(xml_in,"/ThreeptIniParams/FitProps/Ef",Ef_jack);


    }
  catch( const string& error ){
    cerr << "Error reading input file : " << error << endl;
    }


  //=============================
  //===== READ THE KFAC XML =====

    string kf;
    double Xi, XiE, as, asE; 
    string mf_file;
    int L, pts;
    EnsemReal mf; 


    if(divkfac){

    XMLReader xml_kf_in(kfac_xml);
    try
      {
        read(xml_kf_in, "/kfac/pts", pts);


        read(xml_kf_in,"/kfac/L",L);
        read(xml_kf_in,"/kfac/as",as);
        read(xml_kf_in,"/kfac/asE",asE);
        read(xml_kf_in,"/kfac/Xi",Xi);
        read(xml_kf_in,"/kfac/XiE",XiE);

        read(xml_kf_in,"/kfac/mfFile",mf_file);
        read(xml_kf_in,"/kfac/kfacFile",kf);



      }
    catch( const string& error ){
      cerr << "Error reading input file : " << error << endl;
      }

      read(mf_file, mf);


    }



  //============================
  //===== INTERFACE TO EDB =====


  // Read all the correlators from an edb file the edbs need to be the same irreps they are not verified. User has to make sure
  std::map<key_struct::KeyHadronSUNNPartNPtIrrep_t, int> key; 
  std::map<int,string> dir;
  std::map<int, prefactor> pref;
  int key_count = 0;
  int corr_num;
  double mom_sq_i, mom_sq_f;

  vector<vector<vector<pair<int,int>>> >x; vector<vector< vector<ENSEM::EnsemReal> >> y_ensem;
  vector<pair<int,int>> x_t; vector<ENSEM::EnsemReal> y_t_ensem; 
  vector<vector<pair<int,int>>> x_vec; vector< std::vector<ENSEM::EnsemReal> > y_vec;  

  for(int i = 0; i < num; i++ )
  {

    std::cout << "read file = " << filename[i] << std::endl;

    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > corrs = edb_reader::getCorrs(filename[i]);

    std::cout << "end of file =  " << filename[i] << std::endl;


    //=================================================
    //===== LOAD THE DATA IN THE EDB AND THE KFAC =====

      
    for(ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> >::const_iterator  corr = corrs.begin(); corr != corrs.end(); ++corr) 
    {
      
      //name of directory
      std::string name;
      prefactor pf_tmp;

      int np = corr->first.npoint.size();

      vector<XMLArray::Array<int>> key_mom(np); vector<int> key_row(np); vector<string> key_irrep(np);
      string Ei_name,Ef_name;

      for(int k = 1; k <= np; k++ ) //Array1DO has 1-based indices
      {
        key_row[k-1]   =  corr->first.npoint[k].irrep.irrep_mom.row;
        key_mom[k-1]   =  corr->first.npoint[k].irrep.irrep_mom.mom;

        string irp = corr->first.npoint[k].irrep.op.ops[1].name;
        std::size_t found = irp.find_last_of("_");

        key_irrep[k-1] =  irp.substr(found+1);
      
      }

      //==========================================================================================
      //===== DIVIDE BY THE EXP OF SOURCE ENERGY LEVEL * t AND SINK ENERGY LEVEL * (Dt - t)  =====

      /* Find the required MassJackFile for dividing out the exponential of the energies */
      if(corr->first.npoint[1].irrep.creation_op && !corr->first.npoint[3].irrep.creation_op)
      { 
        for(int fn = 1; fn <= Ei_jack.size(); fn++){
          if( (corr->first.npoint[1].irrep.op.ops[1].name+".jack" == Ei_jack[fn]) ){Ei_name = Ei_jack[fn]; break;}
        }

        for(int fn = 1; fn <= Ef_jack.size(); fn++){
          if( (corr->first.npoint[3].irrep.op.ops[1].name+".jack" == Ef_jack[fn]) ){Ef_name = Ef_jack[fn]; break;}
        }
        if(!Ei_name.size()){cerr << "No Ei jackfile provided for the required initial state: " <<
                        corr->first.npoint[1].irrep.op.ops[1].name+".jack" << endl; exit(1);}
        if(!Ef_name.size()){cerr << "No Ef jackfile provided for the required final state: " << 
                        corr->first.npoint[3].irrep.op.ops[1].name+".jack" << endl; exit(1);}

        //name of directory
        name = key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(key_mom[0][0]) + to_string(key_mom[0][1]) + to_string(key_mom[0][2]) + "__" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(key_mom[1][0]) + to_string(key_mom[1][1]) + to_string(key_mom[1][2]) + "__" +
              key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(key_mom[2][0]) + to_string(key_mom[2][1]) + to_string(key_mom[2][2]); 

        mom_sq_i = pow(key_mom[0][0],2) + pow(key_mom[0][1],2) + pow(key_mom[0][2],2);
        mom_sq_f = pow(key_mom[2][0],2) + pow(key_mom[2][1],2) + pow(key_mom[2][2],2); 

      }

      else if(!corr->first.npoint[1].irrep.creation_op && corr->first.npoint[3].irrep.creation_op)
      {

        for(int fn = 1; fn <= Ei_jack.size(); fn++){
          if( (corr->first.npoint[3].irrep.op.ops[1].name+".jack" == Ei_jack[fn]) ){Ei_name = Ei_jack[fn]; break;}
        }

        for(int fn = 1; fn <= Ef_jack.size(); fn++){
          if( (corr->first.npoint[1].irrep.op.ops[1].name+".jack" == Ef_jack[fn]) ){Ef_name = Ef_jack[fn]; break;}
        }
        if(!Ei_name.size()){cerr << "No Ei jackfile provided for the required initial state: " << 
                            corr->first.npoint[3].irrep.op.ops[1].name+".jack" << endl; exit(1);}
        if(!Ef_name.size()){cerr << "No Ef jackfile provided for the required final state: " <<
                          corr->first.npoint[1].irrep.op.ops[1].name+".jack" << endl; exit(1);}

        //name of directory                
        name = key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(key_mom[2][0]) + to_string(key_mom[2][1]) + to_string(key_mom[2][2]) + "__" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(key_mom[1][0]) + to_string(key_mom[1][1]) + to_string(key_mom[1][2]) + "__" +
              key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(key_mom[0][0]) + to_string(key_mom[0][1]) + to_string(key_mom[0][2]); 

        mom_sq_f = pow(key_mom[0][0],2) + pow(key_mom[0][1],2) + pow(key_mom[0][2],2);
        mom_sq_i = pow(key_mom[2][0],2) + pow(key_mom[2][1],2) + pow(key_mom[2][2],2); 

      }


      /* Make an ensemble of the two-point function energies */

      EnsemReal Ei_ensem; read(Ei_name, Ei_ensem); EnsemReal Ef_ensem; read(Ef_name, Ef_ensem); EnsemReal qsq; qsq.resize(Ei_ensem.size());
      double q;

      if(corr->second.size() != Ei_ensem.size()){ 
        cout << "Ei size and the ebd ensem size don't match" << endl;
        exit(1);
      }
      
      if(corr->second.size() != Ef_ensem.size()){ 
        cout << "Ef size and the ebd ensem size don't match" << endl;
        exit(1);
      }
        
      key_struct::KeyHadronSUNNPartNPtIrrep_t tmp_k(key_row, key_mom, key_irrep);

      if(key.find(tmp_k) == key.end()){key.insert(make_pair(tmp_k,key_count)); corr_num = key_count; key_count++;}
      else{corr_num = key.find(tmp_k)->second;}

      //name of directory
      dir.insert(make_pair(corr_num,name));



      int dt = std::abs(corr->first.npoint[1].t_slice - corr->first.npoint[3].t_slice ); 

      for(int t = 0; t < corr->second[0].numElem(); t++){ //numElem gives the number of observables in the ensem so use any index for convinience use 0

        ENSEM::EnsemReal yt; yt.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 


        for(int bin = 0; bin < corr->second.size(); bin++){

          double Ef, Ei;

          Ef = SEMBLE::toScalar(Ef_ensem.elem(bin));
          Ei = SEMBLE::toScalar(Ei_ensem.elem(bin));

          yt.elem(bin) = std::exp(Ef*(dt - t)) * std::exp(Ei*t) * SEMBLE::toScalar(real(peekObs(corr->second[bin], t))); //multiplying by the normalization sqrt(Ei*Ef)*mf

          /* Find the Q^2 */
          if(divkfac && (t == 0)){
            q =  pow(key_mom[1][0],2) + pow(key_mom[1][1],2) + pow(key_mom[1][2],2);
            qsq.elem(bin) = pow(2*PI/(L*Xi), 2) * q - pow(Ef - Ei, 2);
          }
        }

        x_t.push_back(make_pair(dt, t));
        y_t_ensem.push_back( yt );

      }

      string kfpath = name + "/" + kf;

      if(divkfac){
        read(kfpath, pf_tmp.kfac);
        pf_tmp.ei = Ei_ensem; //if the vector is the creation operator
        pf_tmp.ef = Ef_ensem;
        pf_tmp.qsq = qsq;
        if(pref.find(corr_num) == pref.end()){pref.insert(make_pair(corr_num, pf_tmp));}
        }


      if(corr_num >= x.size()){
        x_vec.push_back(x_t); y_vec.push_back(y_t_ensem);
        x.resize(corr_num+1, x_vec); y_ensem.resize(corr_num+1, y_vec);
      }

      else{
        x_vec = x.at(corr_num); y_vec = y_ensem.at(corr_num);
        x_vec.push_back(x_t); y_vec.push_back(y_t_ensem);
        x.at(corr_num) = x_vec; y_ensem.at(corr_num) = y_vec;
      }

      
      key_row.clear(); key_mom.clear(); key_irrep.clear(); x_vec.clear(); y_vec.clear();
      x_t.clear();  y_t_ensem.clear();


    }

    corrs.clear();

  }


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
  for(int num_corr = 0; num_corr < x.size(); num_corr += 1)
  {

    int num_dt = x[num_corr].size();
    vector<Data> data(num_dt+1);
    vector<pair<int,int> >x_total; vector<ENSEM::EnsemReal> y_ensem_total;

    for(int i = 0; i < num_dt; i++ )
    {
      make_pair_int_abscissa_data(x[num_corr][i], y_ensem[num_corr][i], data[i] );

      if(i == 0){x_total = x[num_corr][i]; y_ensem_total = y_ensem[num_corr][i];}
      else
      {
        x_total.insert(std::end(x_total), std::begin(x[num_corr][i]), std::end(x[num_corr][i]));
        y_ensem_total.insert(std::end(y_ensem_total), std::begin(y_ensem[num_corr][i]), std::end(y_ensem[num_corr][i]));        
      }

    }

    make_pair_int_abscissa_data(x_total, y_ensem_total, data[num_dt] );
    
    cout << "loaded data::" << endl << data[num_dt].print_data() << endl;


    //=================================

    //=================================
    //===== REMOVE UNWANTED POINTS ====


    vector<vector<bool> > keep(num_dt+1);

    vector<bool> remove_dt_all = !(data[num_dt].make_active_data()); //all false
    vector<bool> remove_noisy_all = data_below_y_error(data[num_dt], noise); //removes all noisy data points depending on the absolute error in y
    vector<bool> tmin_tmax_all = !(data[num_dt].make_active_data()); //all false;

  
    for(int i = 0; i < num_dt; i++ )
    {

      vector<bool> remove_dt = !(data[i].make_active_data()); //remove the contact terms
      {
        Abscissa* x_dt = new PairIntAbscissa(make_pair(dt_v[i],dt_v[i]));
        Abscissa* x_0 = new PairIntAbscissa(make_pair(dt_v[i],0));
        vector<bool> remove_ends = ( data_at_x(data[i], x_dt ) ) || ( data_at_x(data[i], x_0 ) ); 
        remove_dt = remove_dt || !(remove_ends);

        vector<bool> remove_ends_all = ( data_at_x(data[num_dt], x_dt ) ) || ( data_at_x(data[num_dt], x_0 ) ); 
        remove_dt_all = remove_dt_all || !(remove_ends_all);


        delete x_dt, x_0;
      }
      
      vector<bool> remove_noisy = data_below_y_error( data[i], noise); //removes all noisy data points depending on the absolute error in y
      
      vector<bool> tmin_tmax = !(data[i].make_active_data()); //remove all the points below the tmin and above tmax provided by the user;

      {
        Abscissa* x_tmin = new PairIntAbscissa(make_pair(dt_v[i], tmin_v[i] - 1 ));
        Abscissa* x_tmax = new PairIntAbscissa(make_pair(dt_v[i], tsnk_v[i] + 1 ));
        tmin_tmax = tmin_tmax || data_in_x_range( data[i], make_pair(x_tmin, x_tmax) );
        tmin_tmax_all = tmin_tmax_all || data_in_x_range( data[num_dt], make_pair(x_tmin, x_tmax) );
        delete x_tmin; delete x_tmax;
      }

      keep[i] = ( remove_dt && remove_noisy && tmin_tmax );
      remove_dt.clear(); remove_noisy.clear(); tmin_tmax.clear();

    }

    keep[num_dt] = ( remove_dt_all && remove_noisy_all && tmin_tmax_all );
    
    cout << "acceptable data::" << endl << data[num_dt].print_data(keep[num_dt]) << endl;


    cout << "###############################" << endl;


    remove_dt_all.clear(); remove_noisy_all.clear(); tmin_tmax_all.clear();

    //===================================


    //============================
    //======= CONTROL FILE  ======== 

    
    vector<fit_three_point_control> control(num_dt+1);
    
    for(int i = 0; i < num_dt; i++)
    {

      control[i].tsrc.push_back(tsrc_v[i]);
      control[i].tsnk.push_back(tsnk_v[i]);
      control[i].dt.push_back(dt_v[i]); 
      control[i].tmin_max.push_back(tmin_max_v[i]); 

      control[num_dt].tsrc.push_back(tsrc_v[i]);
      control[num_dt].tsnk.push_back(tsnk_v[i]);
      control[num_dt].dt.push_back(dt_v[i]); 
      control[num_dt].tmin_max.push_back(tmin_max_v[i]);   

      control[i].Nt_min.push_back(Nt_min_v[i]);
      control[i].Ei_min.push_back(Ei_min_v[i]);
      control[i].Ef_min.push_back(Ef_min_v[i]);
      control[i].Ei_max.push_back(Ei_max_v[i]);
      control[i].Ef_max.push_back(Ef_max_v[i]);

      control[num_dt].Nt_min.push_back(Nt_min_v[i]);
      control[num_dt].Ei_min.push_back(Ei_min_v[i]);
      control[num_dt].Ef_min.push_back(Ef_min_v[i]);
      control[num_dt].Ei_max.push_back(Ei_max_v[i]);
      control[num_dt].Ef_max.push_back(Ef_max_v[i]);

      control[i].plots = true ;
    }
    
    control[num_dt].plots = true ;

    //===================================

    for(int i = 0; i < num_dt; i++ )
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
    //======= DIR NAMING ======== 

    /* Creates a subdirectory with the irrep names and the files for the respective fit is stored in this subdirectory */

    std::string name = dir.find(num_corr)->second;
    std::stringstream ss;
    ss << name;
    std::string path = SEMBLE::SEMBLEIO::getPath() += ss.str();
    SEMBLE::SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    cout << path << endl;

    //============================
    //======= DO THE FITS ======== 

    fit_three_point_output output =  fit( data, keep, control, fit_qual, chisq_cut, num_dt);

    //===============================================
    //======= DIVIDE BY THE KINEMATIC FACTOR ======== 

    if(output.success && divkfac){

      
      prefactor pf = pref.find(num_corr)->second;

      // /* write the Q^2 vs F(Q^2) output */
      {

        ENSEM::EnsemComplex fq2 = sqrt(pf.ei* pf.ef) * SEMBLE::toScalar(pow(L,3)/(2.0* Xi)) / pf.kfac * output.F;


        {
          ostringstream outfile; outfile << path << name << "_Q2.jack"; 
          write(outfile.str(), pf.qsq);
        }

        {
          ostringstream outfile; outfile << path << name << "_F.jack"; 
          write(outfile.str(), fq2);
        }

        {
          ostringstream outfile; outfile << path << name << "_Estar.jack"; 
          write(outfile.str(), pf.ei);
        }


        if(fqsqmv.find(name) == fqsqmv.end() ){

          vector<pair<pair< ENSEM::EnsemReal, ENSEM::EnsemReal>, ENSEM::EnsemComplex>> fqsqmv_val;
          fqsqmv_val.push_back(make_pair(make_pair(pf.qsq, pf.ei),fq2));
          fqsqmv.insert(make_pair(name, fqsqmv_val)); 

          }

        else{ fqsqmv.find(name)->second.push_back(make_pair(make_pair(pf.qsq, pf.ei),fq2)); }
        
      }

    }

 
    //=============================
    //======= OUTPUT STUFF ========  


    if(output.success)
    {
      /* write log to file */
      {
        stringstream s; s << path << name << "_three_pt_fit.log"; 
        ofstream out; out.open(s.str().c_str());
        out << output.fit_summary;
        out.close();
      }
      
      /* write mass ensem file */
      {  ostringstream outfile; outfile << path <<  name << "_three_pt_fit.jack";  write(outfile.str(), output.F );  }

      /* write mass variations to file */
      {
        stringstream s; s << path <<  name << "_three_pt_fit.syst"; 
        ofstream out; out.open(s.str().c_str());
        out << output.F_fit_variation;
        out.close();
      }

      // /* write plot data to file */
      {
        stringstream s; s << path << name << "_three_pt_fit.plot"; 
        ofstream out; out.open(s.str().c_str());
        out << "## name= " << name << endl;
        out << output.plot_data;
        out.close();
      }
    }  

    data.clear(); keep.clear(); control.clear(); 
    delete fit_qual;

    cout << "finisihed fitting " << num_corr + 1 << " corrs out of " <<  x.size() << " corrs" << endl;
    cout << "====================================================================================" << endl;

  } // END OF THE LOOP OVER THE CORRS IN THE EBD

  cout << endl << "*******************************************************************************************************************************************" << endl;
  cout << " Finished fiiting all the three point correlation functions in the given edbs. " << endl << endl;

  /* Write the F(Q^2) vs Q^2 plot */

  if(divkfac){
    
    std::string path = SEMBLE::SEMBLEIO::getPath();
    {
      stringstream s; s << path <<  xmlini << "_F_vs_Q2.plot";         
      ofstream out; out.open(s.str().c_str());
      out << "## Q  Q_E  F  F_E" << endl << endl;
      for(auto it = fqsqmv.begin(); it != fqsqmv.end(); it++ )
      {
        //out << it->first.first <<  "  " << it->first.first + it->first.second  << "  " << it->first.first - it->first.second  << "  " << 
              //it->second.first <<  "  " << it->second.first + it->second.second  << "  " << it->second.first - it->second.second  << endl;
        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){  
          pair<double, double> qplot = mean_err(it1->first.first);
          pair<complex<double>, complex<double>> fplot = make_pair( SEMBLE::toScalar( mean(it1->second) ), SEMBLE::toScalar( sqrt( variance( it1->second ) ) ) );
          out << qplot.first << "  " << qplot.second << "  " <<  fplot.first <<  "  " << fplot.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

    {
      stringstream s; s << path <<  xmlini << "_F_vs_mv.plot";         
      ofstream out; out.open(s.str().c_str());
      out << "## E  E_E  F  F_E" << endl << endl;
      for(auto it = fqsqmv.begin(); it != fqsqmv.end(); it++ )
      {
        //out << it->first.first <<  "  " << it->first.first + it->first.second  << "  " << it->first.first - it->first.second  << "  " << 
              //it->second.first <<  "  " << it->second.first + it->second.second  << "  " << it->second.first - it->second.second  << endl;
        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){ 
          pair<double, double> eiplot = mean_err(it1->first.second);
          pair<complex<double>, complex<double>> fplot = make_pair( SEMBLE::toScalar( mean(it1->second) ), SEMBLE::toScalar( sqrt( variance( it1->second ) ) ) );
          out << eiplot.first << "  " << eiplot.second << "  " <<  fplot.first <<  "  " << fplot.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

  }


};
