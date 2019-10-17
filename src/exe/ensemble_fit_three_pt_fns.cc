/*
  read a corr ensem file and do a set of fits to it for various timeslice regions with constant, one and two exponentials and find the best fit
*/

// three_pt_fit
#include "lib/pair_int_abscissa_data_builder.h"
#include "lib/three_point_timeslice_fitting.h"
#include "lib/edb_reader.h"
#include "lib/key_struct.h"
#include "lib/rotations.h"

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

  map<string, vector< pair< pair<ENSEM::EnsemReal, ENSEM::EnsemReal> , ENSEM::EnsemReal >>> fqsq; 
  map<string, vector<ENSEM::EnsemReal> > favg, fmean;
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
  double CurrFac;
  double Zt, Zs;
  string Zt_file, Zs_file; 
  string Zfile;

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
      read(xml_in,"/ThreeptIniParams/Renormalization/Zt",Zt_file);
      read(xml_in,"/ThreeptIniParams/Renormalization/Zs",Zs_file);
      read(xml_in,"/ThreeptIniParams/Renormalization/InvCurrFac2",CurrFac);
      read(xml_in,"/ThreeptIniParams/Renormalization/writeZFile",Zfile);


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
  std::map<int, string> dir;
  std::map<int, int> size_of_base;
  std::map<int, prefactor> pref;
  int key_count = 0;
  int corr_num;
  double mom_sq_i, mom_sq_f;


  /*
  This is the object that is used to make the data object that goes into the fitting routine. The edbs are read and each corr has a Dt and a base corr
  it is a rotation of. We average over all base corrs and the final object will contain the base corrs and the dts
  */  

  vector<vector<vector<pair<int,int>>> >x; vector<vector< vector<ENSEM::EnsemReal> >> y_ensem_rl; vector<vector< vector<ENSEM::EnsemReal> >> y_ensem_im;
  vector<pair<int,int>> x_t; vector<ENSEM::EnsemReal> ytensem_rl; vector<ENSEM::EnsemReal> ytensem_im;  
  vector<vector<pair<int,int>>> x_vec; vector< std::vector<ENSEM::EnsemReal> > y_vec_rl; vector< std::vector<ENSEM::EnsemReal> > y_vec_im;  

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
      std::string base_name, name;
      prefactor pf_tmp;

      int np = corr->first.npoint.size();

      vector<XMLArray::Array<int>> mom(np); vector<int> key_row(np); vector<string> key_irrep(np);
      vector<XMLArray::Array<int>> key_mom(np);
      string Ei_name,Ef_name;

      for(int k = 1; k <= np; k++ ) //Array1DO has 1-based indices
      {
        key_row[k-1]   =  corr->first.npoint[k].irrep.irrep_mom.row;
        mom[k-1]   =  corr->first.npoint[k].irrep.irrep_mom.mom;

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

        mom_sq_i = pow(mom[0][0],2) + pow(mom[0][1],2) + pow(mom[0][2],2);
        mom_sq_f = pow(mom[2][0],2) + pow(mom[2][1],2) + pow(mom[2][2],2); 

        key_mom[0] = Hadron::canonicalOrder(mom[0]);

        Hadron::CubicCanonicalRotation_t ref_angles_z, ref_angles; 


        if(mom_sq_i == 0.0){
          ref_angles_z.alpha = ref_angles_z.beta = ref_angles_z.gamma = 0.0;
          ref_angles.alpha = ref_angles.beta = ref_angles.gamma = 0.0;
        }

        else{
          ref_angles_z = Hadron::cubicCanonicalRotation(mom[0]);
          ref_angles = Hadron::cubicCanonicalRotation(key_mom[0]); 
        }

        key_mom[1] = Rot::EulerRotVec_t(ref_angles_z, ref_angles, mom[1]);
        key_mom[2] = Rot::EulerRotVec_t(ref_angles_z, ref_angles, mom[2]);


        //name of directory

        base_name = key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(key_mom[0][0]) + to_string(key_mom[0][1]) + to_string(key_mom[0][2]) + "__" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(key_mom[1][0]) + to_string(key_mom[1][1]) + to_string(key_mom[1][2]) + "__" +
              key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(key_mom[2][0]) + to_string(key_mom[2][1]) + to_string(key_mom[2][2]);  


        name = key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(mom[0][0]) + to_string(mom[0][1]) + to_string(mom[0][2]) + "__" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(mom[1][0]) + to_string(mom[1][1]) + to_string(mom[1][2]) + "__" +
              key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(mom[2][0]) + to_string(mom[2][1]) + to_string(mom[2][2]); 

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

        key_mom[2] = Hadron::canonicalOrder(mom[2]); 
        
        mom_sq_f = pow(mom[0][0],2) + pow(mom[0][1],2) + pow(mom[0][2],2);
        mom_sq_i = pow(mom[2][0],2) + pow(mom[2][1],2) + pow(mom[2][2],2); 


        Hadron::CubicCanonicalRotation_t ref_angles_z, ref_angles; 


        if(mom_sq_i == 0.0){
          ref_angles_z.alpha = ref_angles_z.beta = ref_angles_z.gamma = 0.0;
          ref_angles.alpha = ref_angles.beta = ref_angles.gamma = 0.0;
        }

        else{
          ref_angles_z = Hadron::cubicCanonicalRotation(mom[2]);
          ref_angles = Hadron::cubicCanonicalRotation(key_mom[2]); 
        }

        key_mom[0] = Rot::EulerRotVec_t(ref_angles_z, ref_angles, mom[0]);
        key_mom[1] = Rot::EulerRotVec_t(ref_angles_z, ref_angles, mom[1]);


        //name of directory    

        base_name = key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(key_mom[2][0]) + to_string(key_mom[2][1]) + to_string(key_mom[2][2]) + "__" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(key_mom[1][0]) + to_string(key_mom[1][1]) + to_string(key_mom[1][2]) + "__" +
              key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(key_mom[0][0]) + to_string(key_mom[0][1]) + to_string(key_mom[0][2]); 


        name = key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(mom[2][0]) + to_string(mom[2][1]) + to_string(mom[2][2]) + "__" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(mom[1][0]) + to_string(mom[1][1]) + to_string(mom[1][2]) + "__" +
              key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(mom[0][0]) + to_string(mom[0][1]) + to_string(mom[0][2]); 

      }

      EnsemReal Ei_ensem; read(Ei_name, Ei_ensem); EnsemReal Ef_ensem; read(Ef_name, Ef_ensem); EnsemReal qsq; qsq.resize(Ei_ensem.size());
      Ei_ensem = rescaleEnsemDown( Ei_ensem );
      Ef_ensem = rescaleEnsemDown( Ef_ensem );
      double q;

      /* Find the Q^2 */
      for(int bin = 0; bin < Ei_ensem.size(); bin++){

        double Ef = SEMBLE::toScalar(Ef_ensem.elem(bin));
        double Ei = SEMBLE::toScalar(Ei_ensem.elem(bin));
        q =  pow(mom[1][0],2) + pow(mom[1][1],2) + pow(mom[1][2],2);
        qsq.elem(bin) = (pow(2*PI/(L*Xi), 2) * q) - pow(Ef - Ei, 2);

      }

      /* Make an ensemble of the two-point function energies */
      std::stringstream path, kfpath;
      path << "Q2_";
      path << std::fixed << std::setprecision(6) << toDouble( mean(qsq) ); 
      path << "/" + base_name + "/" + name + "/" ;
      kfpath << path.str() + kf;
      EnsemComplex kfac_ensem; read(kfpath.str(), kfac_ensem);
      kfac_ensem = rescaleEnsemDown(kfac_ensem);

      if(( toDouble( mean(real(kfac_ensem)) ) == 0.0) && (toDouble( mean(imag(kfac_ensem)) ) == 0.0) ){continue;} // rule out corrs with kfac = 0.0


      if(corr->second.size() != Ei_ensem.size()){ 
        cout << "Ei size and the edb ensem size don't match" << endl;
        exit(1);
      }
      
      if(corr->second.size() != Ef_ensem.size()){ 
        cout << "Ef size and the edb ensem size don't match" << endl;
        exit(1);
      }
        
      key_struct::KeyHadronSUNNPartNPtIrrep_t tmp_k(key_row, key_mom, key_irrep);

      if(key.find(tmp_k) == key.end()){
        key.insert(make_pair(tmp_k,key_count)); corr_num = key_count; dir.insert(make_pair(corr_num,base_name)); key_count++;}

      else{corr_num = key.find(tmp_k)->second;}


      if(size_of_base.find(corr_num) == size_of_base.end()){size_of_base.insert(make_pair(corr_num,1));}
      else{size_of_base.find(corr_num)->second += 1;}


      int dt = std::abs(corr->first.npoint[1].t_slice - corr->first.npoint[3].t_slice ); 

      EnsemVectorComplex yt_post; yt_post.resize(corr->second.size()); yt_post.resizeObs(corr->second[0].numElem()); //to write the files 
      EnsemVectorComplex yt_pre; yt_pre.resize(corr->second.size()); yt_pre.resizeObs(corr->second[0].numElem()); //to write the files  


      for(int t = 0; t < corr->second[0].numElem(); t++){ //numElem gives the number of observables in the ensem so use any index for convinience use 0

        ENSEM::EnsemReal yt_rl; yt_rl.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 
        ENSEM::EnsemReal yt_im; yt_im.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 
        ENSEM::EnsemComplex yt; yt.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 

        for(int bin = 0; bin < corr->second.size(); bin++){//corr num issue

          yt_rl.elem(bin) = SEMBLE::toScalar(real(peekObs(corr->second[bin], t))); 
          yt_im.elem(bin) = SEMBLE::toScalar(imag(peekObs(corr->second[bin], t))); 

          ENSEM::RComplex<double> tmp(SEMBLE::toScalar(yt_rl.elem(bin)),SEMBLE::toScalar(yt_im.elem(bin))); //for some reason cannot read in a complex corr have to stick to this method
          yt.elem(bin)    = tmp;

        }

        pokeObs(yt_pre, yt, t);

        yt_rl = rescaleEnsemDown(yt_rl);
        yt_im = rescaleEnsemDown(yt_im);


        for(int bin = 0; bin < corr->second.size(); bin++){//corr num issue


          double Ef, Ei;

          Ef = SEMBLE::toScalar(Ef_ensem.elem(bin)); 
          Ei = SEMBLE::toScalar(Ei_ensem.elem(bin));

          if(!divkfac){
            yt_rl.elem(bin) = std::exp(Ef*(dt - t))/(sqrt(CurrFac)) * std::exp(Ei*t) * SEMBLE::toScalar(yt_rl.elem(bin)); //multiplying by the exp
            yt_im.elem(bin) = std::exp(Ef*(dt - t))/(sqrt(CurrFac)) * std::exp(Ei*t) * SEMBLE::toScalar(yt_im.elem(bin)); 
                                                 }
          else{

	          complex<double> ff(SEMBLE::toScalar(yt_rl.elem(bin)),SEMBLE::toScalar(yt_im.elem(bin))) ;
            complex<double> kf = SEMBLE::toScalar(kfac_ensem.elem(bin));

            yt_rl.elem(bin) = std::exp(Ef*(dt - t))/(sqrt(CurrFac)) * std::exp(Ei*t) * real(ff/kf); //multiplying by the exp

            yt_im.elem(bin) = std::exp(Ef*(dt - t))/(sqrt(CurrFac)) * std::exp(Ei*t) * imag(ff/kf); //multiplying by exp  
          }

          ENSEM::RComplex<double> tmp(SEMBLE::toScalar(yt_rl.elem(bin)),SEMBLE::toScalar(yt_im.elem(bin))); //for some reason cannot read in a complex corr have to stick to this method
          yt.elem(bin)    = tmp;

        }

        yt_rl = rescaleEnsemUp(yt_rl);
        yt_im = rescaleEnsemUp(yt_im);
        yt    = rescaleEnsemUp(yt);

        pokeObs(yt_post, yt, t);

        x_t.push_back(make_pair(dt, t));
        ytensem_rl.push_back( yt_rl ); ytensem_im.push_back(yt_im);

      }

      {
        ostringstream outfile; outfile << path.str() << "/corr.jack";
        write(outfile.str(), yt_pre);
      }   
      {
        ostringstream outfile; outfile << path.str() << "/normalized_corr.jack";
        write(outfile.str(), yt_post);
      }

      pf_tmp.ei   = Ei_ensem; //if the vector is the creation operator
      pf_tmp.ef   = Ef_ensem;
      pf_tmp.qsq  = qsq;
      if(pref.find(corr_num) == pref.end()){pref.insert(make_pair(corr_num, pf_tmp));}


      if(corr_num >= x.size()){
        x_vec.push_back(x_t); y_vec_rl.push_back(ytensem_rl); y_vec_im.push_back(ytensem_im);
        x.resize(corr_num+1, x_vec); y_ensem_rl.resize(corr_num+1, y_vec_rl); y_ensem_im.resize(corr_num+1, y_vec_im);
      }

      else{
        x_vec = x.at(corr_num); y_vec_rl = y_ensem_rl.at(corr_num); y_vec_im = y_ensem_im.at(corr_num);
        for(int i = 0; i < x_vec.size(); i++){
          if( dt == (x_vec.at(i)).at(0).first ){
            double sizeo = (size_of_base.find(corr_num)->second - 1);
            double sizen = size_of_base.find(corr_num)->second;

            for(int tin = 0; tin < ytensem_rl.size(); tin++){
              y_vec_rl.at(i).at(tin) =  ((SEMBLE::toScalar(sizeo) * rescaleEnsemDown(y_vec_rl.at(i).at(tin)) ) + rescaleEnsemDown(ytensem_rl.at(tin)))/SEMBLE::toScalar(sizen);
              y_vec_im.at(i).at(tin) =  ((SEMBLE::toScalar(sizeo) * rescaleEnsemDown(y_vec_im.at(i).at(tin)) ) + rescaleEnsemDown(ytensem_im.at(tin)))/SEMBLE::toScalar(sizen);

              y_vec_rl.at(i).at(tin) = rescaleEnsemUp(y_vec_rl.at(i).at(tin));
              y_vec_im.at(i).at(tin) = rescaleEnsemUp(y_vec_im.at(i).at(tin));

            }

            break;
          }

          else{continue;}
        }

        y_ensem_rl.at(corr_num) = y_vec_rl; y_ensem_im.at(corr_num) = y_vec_im; 
      }

      
      key_row.clear(); mom.clear(); key_irrep.clear(); x_vec.clear(); y_vec_rl.clear(); y_vec_im.clear();
      x_t.clear();  ytensem_rl.clear(); ytensem_im.clear();


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
    vector<Data> data_rl(num_dt+1); vector<Data> data_im(num_dt+1);
    vector<pair<int,int> >x_total; vector<ENSEM::EnsemReal> y_ensem_total_rl; vector<ENSEM::EnsemReal> y_ensem_total_im;

    for(int i = 0; i < num_dt; i++ )
    {
      make_pair_int_abscissa_data(x[num_corr][i], y_ensem_rl[num_corr][i], data_rl[i] );
      make_pair_int_abscissa_data(x[num_corr][i], y_ensem_im[num_corr][i], data_im[i] );

      if(i == 0){x_total = x[num_corr][i]; 
        y_ensem_total_rl = y_ensem_rl[num_corr][i];
        y_ensem_total_im = y_ensem_im[num_corr][i];}
      else
      {
        x_total.insert(std::end(x_total), std::begin(x[num_corr][i]), std::end(x[num_corr][i]));
        y_ensem_total_rl.insert(std::end(y_ensem_total_rl), std::begin(y_ensem_rl[num_corr][i]), std::end(y_ensem_rl[num_corr][i]));
        y_ensem_total_im.insert(std::end(y_ensem_total_im), std::begin(y_ensem_im[num_corr][i]), std::end(y_ensem_im[num_corr][i]));        
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
    vector<bool> remove_noisy_all_rl = data_below_y_error(data_rl[num_dt], noise); //removes all noisy data points depending on the absolute error in y
    vector<bool> tmin_tmax_all_rl = !(data_rl[num_dt].make_active_data()); //all false;

    vector<bool> remove_dt_all_im = !(data_im[num_dt].make_active_data()); //all false
    vector<bool> remove_noisy_all_im = data_below_y_error(data_im[num_dt], noise); //removes all noisy data points depending on the absolute error in y
    vector<bool> tmin_tmax_all_im = !(data_im[num_dt].make_active_data()); //all false;
  
    for(int i = 0; i < num_dt; i++ )
    {

      vector<bool> remove_dt_rl = !(data_rl[i].make_active_data()); //remove the contact terms
      vector<bool> remove_dt_im = !(data_im[i].make_active_data()); //remove the contact terms
      {
        Abscissa* x_dt = new PairIntAbscissa(make_pair(dt_v[i],dt_v[i]));
        Abscissa* x_0 = new PairIntAbscissa(make_pair(dt_v[i],0));

        vector<bool> remove_ends_rl = ( data_at_x(data_rl[i], x_dt ) ) || ( data_at_x(data_rl[i], x_0 ) ); 
        vector<bool> remove_ends_im = ( data_at_x(data_im[i], x_dt ) ) || ( data_at_x(data_im[i], x_0 ) ); 

        remove_dt_rl = remove_dt_rl || !(remove_ends_rl); remove_dt_im = remove_dt_im || !(remove_ends_im);

        vector<bool> remove_ends_all_rl = ( data_at_x(data_rl[num_dt], x_dt ) ) || ( data_at_x(data_rl[num_dt], x_0 ) ); 
        remove_dt_all_rl = remove_dt_all_rl || !(remove_ends_all_rl);

        vector<bool> remove_ends_all_im = ( data_at_x(data_im[num_dt], x_dt ) ) || ( data_at_x(data_im[num_dt], x_0 ) ); 
        remove_dt_all_im = remove_dt_all_im || !(remove_ends_all_im);

        delete x_dt, x_0;
      }
      
      vector<bool> remove_noisy_rl = data_below_y_error( data_rl[i], noise); //removes all noisy data points depending on the absolute error in y
      vector<bool> remove_noisy_im = data_below_y_error( data_im[i], noise); //removes all noisy data points depending on the absolute error in y     

      vector<bool> tmin_tmax_rl = !(data_rl[i].make_active_data()); //remove all the points below the tmin and above tmax provided by the user;
      vector<bool> tmin_tmax_im = !(data_im[i].make_active_data()); //remove all the points below the tmin and above tmax provided by the user;

      {
        Abscissa* x_tmin = new PairIntAbscissa(make_pair(dt_v[i], tmin_v[i] - 1 ));
        Abscissa* x_tmax = new PairIntAbscissa(make_pair(dt_v[i], tsnk_v[i] + 1 ));

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
      if( count_active(keep_rl[i]) < Nt_min_v[i] )
        { cerr << "fewer than " << Nt_min_v[i] << " timeslices survive your restrictions, no fits on the real part will be acceptable, exiting ..." << endl; exit(1); }
    }

    for(int i = 0; i < num_dt; i++ )
    {
      if( count_active(keep_im[i]) < Nt_min_v[i] )
        { cerr << "fewer than " << Nt_min_v[i] << " timeslices survive your restrictions, no fits on the imag part will be acceptable, exiting ..." << endl; exit(1); }
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
    ss << "Q2_";
    ss << std::fixed << std::setprecision(6) << toDouble( mean(pref.find(num_corr)->second.qsq) ); 
    ss << "/" + name;
    std::string path = SEMBLE::SEMBLEIO::getPath() += ss.str();
    SEMBLE::SEMBLEIO::makeDirectoryPath(path);
    path += std::string("/");

    cout << path << endl;

    //============================
    //======= DO THE FITS ======== 

    fit_three_point_output output_rl =  fit( data_rl, keep_rl, control, fit_qual, chisq_cut, num_dt);
    fit_three_point_output output_im =  fit( data_im, keep_im, control, fit_qual, chisq_cut, num_dt);

    //=============================
    //======= OUTPUT STUFF ======== 

    prefactor pf = pref.find(num_corr)->second;

    if(output_rl.success){


      // /* write the Q^2 vs F(Q^2) output */
      {

        ENSEM::EnsemReal fq2 = output_rl.F; 

        {
          ostringstream outfile; outfile << path << "Real/";
          SEMBLE::SEMBLEIO::makeDirectoryPath(outfile.str());

          
          outfile << name << "_FReal.jack"; 
          write(outfile.str(), fq2);
        }

      }

      /* write log to file */
      {
        stringstream s; s << path << "Real/" << name << "_three_pt_fit_real.log"; 
        ofstream out; out.open(s.str().c_str());
        out << output_rl.fit_summary;
        out.close();
      }
      
      /* write mass ensem file */
      {  ostringstream outfile; outfile << path << "Real/" <<  name << "_three_pt_fit_real.jack";  write(outfile.str(), output_rl.F );  }

      /* write mass variations to file */
      {
        stringstream s; s << path << "Real/" <<  name << "_three_pt_fit_real.syst"; 
        ofstream out; out.open(s.str().c_str());
        out << output_rl.F_fit_variation;
        out.close();
      }

      // /* write plot data to file */
      {
        stringstream s; s << path << "Real/" <<  name << "_three_pt_fit_real.plot"; 
        ofstream out; out.open(s.str().c_str());
        out << "## name= " << name << endl;
        out << output_rl.plot_data;
        out.close();
      }

    }


    if(output_im.success){

      // /* write the Q^2 vs F(Q^2) output */
      {

        ENSEM::EnsemReal fq2 = output_im.F;//*Zv*Zs*CurrFac;

        {
          ostringstream outfile; outfile << path << "Imag/";
          SEMBLE::SEMBLEIO::makeDirectoryPath(outfile.str());

          outfile << name << "_FImag.jack"; 
          write(outfile.str(), fq2);
        }
        
      }

      /* write log to file */
      {
        stringstream s; s << path << "Imag/" << name << "_three_pt_fit_imag.log"; 
        ofstream out; out.open(s.str().c_str());
        out << output_im.fit_summary;
        out.close();
      }
      
      /* write mass ensem file */
      {  ostringstream outfile; outfile << path << "Imag/" <<  name << "_three_pt_fit_imag.jack";  write(outfile.str(), output_im.F );  }

      /* write mass variations to file */
      {
        stringstream s; s << path << "Imag/" << name << "_three_pt_fit_imag.syst"; 
        ofstream out; out.open(s.str().c_str());
        out << output_im.F_fit_variation;
        out.close();
      }

      // /* write plot data to file */
      {
        stringstream s; s << path << "Imag/" << name << "_three_pt_fit_imag.plot"; 
        ofstream out; out.open(s.str().c_str());
        out << "## name= " << name << endl;
        out << output_im.plot_data;
        out.close();
      }

    }

    //=============================================================================
    //======= DECIDE IF THE IMAGINARY PART IS SIGNIFICANT AND OUTPUT STUFF ======== 

    ENSEM::EnsemReal fq2;


    bool out_rl = false;
    bool out_im = false;
    
    if(output_rl.success){
      pair<double,double>  const_rl = mean_err(output_rl.F);
      if( (abs(const_rl.first) - (3.*(const_rl.second)) ) > 0. ){out_rl = true;}
    }
    if(output_im.success){
      pair<double,double>  const_im = mean_err(output_im.F);
      if( (abs(const_im.first) - (3.*(const_im.second)) ) > 0. ){out_im = true;}
    }
    

    if(output_rl.success && output_im.success){
      EnsemReal ratio_ensem = rescaleEnsemUp(  atan( rescaleEnsemDown(output_im.F)/rescaleEnsemDown(output_rl.F) ));
      pair<double,double> ratio = mean_err(ratio_ensem);     

      if( (0.25 > (abs(ratio.first) + ratio.second)) && ( (abs(ratio.first) - ratio.second) > -0.25) ){

        ostringstream outfile; outfile << path << "fit.log";
        ofstream file; // out file stream
        file.open(outfile.str());
        file << "The Arg(F) is: " << ratio.first << " ( "  << ratio.second << " ) " << endl;
        file << "the real part is: " << mean_err(output_rl.F).first << " ( " << mean_err(output_rl.F).second << " ) " << endl;
        file << "the imag part is: " << mean_err(output_im.F).first << " ( " << mean_err(output_im.F).second << " ) " << endl;
        file << "****************************************************************" << endl << endl;
        if(out_rl){
          file << "The value is real with the value of F: " << mean_err(output_rl.F).first << " ( " << mean_err(output_rl.F).second << " ) " << endl;

          file.close();
              
          fq2 = output_rl.F;

          {
            ostringstream outfile; outfile << path << name << "_F.jack"; 
            write(outfile.str(), fq2);
          }

          if(Zfile == "true"?true:false)
          {

            ENSEM::EnsemReal Zv = SEMBLE::toScalar(1.0)/output_rl.F; 

            {
              ostringstream outfile; outfile << path << name << "_Zv.jack"; 
              write(outfile.str(), Zv);

            }        

          }          
          
        }

        else{file << "Decided the value is statistically consistent with zero" << endl;

          file.close();        
                
        }

        
      }    


      else{

        ostringstream outfile; outfile << path << "fit.log";
        ofstream file; // out file stream
        file.open(outfile.str());
        file << "The Arg(F) is: " << ratio.first << " ( "  << ratio.second << " ) " << endl;
        file << "the real part is: " << mean_err(output_rl.F).first << " ( " << mean_err(output_rl.F).second << " ) " << endl;
        file << "the imag part is: " << mean_err(output_im.F).first << " ( " << mean_err(output_im.F).second << " ) " << endl;
        file << "****************************************************************" << endl << endl;

        if(out_im){
          file << "The value is imag with the value of F: " << mean_err(output_im.F).first << " ( " << mean_err(output_im.F).second << " ) " << endl;

          file.close();
              
          fq2 = output_im.F;

          {
            ostringstream outfile; outfile << path << name << "_F.jack"; 
            write(outfile.str(), fq2);
          }

          if(Zfile == "true"?true:false)
          {

            ENSEM::EnsemReal Zv = SEMBLE::toScalar(1.0)/output_im.F; 

            {
              ostringstream outfile; outfile << path << name << "_Zv.jack"; 
              write(outfile.str(), Zv);

            }        

          }
        }

        else{file << "Decided the value is statistically consistent with zero" << endl;

          file.close();     

        }

      }

    }

    else if(output_rl.success){
      ostringstream outfile; outfile << path << "fit.log";
      ofstream file; // out file stream
      file.open(outfile.str());

      
      file << "The fits to the imaginary part failed" << endl;
      file << "The value of the real part is: " << mean_err(output_rl.F).first << " ( " << mean_err(output_rl.F).second << " ) " << endl;
      file << "****************************************************************" << endl << endl;
      if(out_rl){
        file << "The value is real with the value of F: " << mean_err(output_rl.F).first << " ( " << mean_err(output_rl.F).second << " ) " << endl; 

        file.close();
          
        fq2 = output_rl.F;

        {
          ostringstream outfile; outfile << path << name << "_F.jack"; 
          write(outfile.str(), fq2);
        }

        if(Zfile == "true"?true:false)
        {

          ENSEM::EnsemReal Zv = SEMBLE::toScalar(1.0)/output_rl.F; 

          {
            ostringstream outfile; outfile << path << name << "_Zv.jack"; 
            write(outfile.str(), Zv);

          }        

        }
      }
      else{file << "Decided the value is statistically consistent with zero" << endl;
          file.close();     
      }

      cout << "All fits to the imaginary part failed " << endl;
    }

    else if(output_im.success){
      ostringstream outfile; outfile << path << "fit.log";
      ofstream file; // out file stream
      file.open(outfile.str());

      file << "The fits to the real part failed" << endl;
      file << "The value of the imag part is: " << mean_err(output_im.F).first << " ( " << mean_err(output_im.F).second << " ) " << endl;
      file << "****************************************************************" << endl << endl;
      if(out_im){
        file << "The value is imag with the value of F: " <<  mean_err(output_im.F).first << " ( " << mean_err(output_im.F).second << " ) " << endl;
        file.close();
            
        fq2 = output_im.F;

        {
          ostringstream outfile; outfile << path << name << "_F.jack"; 
          write(outfile.str(), fq2);
        }

        if(Zfile == "true"?true:false)
        {

          ENSEM::EnsemReal Zv = SEMBLE::toScalar(1.0)/output_im.F; 

          {
            ostringstream outfile; outfile << path << name << "_Zv.jack"; 
            write(outfile.str(), Zv);

          }        

        }
        
        cout << "All fits to the real part failed " << endl;
      }

      else{file << "Decided the value is statistically consistent with zero" << endl;
          file.close();     
      }
    }

    else{ cout << "All fits to the real part and the imaginary part failed" << endl; }



    //=============================
    //======= OUTPUT STUFF ======== 

    if((output_im.success || output_rl.success) && (out_rl || out_im)){

      SEMBLE::SEMBLEIO::makeDirectoryPath(path + "../jackfiles/"); 

      {
        ostringstream outfile; outfile << path << name << "_Q2.jack"; 
        write(outfile.str(), pf.qsq);
      }

      {
        ostringstream outfile; outfile << path << "../jackfiles/" << "Q2.jack"; 
        write(outfile.str(), pf.qsq);
      }

      {
        ostringstream outfile; outfile << path << name << "_Estar.jack"; 
        write(outfile.str(), pf.ei);
      }

      if(fqsq.find(name) == fqsq.end() ){

        vector<pair<pair< ENSEM::EnsemReal, ENSEM::EnsemReal>, ENSEM::EnsemReal>> fqsq_val;
        fqsq_val.push_back(make_pair(make_pair(pf.qsq, pf.ei),fq2));
        fqsq.insert(make_pair(name, fqsq_val)); 

        }

      else{ fqsq.find(name)->second.push_back(make_pair(make_pair(pf.qsq, pf.ei),fq2)); }

      stringstream s; s <<  std::fixed << std::setprecision(6) << toDouble(mean(pf.qsq));
      string q2 = s.str();

      stringstream ss; ss <<  std::fixed << std::setprecision(6) << toDouble(mean(pf.ei));
      string ei = s.str();

      if(favg.find(q2) == favg.end()){ vector<ENSEM::EnsemReal> val; val.push_back(fq2); favg.insert( make_pair(q2,val) );}
      else{favg.find(q2)->second.push_back(fq2);}

      if(fmean.find(ei) == fmean.end()){ vector<ENSEM::EnsemReal> val; val.push_back(fq2); fmean.insert( make_pair(ei,val) );}
      else{fmean.find(ei)->second.push_back(fq2);}
    }

    data_rl.clear(); keep_rl.clear(); control.clear(); data_im.clear(); keep_im.clear(); 
    delete fit_qual;

    cout << "finisihed fitting " << num_corr + 1 << " corrs out of " <<  x.size() << " corrs" << endl;
    cout << "====================================================================================" << endl;

  } // END OF THE LOOP OVER THE CORRS IN THE EBD

  cout << endl << "*******************************************************************************************************************************************" << endl;
  cout << " Finished fiiting all the three point correlation functions in the given edbs. " << endl << endl;



  /* Write the F(Q^2) vs Q^2 plot */

  if(divkfac){
    
    std::string path = SEMBLE::SEMBLEIO::getPath();
    std::string outdir = path + "output/";
    SEMBLE::SEMBLEIO::makeDirectoryPath(outdir);


    {
      stringstream s; s << outdir <<  xmlini << "_FQ2.out";         
      ofstream out; out.open(s.str().c_str());
      out << "## Q  Q_E  F  F_E" << endl << endl;
      for(auto it = fqsq.begin(); it != fqsq.end(); it++ )
      {

        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){  
          pair<double, double> q_out = mean_err(it1->first.first);
          pair<double, double> f_out = mean_err(it1->second);
          out << q_out.first << "  " << q_out.second << "  " <<  f_out.first <<  "  " << f_out.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

    {
      stringstream s; s << outdir <<  xmlini << "_ZvQ2.out";         
      ofstream out; out.open(s.str().c_str());
      out << "## Q  Q_E  F  F_E" << endl << endl;
      for(auto it = fqsq.begin(); it != fqsq.end(); it++ )
      {

        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){  
          pair<double, double> q_out = mean_err(it1->first.first);
          pair<double, double> f_out = mean_err(SEMBLE::toScalar(1)/it1->second);
          out << q_out.first << "  " << q_out.second << "  " <<  f_out.first <<  "  " << f_out.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

    {
      stringstream s; s << outdir <<  xmlini << "_FMv.out";         
      ofstream out; out.open(s.str().c_str());
      out << "## E  E_E  F  F_E" << endl << endl;
      for(auto it = fqsq.begin(); it != fqsq.end(); it++ )
      {

        out << "## irrep=" << it->first << endl;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){ 
          pair<double, double> e_out = mean_err(it1->first.second);
          pair<double, double> f_out = mean_err(it1->second);
          out << e_out.first << "  " << e_out.second << "  " <<  f_out.first <<  "  " << f_out.second << endl;}
          out << endl << endl;
      }

      out.close();  
    }

    std::string plotdir = path + "plots/";
    SEMBLE::SEMBLEIO::makeDirectoryPath(plotdir); 



    if(Zfile == "true"?true:false){
      stringstream ff; ff << plotdir <<  xmlini << "_FF.plot";   
      stringstream Zvf; Zvf << plotdir <<  xmlini << "_Zv.plot";   
            
      ofstream write_ff; write_ff.open(ff.str().c_str());
      write_ff << "## Q2        F        F_E" << endl << endl;

      ofstream write_zvf; write_zvf.open(Zvf.str().c_str());
      write_zvf << "## Q2        Zv        Zv_E" << endl << endl;

      for(auto it = favg.begin(); it != favg.end(); it++ ){
        ENSEM::EnsemReal avg = it->second.at(0);
        for(auto it1 = it->second.begin()+1; it1 != it->second.end(); it1++){avg += (SEMBLE::toScalar(1/it->second.size()) * (*it1));}


        stringstream s; s << path << "Q2_" << it->first << "/jackfiles/FF.jack";
        write(s.str(),avg);
        stringstream ss; ss << path << "Q2_" << it->first << "/jackfiles/Zv.jack"; ENSEM::EnsemReal Zavg = SEMBLE::toScalar(1)/avg; write(ss.str(), Zavg);  

        pair<double, double> f_avg =  mean_err(avg); 
        pair<double, double> Z_avg =  mean_err(SEMBLE::toScalar(1)/avg);    

        write_ff << it->first << " " << f_avg.first << " " <<  f_avg.second << endl;
        write_zvf << it->first << " " << Z_avg.first << " " <<  Z_avg.second << endl;

      }

      write_ff.close(); 
      write_zvf.close();
    }


    else{
      stringstream ff; ff << plotdir <<  xmlini << "_FF.plot";   
  
            
      ofstream write_ff; write_ff.open(ff.str().c_str());
      write_ff << "## Q2        F        F_E" << endl << endl;


      for(auto it = favg.begin(); it != favg.end(); it++ ){
        ENSEM::EnsemReal avg = it->second.at(0);

        for(auto it1 = it->second.begin()+1; it1 != it->second.end(); it1++){avg += (SEMBLE::toScalar(1/it->second.size()) * (*it1));}

        stringstream s; s << path << "Q2_" << it->first << "/FF.jack";


        write(s.str(),avg);

        pair<double, double> f_avg =  mean_err(avg);     
        write_ff << it->first << " " << f_avg.first << " " <<  f_avg.second << endl;


      }

      write_ff.close(); 
    }


  }

};
