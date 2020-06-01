
#include "edb_reader.h"

//! Get corrs
ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > getCorrs(const std::string& dbase)
{
  typedef Hadron::KeyHadronSUNNPartNPtCorr_t              K;
  typedef ENSEM::VectorComplex                            V;


  // Open DB
  FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<K>,  ADATIO::SerialDBData<V> > database;

  if (database.open(dbase, O_RDONLY, 0400) != 0)
  {
    std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
    exit(1);
  }

  // Read all the keys and values
  std::vector< ADATIO::SerialDBKey<K> >                  keys;
  std::vector< std::vector< ADATIO::SerialDBData<V> > >  vals;

  database.keysAndData(keys, vals);

  // Put into something more useful
  ADAT::MapObject<K, std::vector<V> > dest;
  
  for(int i = 0; i < keys.size(); ++i) {
    std::vector<V> v(vals[i].size());

    for(int n=0; n < v.size(); ++n) {
        v[n] = vals[i][n].data();
    }

    dest.insert(keys[i].key(), v);
  }

  return dest;

}

 //************************************************************************************************

//! Read the edb data and divide off the leading order exponentials along with the kinematic factors
EDBdata_t read_corrs(XMLinidata_t& xml_data){


  EDBdata_t edb_data;

  const double PI = (atan(double(1)) * double(4.0));

  // Read all the correlators from an edb file the edbs need to be the same irreps they are not verified. User has to make sure
  std::map<ThreePtKeys::KeyHadronSUNNPartNPtIrrep_t, int> key; 
  std::map<int,  vector< vector<EnsemComplex>> > kfacs;

  std::map<int, int> size_of_base;
  int key_count = 0;
  size_t corr_num;
  double mom_sq_i, mom_sq_f;


  /*
  This is the object that is used to make the data object that goes into the fitting routine. The edbs are read and each corr has a Dt and a base corr
  it is a rotation of. We average over all base corrs and the final object will contain the base corrs and the dts
  */  

  vector<pair<int,int>> x_t; vector<ENSEM::EnsemReal> ytensem_rl; vector<ENSEM::EnsemReal> ytensem_im;  
  vector<vector<pair<int,int>>> x_vec; vector< std::vector<ENSEM::EnsemReal> > y_vec_rl; vector< std::vector<ENSEM::EnsemReal> > y_vec_im;  

  for(int i = 0; i < xml_data.num; i++ )
  {

    std::cout << "read file = " << xml_data.filename[i] << std::endl;

    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > corrs = getCorrs(xml_data.filename[i]);

    std::cout << "end of file =  " << xml_data.filename[i] << std::endl;

    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > improv_corrs;

    if(xml_data.improv){
    std::cout << "read improv file = " << xml_data.improvname[i] << std::endl;

    improv_corrs = getCorrs(xml_data.improvname[i]);

    std::cout << "end of improv file =  " << xml_data.improvname[i] << std::endl;      
    }
    
    //=================================================
    //===== LOAD THE DATA IN THE EDB AND THE KFAC =====

      
    for(ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> >::const_iterator  corr = corrs.begin(); corr != corrs.end(); ++corr) 
    {
      
      //name of directory
      std::string base_name, name;
      prefactor pf_tmp;


      int np = corr->first.npoint.size();

      vector<XMLArray::Array<int>> mom(np); vector<int> key_row(np); vector<string> key_irrep(np);
      vector<XMLArray::Array<int>> key_mom(np); vector<string> levs(np);
      string Ei_name,Ef_name;

      for(int k = 1; k <= np; k++ ) //Array1DO has 1-based indices
      {
        key_row[k-1]   =  corr->first.npoint[k].irrep.irrep_mom.row;
        mom[k-1]   =  corr->first.npoint[k].irrep.irrep_mom.mom;

        string irp = corr->first.npoint[k].irrep.op.ops[1].name;
        std::size_t found = irp.find_last_of("_");

        key_irrep[k-1] =  irp.substr(found+1);

        size_t s = irp.find("_");
        size_t e = irp.find_first_not_of("_",s);
        levs[k-1] = irp.substr(s + 1, e);


      }

      Hadron::KeyHadronSUNNPartNPtCorr_t improv_key = corr->first;
      std::vector<ENSEM::VectorComplex> improv_val;
      if(xml_data.improv){
        improv_key.npoint[2].irrep.op.ops[1].name = xml_data.photon.name + "_" + key_irrep[1];
        improv_val = improv_corrs[improv_key];}


      //==========================================================================================
      //===== DIVIDE BY THE EXP OF SOURCE ENERGY LEVEL * t AND SINK ENERGY LEVEL * (Dt - t)  =====

      /* Find the required MassJackFile for dividing out the exponential of the energies */
      if(corr->first.npoint[1].irrep.creation_op && !corr->first.npoint[3].irrep.creation_op)
      { 
        for(int fn = 1; fn <= xml_data.Ei_jack.size(); fn++){
          if( (corr->first.npoint[1].irrep.op.ops[1].name+".jack" == xml_data.Ei_jack[fn]) ){Ei_name = xml_data.Ei_jack[fn]; break;}
        }

        for(int fn = 1; fn <= xml_data.Ef_jack.size(); fn++){
          if( (corr->first.npoint[3].irrep.op.ops[1].name+".jack" == xml_data.Ef_jack[fn]) ){Ef_name = xml_data.Ef_jack[fn]; break;}
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

        base_name = key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(key_mom[0][0]) + to_string(key_mom[0][1]) + to_string(key_mom[0][2]) + "_" + levs[0] + "x" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(key_mom[1][0]) + to_string(key_mom[1][1]) + to_string(key_mom[1][2]) + "x" +
              key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(key_mom[2][0]) + to_string(key_mom[2][1]) + to_string(key_mom[2][2]) + "_" + levs[2];  


        name = key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(mom[0][0]) + to_string(mom[0][1]) + to_string(mom[0][2]) + "_" + levs[0] + "x" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(mom[1][0]) + to_string(mom[1][1]) + to_string(mom[1][2]) + "x" +
              key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(mom[2][0]) + to_string(mom[2][1]) + to_string(mom[2][2]) + "_" + levs[2]; 

      }

      else if(!corr->first.npoint[1].irrep.creation_op && corr->first.npoint[3].irrep.creation_op)
      {

        for(int fn = 1; fn <= xml_data.Ei_jack.size(); fn++){
          if( (corr->first.npoint[3].irrep.op.ops[1].name+".jack" == xml_data.Ei_jack[fn]) ){Ei_name = xml_data.Ei_jack[fn]; break;}
        }

        for(int fn = 1; fn <= xml_data.Ef_jack.size(); fn++){
          if( (corr->first.npoint[1].irrep.op.ops[1].name+".jack" == xml_data.Ef_jack[fn]) ){Ef_name = xml_data.Ef_jack[fn]; break;}
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

        base_name = key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(key_mom[2][0]) + to_string(key_mom[2][1]) + to_string(key_mom[2][2]) + "_" + levs[0] + "x" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(key_mom[1][0]) + to_string(key_mom[1][1]) + to_string(key_mom[1][2]) + "x" +
              key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(key_mom[0][0]) + to_string(key_mom[0][1]) + to_string(key_mom[0][2]) + "_" + levs[2];


        name = key_irrep[2] + "r" + to_string(key_row[2]) + "_p" + to_string(mom[2][0]) + to_string(mom[2][1]) + to_string(mom[2][2]) + "_" + levs[0] + "x" + key_irrep[1] + "r" +
              to_string(key_row[1]) + "_p" + to_string(mom[1][0]) + to_string(mom[1][1]) + to_string(mom[1][2]) + "x" +
              key_irrep[0] + "r" + to_string(key_row[0]) + "_p" + to_string(mom[0][0]) + to_string(mom[0][1]) + to_string(mom[0][2]) + "_" + levs[2];

      }

      EnsemReal Ei_ensem; read(xml_data.elab_dir+"/"+Ei_name, Ei_ensem); EnsemReal Ef_ensem; read(xml_data.elab_dir+"/"+Ef_name, Ef_ensem); EnsemReal qsq; qsq.resize(Ei_ensem.size());
      Ei_ensem = rescaleEnsemDown( Ei_ensem );
      Ef_ensem = rescaleEnsemDown( Ef_ensem );
      double q;

      /* Find the Q^2 */
      for(size_t bin = 0; bin < Ei_ensem.size(); bin++){

        double Ef = toScalar(Ef_ensem.elem(bin));
        double Ei = toScalar(Ei_ensem.elem(bin));
        q =  pow(mom[1][0],2) + pow(mom[1][1],2) + pow(mom[1][2],2);
        qsq.elem(bin) = (pow(2*PI/(xml_data.L*xml_data.Xi), 2) * q) - pow(Ef - Ei, 2);

      }

      /* Make an ensemble of the two-point function energies */
      std::stringstream path; 
      vector<std::stringstream> kfpath(xml_data.num_matrix_elem);
      vector<EnsemComplex> kfac_ensem(xml_data.num_matrix_elem);
      double check_zero;

      path << "Q2_";
      path << std::fixed << std::setprecision(6) << toDouble( mean(qsq) ); 
      path << "/" + base_name + "/" + name + "/" ;

      for(int matrix_num = 0; matrix_num < xml_data.num_matrix_elem; matrix_num++){
        kfpath[matrix_num] << path.str() + xml_data.matrix_name[xml_data.num];
        read(kfpath[matrix_num].str(), kfac_ensem[matrix_num]);
        kfac_ensem[matrix_num] = rescaleEnsemDown(kfac_ensem[matrix_num]);
        check_zero +=  pow(toDouble( mean(real(kfac_ensem[matrix_num])) ), 2) + pow(toDouble( mean(imag(kfac_ensem[matrix_num])) ), 2);
        
      }

      if(( check_zero == 0.0) ){continue;} // rule out corrs with kfac = 0.0


      if(corr->second.size() != Ei_ensem.size()){ 
        cout << "Ei size and the edb ensem size don't match" << endl;
        exit(1);
      }
      
      if(corr->second.size() != Ef_ensem.size()){ 
        cout << "Ef size and the edb ensem size don't match" << endl;
        exit(1);
      }

      if(xml_data.improv){if(corr->second.size() != improv_val.size()){ 
        cout << "the edb sizes of corr and improv don't match" << endl;
        exit(1);
      }}

        
      ThreePtKeys::KeyHadronSUNNPartNPtIrrep_t tmp_k(key_row, key_mom, key_irrep);

      if(key.find(tmp_k) == key.end()){
        key.insert(make_pair(tmp_k,key_count)); corr_num = key_count; edb_data.dir.insert(make_pair(corr_num,base_name)); key_count++;}

      else{corr_num = key.find(tmp_k)->second;}


      if(size_of_base.find(corr_num) == size_of_base.end()){size_of_base.insert(make_pair(corr_num,1));}
      else{size_of_base.find(corr_num)->second += 1;}


      int dt = std::abs(corr->first.npoint[1].t_slice - corr->first.npoint[3].t_slice ); 

      EnsemVectorComplex yt_post; yt_post.resize(corr->second.size()); yt_post.resizeObs(corr->second[0].numElem()); //to write the files 
      EnsemVectorComplex yt_pre; yt_pre.resize(corr->second.size()); yt_pre.resizeObs(corr->second[0].numElem()); //to write the files  

      EnsemVectorComplex dyt_post; dyt_post.resize(corr->second.size()); dyt_post.resizeObs(corr->second[0].numElem()); //to write the improv files 
      EnsemVectorComplex dyt_pre; dyt_pre.resize(corr->second.size()); dyt_pre.resizeObs(corr->second[0].numElem()); //to write the improv files  

      for(int t = 0; t < corr->second[0].numElem(); t++){ //numElem gives the number of observables in the ensem so use any index for convinience use 0

        ENSEM::EnsemReal yt_rl; yt_rl.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 
        ENSEM::EnsemReal yt_im; yt_im.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 
        ENSEM::EnsemComplex yt; yt.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 

        ENSEM::EnsemReal dyt_rl; dyt_rl.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 
        ENSEM::EnsemReal dyt_im; dyt_im.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 
        ENSEM::EnsemComplex dyt; dyt.resize(corr->second.size()); //resize the ensem to the size of the vector = number of bins 

        for(size_t bin = 0; bin < corr->second.size(); bin++){//corr num issue

          yt_rl.elem(bin) = toScalar(real(peekObs(corr->second[bin], t))); 
          yt_im.elem(bin) = toScalar(imag(peekObs(corr->second[bin], t))); 

          ENSEM::RComplex<double> tmp(toScalar(yt_rl.elem(bin)),toScalar(yt_im.elem(bin))); //for some reason cannot read in a complex corr have to stick to this method
          yt.elem(bin)    = tmp;

          if(xml_data.improv){//if improv is present
          dyt_rl.elem(bin) = toScalar(real(peekObs(improv_val[bin], t))); 
          dyt_im.elem(bin) = toScalar(imag(peekObs(improv_val[bin], t))); 

          ENSEM::RComplex<double> dtmp(toScalar(dyt_rl.elem(bin)),toScalar(dyt_im.elem(bin))); //for some reason cannot read in a complex corr have to stick to this method
          dyt.elem(bin)    = dtmp;            
          }
        }

        pokeObs(yt_pre, yt, t);

        yt_rl = rescaleEnsemDown(yt_rl);
        yt_im = rescaleEnsemDown(yt_im);

        if(xml_data.improv){
          pokeObs(dyt_pre, dyt, t);

          dyt_rl = rescaleEnsemDown(dyt_rl);
          dyt_im = rescaleEnsemDown(dyt_im);          
        }



        for(size_t bin = 0; bin < corr->second.size(); bin++){//corr num issue


          double Ef, Ei;
          complex<double> ff(toScalar(yt_rl.elem(bin)),toScalar(yt_im.elem(bin))) ;

          Ef = toScalar(Ef_ensem.elem(bin)); 
          Ei = toScalar(Ei_ensem.elem(bin));

          if(xml_data.improv){
            complex<double> dff(toScalar(dyt_rl.elem(bin)),toScalar(dyt_im.elem(bin))) ;
            complex<double> cf;
            cf = (-2.43787/4.0) * (Ef - Ei);
            ff += (cf*dff);}


          if(xml_data.divkfac){complex<double> kf = toScalar(kfac_ensem[0].elem(bin)); ff = ff/kf;}

          for(int ph= 1; ph <= xml_data.rephase.size(); ph++ ){if(xml_data.rephase[ph] == Ei_name || xml_data.rephase[ph] == Ef_name){ ff = -1.0*ff;}}
          
          yt_rl.elem(bin) = std::exp(Ef*(dt - t))/(sqrt(xml_data.CurrFac)) * std::exp(Ei*t) * real(ff); //multiplying by the exp

          yt_im.elem(bin) = std::exp(Ef*(dt - t))/(sqrt(xml_data.CurrFac)) * std::exp(Ei*t) * imag(ff); //multiplying by exp 
             

          ENSEM::RComplex<double> tmp(toScalar(yt_rl.elem(bin)),toScalar(yt_im.elem(bin))); //for some reason cannot read in a complex corr have to stick to this method
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

      if(xml_data.improv){
      {
        ostringstream outfile; outfile << path.str() << "/improv_corr.jack";
        write(outfile.str(), dyt_pre);
      }

	}

      pf_tmp.ei   = rescaleEnsemUp(Ei_ensem); //if the vector is the creation operator
      pf_tmp.ef   = rescaleEnsemUp(Ef_ensem);
      pf_tmp.qsq  = rescaleEnsemUp(qsq);
      if(edb_data.pref.find(corr_num) == edb_data.pref.end()){edb_data.pref.insert(make_pair(corr_num, pf_tmp));}


      if(corr_num >= edb_data.x.size()){
        x_vec.push_back(x_t); y_vec_rl.push_back(ytensem_rl); y_vec_im.push_back(ytensem_im);
        edb_data.x.resize(corr_num+1, x_vec); edb_data.y_ensem_rl.resize(corr_num+1, y_vec_rl); edb_data.y_ensem_im.resize(corr_num+1, y_vec_im);
      }

      else{
        x_vec = edb_data.x.at(corr_num); y_vec_rl = edb_data.y_ensem_rl.at(corr_num); y_vec_im = edb_data.y_ensem_im.at(corr_num);
        for(size_t i = 0; i < x_vec.size(); i++){
          if( dt == (x_vec.at(i)).at(0).first ){
            double sizeo = (size_of_base.find(corr_num)->second - 1);
            double sizen = size_of_base.find(corr_num)->second;

            for(size_t tin = 0; tin < ytensem_rl.size(); tin++){
              y_vec_rl.at(i).at(tin) =  ((toScalar(sizeo) * rescaleEnsemDown(y_vec_rl.at(i).at(tin)) ) + rescaleEnsemDown(ytensem_rl.at(tin)))/toScalar(sizen);
              y_vec_im.at(i).at(tin) =  ((toScalar(sizeo) * rescaleEnsemDown(y_vec_im.at(i).at(tin)) ) + rescaleEnsemDown(ytensem_im.at(tin)))/toScalar(sizen);

              y_vec_rl.at(i).at(tin) = rescaleEnsemUp(y_vec_rl.at(i).at(tin));
              y_vec_im.at(i).at(tin) = rescaleEnsemUp(y_vec_im.at(i).at(tin));

            }

            break;
          }

          else{continue;}
        }

        edb_data.y_ensem_rl.at(corr_num) = y_vec_rl; edb_data.y_ensem_im.at(corr_num) = y_vec_im; 
      }

      
      key_row.clear(); mom.clear(); key_irrep.clear(); x_vec.clear(); y_vec_rl.clear(); y_vec_im.clear();
      x_t.clear();  ytensem_rl.clear(); ytensem_im.clear();


    }

    corrs.clear();

  }


  return edb_data;


};

 //************************************************************************************************
 
//!Utility functions
void removeSubstrs(std::string& s, const std::string& p) {

   std::string::size_type n = p.length();

   for (std::string::size_type i = s.find(p);
        i != std::string::npos;
        i = s.find(p))
      s.erase(i, n);
};
