#include "xml_reader.h"

// routine to read the input xml //
XMLinidata_t read_xml_ini(XMLReader& xml_in){

  XMLinidata_t xml_data;

  try
      {
        read(xml_in,"/ThreeptIniParams/inputProps/NumdbFiles",xml_data.num);
        read(xml_in,"/ThreeptIniParams/NPoints/Improvement/active",xml_data.improv);

        /* Have to sort the elems if not already descending order in Dt */
        std::map<int,int> sort_dt;
        int dt, j, tmin, tmax, tsrc, tsnk, Nt_min;
        double Ei_min, Ef_min, Ei_max, Ef_max;

        for(int i = 1; i <= xml_data.num; i++ ){ 
          read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(i)+"]/dt", dt);
          sort_dt.insert(make_pair(dt,i));
        }


        string file, imp;

        for(std::map<int,int>::iterator it=sort_dt.begin(); it!=sort_dt.end(); ++it){ 
          j = it->second;
          read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/name", file);
          if(xml_data.improv){read(xml_in, "/ThreeptIniParams/inputProps/dbFnames/elem["+std::to_string(j)+"]/improv", imp);}
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

          xml_data.filename.push_back(file); if(xml_data.improv){xml_data.improvname.push_back(imp);} 
          xml_data.dt_v.push_back(dt); xml_data.tmin_v.push_back(tmin); xml_data.tmax_v.push_back(tmax); xml_data.tsrc_v.push_back(tsrc); xml_data.tsnk_v.push_back(tsnk);
          xml_data.tmin_max_v.push_back(tmin_max); xml_data.Nt_min_v.push_back(Nt_min); xml_data.Ei_min_v.push_back(Ei_min); xml_data.Ef_min_v.push_back(Ef_min);
          xml_data.Ei_max_v.push_back(Ei_max); xml_data.Ef_max_v.push_back(Ef_max);
        }

        read(xml_in,"/ThreeptIniParams/FitProps/fit_qual/type",xml_data.qual_type);
        read(xml_in,"/ThreeptIniParams/FitProps/chisq_cutoff",xml_data.chisq_cut);
        read(xml_in,"/ThreeptIniParams/FitProps/noise_cutoff",xml_data.noise);
        read(xml_in,"/ThreeptIniParams/FitProps/divKfac",xml_data.divkfac);
        read(xml_in,"/ThreeptIniParams/FitProps/kfac",xml_data.kfac_xml);
        read(xml_in,"/ThreeptIniParams/FitProps/Ei",xml_data.Ei_jack);      
        read(xml_in,"/ThreeptIniParams/FitProps/Ef",xml_data.Ef_jack);
        read(xml_in,"/ThreeptIniParams/FitProps/Rephase",xml_data.rephase);;
        read(xml_in, "/ThreeptIniParams/FitProps/ElabDir", xml_data.elab_dir);

        read(xml_in,"/ThreeptIniParams/Renormalization/Z",xml_data.Z_file);
        read(xml_in,"/ThreeptIniParams/Renormalization/InvCurrFac2",xml_data.CurrFac);
        read(xml_in,"/ThreeptIniParams/Renormalization/writeZFile",xml_data.Zfile);


        string lorentz_index;
        vector<pair<string,Array1dO<string>>> threept(3);
        
        read(xml_in,"/ThreeptIniParams/NPoints/LorentzIndex",lorentz_index);

        for(int i = 1; i <= 3; i++){
          string tmp_op; Array1dO<string> tmp_lev;
          read(xml_in,"/ThreeptIniParams/NPoints/Ops/elem["+std::to_string(i)+"]/name",tmp_op);
          read(xml_in,"/ThreeptIniParams/NPoints/Ops/elem["+std::to_string(i)+"]/levels",tmp_lev);
          threept[i-1] = make_pair(tmp_op,tmp_lev);      
        }

        if(xml_data.improv){

          read(xml_in,"/ThreeptIniParams/NPoints/Improvement/name",xml_data.photon.name);

          read(xml_in,"/ThreeptIniParams/NPoints/Improvement/levels",xml_data.photon.levels);        
          read(xml_in,"/ThreeptIniParams/NPoints/Improvement/nu",xml_data.photon.nu_s);
          read(xml_in,"/ThreeptIniParams/NPoints/Improvement/xi_0",xml_data.photon.xi0);
          read(xml_in,"/ThreeptIniParams/NPoints/Improvement/a_s",xml_data.photon.a_s);

        }

        read(xml_in,"/ThreeptIniParams/LattProps/Xi", xml_data.Xi);
        read(xml_in,"/ThreeptIniParams/LattProps/XiError", xml_data.XiE);
        read(xml_in,"/ThreeptIniParams/LattProps/L", xml_data.L);


      }
    catch( const string& error ){
      cerr << "Error reading input file : " << error << endl;
      }


    //=============================
    //===== READ THE KFAC XML =====

      if(xml_data.divkfac){

      XMLReader xml_kf_in(xml_data.kfac_xml);
      try
        {
          read(xml_kf_in, "/kfac/pts", xml_data.pts);    
          read(xml_kf_in,"/kfac/kfacFile",xml_data.matrix_name);
          xml_data.num_matrix_elem = xml_data.matrix_name.size();
      
        }

      catch( const string& error ){
        cerr << "Error reading input file : " << error << endl;
        }
      }

  return xml_data;

};

