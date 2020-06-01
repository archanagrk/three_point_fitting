#ifndef __XML_READER__
#define __XML_READER__


// adat
#include "hadron/ensem_filenames.h"



#include "key_struct.h"


using namespace ThreePtKeys;
using namespace std;

 //************************************************************************************************
// structs

 struct improv_props
 {
   string name;
   ADAT::Array1dO<string> levels;
   ENSEM::Complex coeff;
   double nu_s, a_s;
   double xi0;
 };

struct XMLinidata_t{

        int num; 
        bool improv; 
        vector<string> filename; 
        vector<string> improvname;
        string elab_dir;
        string kfac_xml;
        bool divkfac;
        improv_props photon;
        int pts;
        int num_matrix_elem;
        ADAT::Array1dO<string> matrix_name;
        

        vector<int> dt_v, tmin_v, tmax_v, tsrc_v, tsnk_v, Nt_min_v;
        vector<double> Ei_min_v, Ef_min_v, Ei_max_v, Ef_max_v;
        vector<std::tuple<int,int,int>> tmin_max_v;
        ADAT::Array1dO<string> Ei_jack, Ef_jack, rephase;

        double noise, chisq_cut;
        string qual_type;

        double CurrFac;
        double Zt;
        double Zs;
        string Z_file;
        string Zfile;

        double Xi;
        double XiE; 
        int L;

};



 //************************************************************************************************
// routine to read the input xml //
XMLinidata_t read_xml_ini(XMLReader& xml_in);


#endif