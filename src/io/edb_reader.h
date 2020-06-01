#ifndef __EDB_READER_H__
#define __EDB_READER_H__

/*! \file
 * \brief Read a correlator edb
 */

//adat

#include <AllConfStoreDB.h>
#include <io/key_val_db.h>
#include <adat/map_obj.h>
#include <hadron/hadron_sun_npart_npt_corr.h>

// semble
#include"semble/semble_file_management.h"
#include "semble/semble_meta.h"


#include "xml_reader.h"
#include "lib/rotations.h"


using namespace SEMBLE;

/* structs */
 //************************************************************************************************
struct EDBdata_t{
        
        vector<vector< vector<ENSEM::EnsemReal> >> y_ensem_rl, y_ensem_im;
        std::map<int, string> dir;
        std::map<int, prefactor> pref;
        vector<vector<vector<pair<int,int>>> >x; 
};


/* functions */
 //************************************************************************************************

ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > getCorrs(const std::string& dbase);


 //************************************************************************************************
// routine to read the corrs from edb files //
EDBdata_t read_corrs(XMLinidata_t& xml_data);

 //************************************************************************************************
void removeSubstrs(std::string& s, const std::string& p);



#endif