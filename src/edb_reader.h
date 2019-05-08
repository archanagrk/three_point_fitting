#ifndef __EDB_READER_H__
#define __EDB_READER_H__

/*! \file
 * \brief Read a correlator edb
 */

#include <AllConfStoreDB.h>

#include <io/key_val_db.h>
#include <adat/map_obj.h>
#include <hadron/hadron_sun_npart_npt_corr.h>


/* std library to create dir */
#include <filesystem>

using namespace std;
//namespace fs = std::filesystem;

namespace edb_reader{

    /* functions */
    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > getCorrs(const std::string& dbase);

}


void removeSubstrs(std::string& s, const std::string& p);



#endif