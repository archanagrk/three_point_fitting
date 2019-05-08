#ifndef __EDB_READER_H__
#define __EDB_READER_H__
/*! \file
 * \brief Generic operations on a ensem db file
 */


#include "adat/singleton.h"
#include "adat/funcmap.h"

#include "io/adat_xmlio.h"
#include "hadron/ensem_filenames.h"
#include "io/key_val_db.h"
#include "hadron/dbtype.h"

#include "AllConfStoreDB.h"
#include "DBString.h"

#include "formfac/hadron_1pt_corr.h"
#include "formfac/hadron_2pt_corr.h"
#include "formfac/hadron_3pt_corr.h"
#include "formfac/discoblocks.h"

#include "hadron/hadron_npart_npt_conn_graph.h"
#include "hadron/hadron_sun_npart_npt_corr.h"
#include "hadron/su2_corr/hadron_npart_npt_corr.h"
#include "hadron/hadron_npart_npt_conn_graph.h"
#include "hadron/meson_elem_type.h"
#include "hadron/prop_elem_type.h"
#include "hadron/hadron_timeslice.h"

using namespace std;
using namespace ENSEM;
using namespace ADATXML;
using namespace ADATIO;
using namespace Util;
using namespace FILEDB;
using namespace FF;
using namespace Hadron;

namespace edb_reader
{
    int usage (int argc, char** argv);
    void printMeta(const string& dbase);

    template<typename K>
    void printKeys(const string& dbase);

    template<typename K>
    void printKeysXML(const string& xml_file, const string& dbase);

    template<typename K, typename V> V 
    printKeyValue(const K& ky, AllConfStoreDB< SerialDBKey<K>,  SerialDBData<typename EnsemScalar<V>::Type_t> >& database);

    template<typename K, typename V>
    void printEnsem(const string& xml_key_file, const string& dbase);

    void printMeson(const string& xml_key_file, const string& dbase);
    void printProp(const string& xml_key_file, const string& dbase);

    template<typename K, typename V>
    int workHandler(int argc, char** argv);

    template<typename K, typename V>
    int mesonWorkHandler(int argc, char** argv);


    template<typename K, typename V>
    int propWorkHandler(int argc, char** argv);

    struct DumbDisambiguator {};
}


//! Hold possible functions
typedef SingletonHolder< 
    FunctionMap<edb_reader::DumbDisambiguator,
        int,
        std::string,
        TYPELIST_2(int, char**),
        int (*)(int, char**),
        StringFunctionMapError> >
TheKeyHandlerObjFuncMap;

namespace edb_reader
{
    //! Registration hack
    bool registerAll();
}



#endif
