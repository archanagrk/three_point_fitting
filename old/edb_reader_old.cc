#include "edb_reader.h"

/* Support for EDBs. Most of the functions are from Robert's dbutil with minor tweaks */

namespace edb_reader
{
  int usage (int argc, char* argv[])
  {
    cerr << "Usage: " << argv[0] << " <dbasename> <operation> [<operation args> ...]" << endl;
    return -1;
  }


    //! Print the metadata
    void printMeta(const string& dbase)
    {
    // Open DB solely to get the user data
    ConfDataStoreDB<StringKey,UserData> database;
    if (database.open(dbase, O_RDONLY, 0400) != 0)
    {
        std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
        exit(1);
    }

    std::string meta_str;
    database.getUserdata(meta_str);

    std::cout << __func__ << ": metadata:\n" << meta_str << std::endl;
    }


    //! Print all the keys
    template<typename K>
    void printKeys(const string& dbase)
    {
    // Open DB
    AllConfStoreDB< SerialDBKey<K>, UserData> database;
    if (database.open(dbase, O_RDONLY, 0400) != 0)
    {
        std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
        exit(1);
    }

    std::vector< SerialDBKey<K> > keys;
    database.keys(keys);

    for(int i=0; i < keys.size(); ++i)
    {
        std::cout << keys[i].key() << "\n";
    }
    }


    //! Print all the keys to a xml file
    template<typename K>
    void printKeysXML(const string& xml_file, const string& dbase)
    {
    // Open DB
    AllConfStoreDB< SerialDBKey<K>, UserData> database;
    if (database.open(dbase, O_RDONLY, 0400) != 0)
    {
        std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
        exit(1);
    }

    std::vector< SerialDBKey<K> > keys;
    database.keys(keys);

    // Print keys
    XMLFileWriter xml_out(xml_file);
    push(xml_out, "Keys");
    for(int i=0; i < keys.size(); ++i)
    {
        write(xml_out, "elem", keys[i].key());
    }
    pop(xml_out);

    xml_out.close();
    }


    //! Get a key/value
    template<typename K, typename V>
    V printKeyValue(const K& ky, 
            AllConfStoreDB< SerialDBKey<K>,  SerialDBData<typename EnsemScalar<V>::Type_t> >& database)
    {
    typedef typename EnsemScalar<V>::Type_t SV;

    SerialDBKey<K> key;
    key.key() = ky;

    std::vector< SerialDBData<SV> > vals;
    int ret;
    if ((ret = database.get(key, vals)) != 0)
    {
        std::cerr << __func__ << ": key not found\n" << ky;
        exit(1);
    }

    V eval;
    eval.resize(vals.size());
    eval.resizeObs(vals[0].data().numElem());

    for(int i=0; i < vals.size(); ++i)
    {
        SV sval = vals[i].data();
        pokeEnsem(eval, sval, i);
    }

    return eval;
    }


    //! Print ensemble files
    template<typename K, typename V>
    void printEnsem(const string& xml_key_file, const string& dbase)
    {
    typedef typename EnsemScalar<V>::Type_t SV;

    // Open DB
    AllConfStoreDB< SerialDBKey<K>,  SerialDBData<SV> > database;
    if (database.open(dbase, O_RDONLY, 0400) != 0)
    {
        std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
        exit(1);
    }

    try
    {
        XMLReader xml_in(xml_key_file);
        Array<K> keys;
        read(xml_in, "/Keys", keys);

        // Write ensemble file
        for(int i=0; i < keys.size(); ++i)
        {
        write(ensemFileName(keys[i]), printKeyValue<K,V>(keys[i], database));
        }
    }
    catch(const std::string& e) 
    {
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
    }
    catch(std::exception& e) 
    {
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
    }
    }


    //! Print meson elemental
    void printMeson(const string& xml_key_file, const string& dbase)
    {
    typedef KeyMesonElementalOperator_t K;
    typedef ValMesonElementalOperator_t V;
    
    try
    {
        // Open DB
        ConfDataStoreDB< SerialDBKey<K>,  SerialDBData<V> > meson_db;
        if (meson_db.open(dbase, O_RDONLY, 0400) != 0)
        {
        std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
        exit(1);
        }

        std::string meta_data;
        if (meson_db.getUserdata(meta_data) != 0)
        {
        std::cerr << __func__ << ": some error reading meta data out of db=" << dbase << std::endl;
        exit(1);
        }

        // Read the desired keys
        XMLReader xml_in(xml_key_file);
        std::vector<K> keys;
        read(xml_in, "/Keys", keys);

        // Write binary file
        for(auto k = keys.begin(); k != keys.end(); ++k)
        {
        // The key
        SerialDBKey<KeyMesonElementalOperator_t> key(*k);

        // Find the object
        SerialDBData<ValMesonElementalOperator_t> val;
        if (meson_db.get(key, val) != 0)
        {
        std::cerr << __func__ << ": key not found: key=" << key.key() << std::endl;
        exit(1);
        }

        // Build filename
        std::string file_name;

        if (1)
        {
        std::ostringstream os;
        os << "meson";
        os << ".t" << key.key().t_slice;
        os << ".d";
        for(int i = 0; i < key.key().displacement.size(); ++i)
        {
        os << key.key().displacement[i];
        }

        os << ".p";
        for(int i = 0; i < key.key().mom.size(); ++i)
        {
        os << key.key().mom[i];
        }

        os << ".dat";
        
        file_name = os.str();
        }

        // Write the file
        BinaryFileWriter bin(file_name);

        int nrows = val.data().op.nrows();
        int ncols = val.data().op.nrows();

        write(bin, nrows);
        write(bin, ncols);

        // Write the size
        for(int cl = 0; cl < val.data().op.nrows(); ++cl)
        {
        for(int cr = 0; cr < val.data().op.ncols(); ++cr)
        {
        write(bin, val.data().op(cl,cr));
        } // cr
        } // cl
        } // key
    }
    catch(const std::string& e) 
    {
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
    }
    catch(std::exception& e) 
    {
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
    }
    }


    //! Print prop elemental
    void printProp(const string& xml_key_file, const string& dbase)
    {
    typedef KeyPropElementalOperator_t K;
    typedef ValPropElementalOperator_t V;
    
    try
    {
        // Open DB
        ConfDataStoreDB< SerialDBKey<K>,  SerialDBData<V> > prop_db;
        if (prop_db.open(dbase, O_RDONLY, 0400) != 0)
        {
        std::cerr << __func__ << ": error opening dbase= " << dbase << std::endl;
        exit(1);
        }

        std::string meta_data;
        if (prop_db.getUserdata(meta_data) != 0)
        {
        std::cerr << __func__ << ": some error reading meta data out of db=" << dbase << std::endl;
        exit(1);
        }

        // Read the desired keys
        XMLReader xml_in(xml_key_file);
        std::vector<K> keys;
        read(xml_in, "/Keys", keys);

        // Write binary file
        for(auto k = keys.begin(); k != keys.end(); ++k)
        {
        // The key
        SerialDBKey<KeyPropElementalOperator_t> key(*k);

        // Find the object
        SerialDBData<ValPropElementalOperator_t> val;
        if (prop_db.get(key, val) != 0)
        {
        std::cerr << __func__ << ": key not found: key=" << key.key() << std::endl;
        exit(1);
        }

        // Build filename
        std::string file_name;

        if (1)
        {
        std::ostringstream os;
        os << "prop";
        os << ".m_" << key.key().mass_label;
        os << ".t_" << key.key().t_slice;
        os << ".t0_" << key.key().t_source;
        os << ".sl_" << key.key().spin_l;
        os << ".sr_" << key.key().spin_r;
        os << ".dat";
        
        file_name = os.str();
        }

        // Write the file
        BinaryFileWriter bin(file_name);

        int nrows = val.data().mat.nrows();
        int ncols = val.data().mat.nrows();

        // Write the size
        write(bin, nrows);
        write(bin, ncols);

        // Write the objects
        for(int cl = 0; cl < val.data().mat.nrows(); ++cl)
        {
        for(int cr = 0; cr < val.data().mat.ncols(); ++cr)
        {
        write(bin, val.data().mat(cl,cr));
        } // cr
        } // cl
        } // key
    }
    catch(const std::string& e) 
    {
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
    }
    catch(std::exception& e) 
    {
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
    }
    }


    //! Handle all the work
    template<typename K, typename V>
    int workHandler(int argc, char* argv[])
    {
    std::string dbase = argv[0];
    std::string op = argv[1];

    // A simplified version of handler
    if (op == "meta")
    {
        printMeta(dbase);
    }
    else if (op == "keys")
    {
        printKeys<K>(dbase);
    }
    else if (op == "keysxml")
    {
        if (argc != 3)
        {
        return edb_reader::usage(argc, argv);
        }

        printKeysXML<K>(argv[2], dbase);
    }
    else if (op == "get")
    {
        if (argc != 3)
        {
        return edb_reader::usage(argc, argv);
        }

        printEnsem<K,V>(argv[2], dbase);
    }
    else
    {
        std::cerr << argv[0] << ": unsupported operation= " << op << endl;
        return -1;
    }

    return 0;
    }


    //! Handle all the work
    template<typename K, typename V>
    int mesonWorkHandler(int argc, char* argv[])
    {
    std::string dbase = argv[0];
    std::string op = argv[1];

    // A simplified version of handler
    if (op == "meta")
    {
        printMeta(dbase);
    }
    else if (op == "keys")
    {
        printKeys<K>(dbase);
    }
    else if (op == "keysxml")
    {
        if (argc != 3)
        {
        return edb_reader::usage(argc, argv);
        }

        printKeysXML<K>(argv[2], dbase);
    }
    else if (op == "get")
    {
        if (argc != 3)
        {
        return edb_reader::usage(argc, argv);
        }

        printMeson(argv[2], dbase);
    }
    else
    {
        std::cerr << argv[0] << ": unsupported operation= " << op << endl;
        return -1;
    }

    return 0;
    }


    //! Handle all the work
    template<typename K, typename V>
    int propWorkHandler(int argc, char* argv[])
    {
    std::string dbase = argv[0];
    std::string op = argv[1];

    // A simplified version of handler
    if (op == "meta")
    {
        printMeta(dbase);
    }
    else if (op == "keys")
    {
        printKeys<K>(dbase);
    }
    else if (op == "keysxml")
    {
        if (argc != 3)
        {
        return edb_reader::usage(argc, argv);
        }

        printKeysXML<K>(argv[2], dbase);
    }
    else if (op == "get")
    {
        if (argc != 3)
        {
        return edb_reader::usage(argc, argv);
        }

        printProp(argv[2], dbase);
    }
    else
    {
        std::cerr << argv[0] << ": unsupported operation= " << op << endl;
        return -1;
    }

    return 0;
    }

    bool registerAll()
    {
        // Register all the functions
        bool success = true; 

        // Register the functions
        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("hadron1Pt"), 
                                        edb_reader::workHandler<KeyHadron1PtCorr_t,EnsemEnsemReal>);
        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("hadron2Pt"), 
                                        edb_reader::workHandler<KeyHadron2PtCorr_t,EnsemEnsemReal>);
        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("hadron3Pt"), 
                                        edb_reader::workHandler<KeyHadron3PtCorr_t,EnsemEnsemReal>);
        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("hadronNPartNPtCorr"), 
                                        edb_reader::workHandler<KeyHadronNPartNPtCorr_t,EnsemEnsemReal>);
        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("hadronSUNNPartNPtCorr"), 
                                        edb_reader::workHandler<KeyHadronSUNNPartNPtCorr_t,EnsemEnsemReal>);
        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("hadronNPartNPtConnGraph"), 
                                        edb_reader::workHandler<KeyHadronNPartNPtConnGraph_t,EnsemEnsemReal>);

        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("mesonElemOp"), 
                                        edb_reader::mesonWorkHandler<KeyMesonElementalOperator_t,ValMesonElementalOperator_t>);
        success &= TheKeyHandlerObjFuncMap::Instance().registerFunction(string("propElemOp"), 
                                        edb_reader::propWorkHandler<KeyPropElementalOperator_t,ValPropElementalOperator_t>);
        
        return success;
    }    

}
