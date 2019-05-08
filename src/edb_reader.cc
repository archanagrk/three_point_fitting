
#include "edb_reader.h"

//! Get corrs

ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, std::vector<ENSEM::VectorComplex> > edb_reader::getCorrs(const std::string& dbase)
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


void removeSubstrs(std::string& s, const std::string& p) {

   std::string::size_type n = p.length();

   for (std::string::size_type i = s.find(p);
        i != std::string::npos;
        i = s.find(p))
      s.erase(i, n);
}




