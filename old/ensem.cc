
#include "edb_reader.h"

//! Get corrs

ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, ENSEM::EnsemVectorComplex > edb_reader::getCorrs(const std::string& dbase)
{
  typedef Hadron::KeyHadronSUNNPartNPtCorr_t              K;
  typedef ENSEM::VectorComplex                            V;
  typedef ENSEM::EnsemVectorComplex                       EV;


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
  ADAT::MapObject<K, EV > dest;

  for(int i = 0; i < keys.size(); ++i) {
    EV val;
    val.resize(vals[i].size());
    val.resizeObs(vals[i][0].data().elem().size());

    //for(int n=0; n < vals[i].size(); ++n) { 
        //val.elem(i) = vals[i][n].data().elem();
    //}

    dest.insert(keys[i].key(), val);
  }

  return dest;

}


edb_reader::KeyHadronSUNNPartNPtIrrep_t::KeyHadronSUNNPartNPtIrrep_t(vector<Hadron::KeyCGCIrrepMom_t> key_in){
    vector<int> tmp_row;
    vector<XMLArray::Array<int>> tmp_mom;

    for(int i = 0; i < key_in.size(); i++){
        tmp_row.push_back(key_in[i].row); tmp_mom.push_back(key_in[i].mom);
    }

    mom = tmp_mom;
    row = tmp_row;
}


void removeSubstrs(std::string& s, const std::string& p) {

   std::string::size_type n = p.length();

   for (std::string::size_type i = s.find(p);
        i != std::string::npos;
        i = s.find(p))
      s.erase(i, n);
}



//std::string ensemFileName(const Hadron::KeyHadronSUNNPartNPtCorr_t& key);

