#ifndef __KEY_STRUCT_H__
#define __KEY_STRUCT_H__

// -*- C++ -*-
/*! \file
 * \A Key to compare 3pt functions
 */

#include <complex>
#include "ensem/ensem.h"



using namespace std;


namespace key_struct
{
    using namespace ENSEM;

    /* compare function for the key */
    //bool CompAr(XMLArray::Array<int> a1, XMLArray::Array<int> a2);

    /*struct */
    struct KeyHadronSUNNPartNPtIrrep_t
    {
        KeyHadronSUNNPartNPtIrrep_t() : row(0) {}
        KeyHadronSUNNPartNPtIrrep_t(vector<int>  row_, vector<XMLArray::Array<int>>& mom_, vector<string>& irrep_) : row(row_), mom(mom_), irrep(irrep_) {}

        vector<int>                     row;
        vector<XMLArray::Array<int>>    mom;
        vector<string>                  irrep;
    };


    //---------------------------------
    //! Used for error output
    std::ostream& operator<<(std::ostream& os, const KeyHadronSUNNPartNPtIrrep_t& op);

    //----------------------------------------------------------------------------
    //! Read a key
    void read(XMLReader& xml, const std::string& path, KeyHadronSUNNPartNPtIrrep_t& param);

    //! Write a key
    void write(XMLWriter& xml, const std::string& path, const KeyHadronSUNNPartNPtIrrep_t& param);

    //---------------------------------
    //! Binary reader
    void read(BinaryReader& bin, KeyHadronSUNNPartNPtIrrep_t& param);

    //! Binary writer
    void write(BinaryWriter& bin, const KeyHadronSUNNPartNPtIrrep_t& param);

} //end namespace key_struct

//************************************************************************

namespace std
{

    inline bool CompAr(XMLArray::Array<int> a1, XMLArray::Array<int> a2){
        int lhs = 0; int rhs = 0;
        for(int i = 0; i < a1.size(); i++){
            lhs += a1[a1.size() - 1 -i]*(pow(10,i));
            rhs += a2[a1.size() - 1 -i]*(pow(10,i)); 
        }
        if(lhs < rhs){return true;}
        
        return false;
    }

    template<> struct less<key_struct::KeyHadronSUNNPartNPtIrrep_t>
    {
       bool operator() (const key_struct::KeyHadronSUNNPartNPtIrrep_t& lhs, const key_struct::KeyHadronSUNNPartNPtIrrep_t& rhs) const
       {

        if(std::CompAr(lhs.mom[0],rhs.mom[0])){return true;}
        else if(std::CompAr(rhs.mom[0],lhs.mom[0])){return false;}

        else if(std::CompAr(lhs.mom[1],rhs.mom[1])){return true;}
        else if(std::CompAr(rhs.mom[1],lhs.mom[1])){return false;}

        else if(std::CompAr(lhs.mom[2],rhs.mom[2])){return true;}
        else if(std::CompAr(rhs.mom[2],lhs.mom[2])){return false;}

        else if(lhs.row[0] < rhs.row[0]){return true;}
        else if(rhs.row[0] < lhs.row[0]){return false;}

        else if(lhs.row[1] < rhs.row[1]){return true;}
        else if(rhs.row[1] < lhs.row[1]){return false;}

        else if(lhs.row[2] < rhs.row[2]){return true;}
        else if(rhs.row[2] < lhs.row[2]){return false;}

        else if(lhs.irrep[0] < rhs.irrep[0]){return true;}
        else if(rhs.irrep[0] < lhs.irrep[0]){return false;}

        else if(lhs.irrep[1] < rhs.irrep[1]){return true;}
        else if(rhs.irrep[1] < lhs.irrep[1]){return false;}

        else if(lhs.irrep[2] < rhs.irrep[2]){return true;}
        else if(rhs.irrep[2] < lhs.irrep[2]){return false;}


        else{return false;}


       }
    };
}
//************************************************************************
 /* STRUCT FOR HOLDING THE PREFACTOR */
//************************************************************************

namespace key_struct{

    struct prefactor{
        ENSEM::EnsemReal qsq;
        ENSEM::EnsemComplex kfac;
        ENSEM::EnsemReal ei;
        ENSEM::EnsemReal ef;
    };
};

//************************************************************************

#endif

