#ifndef __FF_QSQ_READER_H__
#define __FF_QSQ_READER_H__

/*! \file
 * \brief Read the qsq list file
 */

#include<stdio.h>
#include <iostream> 
#include <iterator> 
#include <map> 
#include <vector>
#include <fstream>
#include <sstream>
#include "ensem/ensem.h"


using namespace std;


/* functions */

namespace readlist{
    map<double ,ENSEM::EnsemReal> make_ff_qsq(const string& listfile);
}



#endif