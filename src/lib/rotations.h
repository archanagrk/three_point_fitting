#ifndef __ROTMATRX__
#define __ROTMATRX__

//std lib
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "math.h"
#include <stdio.h>
#include <vector>

//itpp
#include <itpp/stat/misc_stat.h>
#include <itpp/base/timing.h>
#include <itpp/itbase.h>

//adat
#include <adat/handle.h>
#include "hadron/irreps_cubic_factory.h"
#include "hadron/irreps_cubic_helicity_factory.h"
#include "hadron/irrep_util.h" //adat_devel
#include "ensem/ensem.h"


using namespace std;



namespace Rot {
  
/* rotation matrix for vectors */
itpp::mat EulerRot_t(double alpha, double beta, double gamma); 
XMLArray::Array<int> EulerRotVec_t(Hadron::CubicCanonicalRotation_t ref_angles_z, Hadron::CubicCanonicalRotation_t ref_angles, XMLArray::Array<int>  in); 

}


#endif