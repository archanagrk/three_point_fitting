#ifndef __MAKE_ELAB_H__
#define __MAKE_ELAB_H__


//adat
#include "io/adat_io.h"
#include "io/adat_xml_group_reader.h"
#include "ensem/ensem.h"

//std lib
#include "math.h"
#include <stdio.h>
#include <vector>
#include <iostream>

//semble
#include "semble/semble_meta.h"

//namespaces
using namespace std;
using namespace ADAT;
using namespace ADATXML;
using namespace ENSEM;

//functions

EnsemReal make_elab(EnsemReal ElabO,double momsq);

#endif