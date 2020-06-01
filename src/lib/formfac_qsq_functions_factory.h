#ifndef __FF_QSQ_FUNCTIONS_FACTORY_H__
#define __FF_QSQ_FUNCTIONS_FACTORY_H__

/* needs the definition of the factory fromthe library */
#include "fitting_lib/functions_factory.h"


/* locally specified functions */
#include "formfac_qsq_functions.h"


namespace FFunctionEnv
{ 
  bool registerAll();
}



#endif
