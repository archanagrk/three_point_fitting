#ifndef __Z_TIMESLICE_FUNCTIONS_FACTORY_H__
#define __Z_TIMESLICE_FUNCTIONS_FACTORY_H__

/* needs the definition of the factory fromthe library */
#include "fitting_lib/functions_factory.h"


/* locally specified functions */
#include "z_timeslice_functions.h"


namespace ZFunctionEnv
{ 
  bool registerAll();
}



#endif
