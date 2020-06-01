#include "z_timeslice_functions_factory.h"


namespace{ 
  bool registered = false;
}

namespace{
  Function* createZtimesliceNExp( XMLReader& xml, const string& path ){ return new ZtimesliceNExp( z_timeslice_nexp_params(xml, path) ); } 
}


namespace ZFunctionEnv
{ 
  
  bool registerAll(){

    /* register the functions in the library */
    bool success = FunctionEnv::registerAll();

    /* add the locally defined functions */
    if(!registered){
      
      /* N.B. still TheFunctionFactory as we're adding to the factory defined in the library */
      success &= TheFunctionFactory::Instance().registerObject( "z_timeslice_nexp",
								createZtimesliceNExp);  
      registered = true;
    }
    
    return success;
  };


}
