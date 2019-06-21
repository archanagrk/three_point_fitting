#include "three_point_timeslice_corr_functions_factory.h"


namespace{ 
  bool registered = false;
}

namespace{
 Function* createThreePointtimesliceCorrNExp( XMLReader& xml, const string& path ){ return new ThreePointtimesliceCorrNExp( three_point_timeslice_corr_nexp_params(xml, path) ); }  
}


namespace ThreePointFunctionEnv
{ 
  
  bool registerAll(){

    /* register the functions in the library */
    bool success = FunctionEnv::registerAll();

    /* add the locally defined functions */
    if(!registered){
      
      /* N.B. still TheFunctionFactory as we're adding to the factory defined in the library */
      success &= TheFunctionFactory::Instance().registerObject( "three_point_timeslice_corr_nexp",
				         createThreePointtimesliceCorrNExp);  
      registered = true;
    }
    
    return success;
  };


}
