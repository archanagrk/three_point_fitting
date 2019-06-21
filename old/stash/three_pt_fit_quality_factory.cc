#include "three_pt_fit_quality_factory.h"


namespace{ 
  bool registered = false;
}

namespace{
  FitQuality*        createGeneric3pt( XMLReader& xml, const string& path ){ return new Generic3pt( gen3pt_params(xml, path) ); }  
}


namespace ThreePtFitQualityEnv
{ 
  
  bool registerAll(){

    /* register the functions in the library */
    bool success = FitQualityEnv::registerAll();

    /* add the locally defined functions */
    if(!registered){
      
      /* N.B. still TheFunctionFactory as we're adding to the factory defined in the library */
      success &= TheFitQualityFactory::Instance().registerObject( "generic",
								  createGeneric3pt);
      registered = true;
    }
    
    return success;
  };


}
