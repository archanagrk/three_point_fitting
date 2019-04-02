#include "three_pt_fit_quality_factory.h"


namespace{ 
  bool registered = false;
}

namespace{
  FitQuality*  createGeneric3pt( XMLReader& xml, const string& path ){ return new Generic3pt( gen3pt_params(xml, path) ); } 
  FitQuality*  createQNGen( XMLReader& xml, const string& path ){ return new QNGen( gen3pt_params(xml, path) ); } 
}


namespace ThreePtFitQualityEnv
{ 
  
  bool registerAll(){

    /* register the functions in the library */
    bool success = FitQualityEnv::registerAll();

    /* add the locally defined functions */
    if(!registered){
      
      /* N.B. still TheFitQualityFactory as we're adding to the factory defined in the library */
      success &= TheFitQualityFactory::Instance().registerObject( "gen_3_pt",
								  createGeneric3pt);

      success &= TheFitQualityFactory::Instance().registerObject( "qn_gen",
								  createQNGen);
      registered = true;
    }
    
    return success;
  };


}
