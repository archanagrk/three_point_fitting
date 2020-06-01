#include "formfac_qsq_functions_factory.h"


namespace{ 
  bool registered = false;
}

namespace{
  Function* createFFqsqVMD( XMLReader& xml, const string& path ){ return new FFqsqVMD( ff_qsq_params(xml, path) ); } 
  Function* createFFqsqGauss( XMLReader& xml, const string& path ){ return new FFqsqGauss( ff_qsq_params(xml, path) ); }   
  Function* createFFqsqZExp( XMLReader& xml, const string& path ){ return new FFqsqZExp( ff_qsq_params(xml, path) ); }   

}


namespace FFunctionEnv
{ 
  
  bool registerAll(){

    /* register the functions in the library */
    bool success = FunctionEnv::registerAll();

    /* add the locally defined functions */
    if(!registered){
      
      /* N.B. still TheFunctionFactory as we're adding to the factory defined in the library */
      success &= TheFunctionFactory::Instance().registerObject( "ff_qsq_vmd",
								createFFqsqVMD);  

      success &= TheFunctionFactory::Instance().registerObject( "ff_qsq_gaussian",
								createFFqsqGauss);  

      success &= TheFunctionFactory::Instance().registerObject( "ff_qsq_z_exp",
								createFFqsqZExp); 
                
      registered = true;
    }
    
    return success;
  };


}
