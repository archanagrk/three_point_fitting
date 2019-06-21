#include "three_pt_fit_quality.h"

//*****************************************************************************************************************


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Generic 3pt error estimator
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* an xml control struct for factory construction */
gen3pt_params:: gen3pt_params(XMLReader& xml_in, const string& path)
{
  try {
    XMLReader paramtop(xml_in, path); 
    read(paramtop, "power", power);
    read(paramtop, "power_mass", power_mass);
    read(paramtop, "power_F", power_F);
    read(paramtop, "power_time", power_time);
    read(paramtop, "multiply_mass", multiply_mass);
  }
  catch(const std::string& e) 
    {  std::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;  exit(1);   }
};

double Generic3pt::operator()( const AvgFit& fit ) const{

  chisq_dof xx = fit.get_chisq();
  double y = double( xx.n_data - xx.n_reset_data_sv - xx.n_pars + xx.n_fixed_pars ) /  xx.chisq;

  /* raise this to the specified power */
  y = pow(y, power);

  vector<bool> active = fit.get_active_data();

  int active_t = 0;

  for(int i = 0; i < active.size(); i++){ if(active[i]){active_t++;} else{continue;} }


  double ratio_timeslice = double(active_t)/double(active.size());

  param_value F = ( (fit.get_result()).par_values.find("F") )->second;

  y *= pow(ratio_timeslice, power_time)*pow(F.value/F.error, power_F);


  /* and optionally multiply by the mass parameter */

  if( ( fit.get_fit_name().find("csrc") != std::string::npos ) || ( fit.get_fit_name().find("c2") != std::string::npos ) && multiply_mass ){
    param_value mass_i = ( (fit.get_result()).par_values.find("mi") )->second;
    param_value F_i = ( (fit.get_result()).par_values.find("Fi") )->second;
    
    y *= pow(mass_i.value/mass_i.error, power_mass)*pow(F_i.value/F_i.error, power_F);
  }
  
    if( ( fit.get_fit_name().find("csnk") != std::string::npos ) || ( fit.get_fit_name().find("c2") != std::string::npos ) && multiply_mass ){
    param_value mass_f = ( (fit.get_result()).par_values.find("mf") )->second;
    param_value F_f = ( (fit.get_result()).par_values.find("Ff") )->second;
    
    y *= pow(mass_f.value/mass_f.error, power_mass)*pow(F_f.value/F_f.error, power_F);
  }
  return y;
}


//*****************************************************************************************************************
