#include "three_pt_fit_quality.h"

//*****************************************************************************************************************

bool cmp(double a, double b){

  if(int( round(a*100) - round(b*100) ) < 9){return true;}
  else if(int( round(a*100) + round(b*100) )  < 9){return true;}
  else{return false;} 
}

const double PI = (atan(double(1)) * double(4.0));


/* an xml control struct for factory construction */
gen3pt_params:: gen3pt_params(XMLReader& xml_in, const string& path)
{
  try {
    XMLReader paramtop(xml_in, path); 
    read(paramtop, "power", power);
    read(paramtop, "powerT", power_time);
    read(paramtop, "powerExp", power_exp);
    read(paramtop, "expFac", multiply_exp);
    read(paramtop, "dev", deviation);
  }
  catch(const std::string& e) 
    {  std::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;  exit(1);   }
};

//*****************************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Generic error estimator
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//*****************************************************************************************************************

double Generic3pt::operator()( const AvgFit& fit ) const{

  chisq_dof xx = fit.get_chisq();
  double y = double( xx.n_data - xx.n_reset_data_sv - xx.n_pars + xx.n_fixed_pars ) /  xx.chisq;

//*****************************************************************************************************************
  /* raise time to the specified power */

  y = pow(y, power);

  vector<bool> active = fit.get_active_data();
  int active_t = 0;

  for(int i = 0; i < active.size(); i++){ if(active[i]){active_t++;} else{continue;} }

  double ratio_timeslice = double(active_t)/double(active.size());

  y *= pow(ratio_timeslice, power_time);

//*****************************************************************************************************************
  /* remove fits with huge correlation between pars */

  vector<string> pars = fit.get_function()->get_par_list();
  map<pair<string,string>,double> par_corr = fit.get_result().par_corr;
  
  for(auto it = par_corr.begin(); it != par_corr.end(); it++){
    if((it->second > 0.8 || it->second < -0.7) && (it->first.first != it->first.second)){return 0;}
  }

//*****************************************************************************************************************
  /* If the function is too different from the data */

  /* Number of Dt */
  string name = fit.get_fit_name();
  size_t n = std::count(name.begin(), name.end(), 'x');
  n++;

  /* Area under the data through quadrature */

  const vector<double>& mean_y_vals = fit.get_data().get_all_y_mean();
  const vector<double>& mean_y_err = fit.get_data().get_all_y_err();
  vector<Abscissa*> x_vals = fit.get_data().get_all_x();

  vector<double> area(n, 0.0);
  int dt_num = 0;
  for(int i = 1; i < mean_y_vals.size()-2; i++ )
  {
    PairIntAbscissa* x1 = static_cast<PairIntAbscissa*>(x_vals.at(i));
    PairIntAbscissa* x2 = static_cast<PairIntAbscissa*>(x_vals.at(i+2));

    if(x1->get_x().first == x2->get_x().first){
      double y1 = mean_y_vals.at(i);
      double y2 = mean_y_vals.at(i+1);

      area.at(dt_num) += 0.5*(abs(y1) + abs(y2));
    }
    else{dt_num+3;}
  }

  /* Area under the function */

  vector<double> area_fn(n, 0.0);
  vector<int> Dts;
  ThreePointtimesliceCorrFunction* fn =  static_cast<ThreePointtimesliceCorrFunction*>(fit.get_function());

  /* Pars as mapstringdouble for the function input */
  mapstringdouble params;
  map<string, param_value>  par_vals = (fit.get_result()).par_values;

  for(auto par_val = par_vals.begin(); par_val != par_vals.end(); par_val++ ){
    params.insert(make_pair(par_val->first , par_val->second.value) );
  }

  //The constant which every fit form has
  param_value F = ( par_vals.find("F") )->second;

    dt_num = 0;
    for(int i = 1; i < mean_y_vals.size()-2; i++ )
    {

      PairIntAbscissa* x1 = static_cast<PairIntAbscissa*>(x_vals.at(i));
      PairIntAbscissa* x2 = static_cast<PairIntAbscissa*>(x_vals.at(i+2));

      if(x1->get_x().first == x2->get_x().first){
        double y1 = (*fn)(*x1,params);
        double y2 = (*fn)(*x2,params);

        area_fn.at(dt_num) += 0.5*(abs(y1) + abs(y2));
      }
      else{dt_num+3;}
    }

    /* Error in data fits */

    vector<double> dev_data(n, 0.0);
    vector<double> err(n, 0.0);
    vector<double> dev_F(n,0.0);


    dt_num = 0;
    for(int i = 1; i < mean_y_vals.size()-2; i++ )
    {

      PairIntAbscissa* x1 = static_cast<PairIntAbscissa*>(x_vals.at(i));
      PairIntAbscissa* x2 = static_cast<PairIntAbscissa*>(x_vals.at(i+2));

      if(x1->get_x().first == x2->get_x().first){

        double y1_fn = (*fn)(*x1,params);
        double y1 = mean_y_vals.at(i);

        double e1 = mean_y_err.at(i);

        double e1_fn = 0;
        for(auto p = par_vals.begin(); p != par_vals.end(); p++){ 

          if(p->first.find("Ei") != std::string::npos){
            param_value Fi = ( par_vals.find("Fi_Dt"+p->first.substr(p->first.find("t") + 1) ) )->second;
            e1_fn += (Fi.error-( p->second.error * x1->get_x().second * Fi.value ))* exp(-p->second.value * x1->get_x().second);
            }

          else if(p->first.find("Ef") != std::string::npos){
            param_value Ff = ( par_vals.find("Ff_Dt"+p->first.substr(p->first.find("t") + 1) ) )->second;
            e1_fn += (Ff.error-( p->second.error * x1->get_x().second * Ff.value ))* exp(-p->second.value * x1->get_x().second);
            }

          //else if(p->first == "F"){e1_fn += p->second.error; dev_F.at(dt_num) = abs( (p->second.value - y1)/(p->second.error - e1) );}
          else if(p->first == "F"){e1_fn += p->second.error; dev_F.at(dt_num) += abs( y1/(p->second.value - y1) );}

          else{continue;}
          }


        dev_data.at(dt_num) += abs( y1/(y1_fn - y1) );
        //err.at(dt_num) += abs(e1_fn - e1);
        err.at(dt_num) += abs(e1/e1_fn);

      }
      else{dt_num+3;}
    }




  /* 
  
  Multiply by |data|/|Model - Data| - to make sure the points in the model are not too far away from data points
  Multiply by |data|/|F - Data| - to make sure that F is not too far away from the data

  */
  
  for(int i = 0; i < area_fn.size() ; i++)
  {

    std::size_t found  = name.find("x");
    string this_name = name.substr(0,found);

    int Dt;
    string dt_i = this_name.substr(this_name.find("Dt") + 2);
    std::size_t dt_f = dt_i.find("_");
    std::istringstream str2int(  dt_i.substr( 0, dt_f ) );
    str2int >> Dt;  
    Dts.push_back(Dt);


    // If the the constant is too away from the data then return 0

    // if(abs(abs(F.value*Dt) - area.at(i)) > (deviation*area.at(i))){
      //cout << i << ":arear:" << abs(abs(area_fn.at(i)) - area.at(i)) << ":" << area.at(i) << endl;
      //cout << "warning: the constant in the fit is more than " << deviation << " sigma away"<< endl; 
      //y *= area.at(i)/abs(abs(F.value*Dt) - area.at(i));
      // }

    //if(!(10 * abs(abs(F.value*Dt) - area.at(i)) > abs(abs(area_fn.at(i)) - area.at(i)))){
      //cout << i << ":arear:" << abs(abs(area_fn.at(i)) - area.at(i)) << ":" << area.at(i) << endl;
      //cout << "warning: the constant is a better fit than exp"<< endl; 
      //y *= abs(abs(area_fn.at(i)) - area.at(i))/abs(abs(F.value*Dt) - area.at(i));
      //}


    // if the areas deviate a lot 
    // if(abs(abs(area_fn.at(i)) - area.at(i)) > (deviation*area.at(i))){
    //   //cout << i << ":area:" << abs(abs(area_fn.at(i)) - area.at(i)) << ":" << area.at(i) << endl;
    //   cout << "warning: the src exp in the fit is more than " << deviation << " away"<< endl; 
      //y *= area.at(i) / abs(abs(area_fn.at(i)) - area.at(i));
      // }




    if(found != std::string::npos){name = name.substr(found+1);}

    // area_fn.at(i) = abs(area_fn.at(i)); 

    y *= pow(dev_data.at(i),2) * dev_F.at(i) * err.at(i);


  }  


  //*****************************************************************************************************************
  /* and optionally multiply by the mass parameter */

  vector<bool> src(n, false); int i = 0;
  vector<bool> snk(n, false); int f = 0;

  vector<double> slope_i(n,0.0), error_i(n,0.0);
  vector<double> slope_f(n,0.0), error_f(n,0.0);

  for(int p = 0; p < pars.size(); p++){

    if( ( pars.at(p).find("Ei") != std::string::npos ) && multiply_exp ){

      int Dt;
      std::istringstream str2int(pars.at(p).substr(pars.at(p).find("t") + 1));
      str2int >> Dt;

      param_value E_i = ( par_vals.find(pars.at(p)) )->second;
      param_value F_i = ( par_vals.find("Fi_Dt"+pars.at(p).substr(pars.at(p).find("t") + 1) ) )->second;
      
      //y *= pow( sin(abs(F_i.value/(E_i.value*F.value*Dt) *(exp(E_i.value*Dt) - 1))*2/PI), power_exp);
      y *= pow( (E_i.value*F_i.value)/ (E_i.error*F_i.error), power_exp);
      src.at(i) = true; i++;

    }


    if( ( pars.at(p).find("Ef") != std::string::npos ) && multiply_exp ){

      int Dt;
      std::istringstream str2int(pars.at(p).substr(pars.at(p).find("t") + 1));
      str2int >> Dt;

      param_value E_f = ( par_vals.find(pars.at(p)) )->second;
      param_value F_f = ( par_vals.find("Ff_Dt"+pars.at(p).substr(pars.at(p).find("t") + 1) ) )->second;
      
      //y *= pow( sin(abs(F_f.value/(E_f.value*F.value*Dt)*(exp(E_f.value*Dt) - 1))*2/PI), power_exp);
      y *= pow( (F_f.value*E_f.value)/ (F_f.error*E_f.error) , power_exp);
      snk.at(f) = true; f++;

    }
  }

  if(multiply_exp){

    int slope = 0;
    int start_point, midpoint;
    int dt_num;

    for(dt_num = 0; dt_num < mean_y_vals.size()-1; dt_num++){

      PairIntAbscissa* x1 = static_cast<PairIntAbscissa*>(x_vals.at(dt_num));
      PairIntAbscissa* x2 = static_cast<PairIntAbscissa*>(x_vals.at(dt_num+1));

      if(x1->get_x().first == x2->get_x().first ){continue;}

      else{

        start_point = 1;
        for(int pt = 0; pt < slope; pt++){start_point += Dts.at(pt) + 1; }
        midpoint = start_point + Dts.at(slope)/2 -1;

        error_i.at(slope) = sqrt(pow(mean_y_err.at(start_point),2) + pow(mean_y_err.at(midpoint -1),2));
        error_f.at(slope) = sqrt(pow(mean_y_err.at(dt_num -1),2) + pow(mean_y_err.at(midpoint -1),2));
        slope_i.at(slope) = abs(mean_y_vals.at(start_point) - mean_y_vals.at(midpoint)/(start_point - midpoint ));
        slope_f.at(slope) = abs(mean_y_vals.at(dt_num -1) - mean_y_vals.at(midpoint)/(dt_num -1 - midpoint )); 
        slope++;  
      }

    }

    start_point = 1;
    for(int pt = 0; pt < slope; pt++){start_point += Dts.at(pt) +1; }
    midpoint = start_point + Dts.at(slope)/2 -1;

    error_i.at(slope) = sqrt(pow(mean_y_err.at(start_point),2) + pow(mean_y_err.at(midpoint),2));
    error_f.at(slope) = sqrt(pow(mean_y_err.at(dt_num -1),2) + pow(mean_y_err.at(midpoint),2));
    slope_i.at(slope) = abs(mean_y_vals.at(start_point) - mean_y_vals.at(midpoint)/(start_point - midpoint ));
    slope_f.at(slope) = abs(mean_y_vals.at(dt_num -1) - mean_y_vals.at(midpoint)/(dt_num -1 - midpoint )); 


    for(int p = 0; p < n; p++){
      if(!src.at(p)){y *= pow(slope_i.at(p)/error_i.at(p),power_exp);}
      if(!snk.at(p)){y *= pow(slope_f.at(p)/error_f.at(p),power_exp);}
    }
    
  }
  return y;
}


//*****************************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// QN Gen error estimator
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//*****************************************************************************************************************
double QNGen::operator()( const AvgFit& fit ) const {
  
  chisq_dof xx = fit.get_chisq();
  
  int n_dof = xx.n_data - xx.n_reset_data_sv - xx.n_pars + xx.n_fixed_pars;
  
  NR::DP c = xx.chisq / 2.0;
  NR::DP n = double(n_dof) / 2.0;

  double y;

  if(n > 0.0){ y = double(n_dof)*double(NR::gammq(n,c));}
  else{        y = -1.0; }   

    
  /* raise this to the specified power */
  y = pow(y, power);

  vector<bool> active = fit.get_active_data();

  int active_t = 0;

  for(int i = 0; i < active.size(); i++){ if(active[i]){active_t++;} else{continue;} }


  double ratio_timeslice = double(active_t)/double(active.size());

  map<string, param_value>  par_vals = (fit.get_result()).par_values;

  param_value F = ( par_vals.find("F") )->second;

  y *= pow(ratio_timeslice, power_time)*pow(F.value/F.error, power_exp);


  /* and optionally multiply by the mass parameter */

  if( ( fit.get_fit_name().find("csrc") != std::string::npos ) || ( fit.get_fit_name().find("c2") != std::string::npos ) && multiply_exp ){
    param_value mass_i = ( par_vals.find("mi") )->second;
    param_value F_i = ( par_vals.find("Fi") )->second;
    
    y *= pow(mass_i.value/mass_i.error, power_exp)*pow(F_i.value/F_i.error, power_exp);

  }
  
    if( ( fit.get_fit_name().find("csnk") != std::string::npos ) || ( fit.get_fit_name().find("c2") != std::string::npos ) && multiply_exp ){
    param_value mass_f = ( par_vals.find("mf") )->second;
    param_value F_f = ( par_vals.find("Ff") )->second;
    
    y *= pow(mass_f.value/mass_f.error, power_exp)*pow(F_f.value/F_f.error, power_exp);

  }
  
  return y; 
  
   
};
