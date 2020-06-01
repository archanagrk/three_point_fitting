#ifndef __FF_QSQ_FUNCTIONS_H__
#define __FF_QSQ_FUNCTIONS_H__

#include "fitting_lib/functions.h"
#include "fitting_lib/tools.h"
#include "fitting_lib/params.h"



//**************************************************************************************
/* class derived from Function
   only difference is the addition of 
   an operator()(double qsq, ...)
   static_cast<FFqsqFunction*>(function) can be used to convert Function pointers
*/
//***************************************************************************************

class FFqsqFunction : public Function {
 public:
  virtual ~FFqsqFunction(){};

  /* still not specified */
  virtual string name() const = 0;                                                        /* you must provide a function name string */
  virtual double operator()( const Abscissa& x, const mapstringdouble& pars ) const = 0;  /* you must provide a function returing a double */

  /* the actual function */
  virtual double operator()( double qsq, const mapstringdouble& pars ) const = 0;

};




//~~~~~~~~~~~~~~~~~~~~~~~~~~~
// qsq 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* an xml control struct for factory construction */
struct ff_qsq_params{
  ff_qsq_params(){};
  ff_qsq_params(XMLReader& xml_in, const string& path);
  
  double m_p;
  double qsq_max;
  bool topt;
  double tcut;
  int n_max;

};



  //************************************************************************************************
/* fit function for VMD */
class FFqsqVMD : public FFqsqFunction{
 public:
  FFqsqVMD(double m_p );                  /* construct with explicit double */
  FFqsqVMD( ff_qsq_params p );    /* construct by factory approach    */
  
  string name() const { return the_name; };
  double operator()( const Abscissa& x, const mapstringdouble& pars ) const;

  /* the actual function that does the work */
  double operator()( double qsq, const mapstringdouble& pars ) const ;
  
 private:
  string the_name; 
  double m_p;
  void initialize(); /* builds the name & logs the param names */
};



  //************************************************************************************************
/* fit function for Gaussian */
class FFqsqGauss : public FFqsqFunction{
 public:
  FFqsqGauss(double m_p_ );          /* construct with explicit double */
  FFqsqGauss( ff_qsq_params p );    /* construct by factory approach    */
  
  string name() const { return the_name; };
  double operator()( const Abscissa& x, const mapstringdouble& pars ) const;

  /* the actual function that does the work */
  double operator()( double qsq, const mapstringdouble& pars ) const ;
  
 private:
  string the_name; 
  double m_p;
  void initialize(); /* builds the name & logs the param names */
};



  //************************************************************************************************
/* fit function for Z Expansion */
class FFqsqZExp : public FFqsqFunction{
 public:
  FFqsqZExp(double m_p_, double qsq_max_, bool topt_, double tcut_, int n_max_);          /* construct with explicit double */
  FFqsqZExp( ff_qsq_params p );    /* construct by factory approach    */
  
  string name() const { return the_name; };
  double operator()( const Abscissa& x, const mapstringdouble& pars ) const;

  int get_nmax() const { return n_max; };
  double get_qsqmax() const { return qsq_max; };

  /* the actual function that does the work */
  double operator()( double qsq, const mapstringdouble& pars ) const ;
  
 private:
  string the_name; 
  double m_p, qsq_max, tcut;
  bool topt;
  int n_max;
  void initialize(); /* builds the name & logs the param names */
};




#endif
