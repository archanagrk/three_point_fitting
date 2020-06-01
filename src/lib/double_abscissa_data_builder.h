#ifndef __DOUBLE_ABSCISSA_DATA_BUILDERS_H__
#define __DOUBLE_ABSCISSA_DATA_BUILDERS_H__

#include "fitting_lib/data.h"
#include "fitting_lib/double_abscissa.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// optional functions to build data objects
// hiding the 'new' of the Abscissa objects
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


inline void make_double_abscissa_data( const vector<double>& x_double, const vector<EnsemReal>& y_ensem,
				    Data& data ) /* you should pass in an empty Data object */ 
{  
  vector< Abscissa* > x; 
  
  for(vector<double>::const_iterator it = x_double.begin(); it != x_double.end(); it++){
    x.push_back( new DoubleAbscissa( *it ) );
  }

  data.make(x, y_ensem);
}



#endif
