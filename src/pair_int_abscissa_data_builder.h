#ifndef __PAIR_INT_ABSCISSA_DATA_BUILDERS_H__
#define __PAIR_INT_ABSCISSA_DATA_BUILDERS_H__

#include "fitting_lib/data.h"
#include "pair_int_abscissa.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// optional functions to build data objects
// hiding the 'new' of the Abscissa objects
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


inline void make_pair_int_abscissa_data( const vector<pair<int,int>>& x_pair_int, const vector<EnsemReal>& y_ensem,
				    Data& data ) /* you should pass in an empty Data object */ 
{  
  vector< Abscissa* > x; 
  
  for(vector<pair<int,int>>::const_iterator it = x_pair_int.begin(); it != x_pair_int.end(); it++){
    x.push_back( new PairIntAbscissa( *it ) );
  }

  data.make(x, y_ensem);
}



#endif
