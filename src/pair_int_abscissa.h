#ifndef __PAIR_INT_ABSCISSA_H__
#define __PAIR_INT_ABSCISSA_H__

#include "fitting_lib/abscissa.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// simple abscissa derived class example
// keying y-values by a pair of x-values
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PairIntAbscissa : public Abscissa
{
 public:
  PairIntAbscissa(std::pair<int, int> x_){ x = x_; }
  
  /* members required by the Abscissa base class */
  string abscissa_type()    const { return "pair_int_abscissa"; }
  void   print(ostream& os) const { os << x.first << "," << x.second; }
  
  inline bool is_less_than(const Abscissa& rhs) const; 
  inline bool is_equal_to( const Abscissa& rhs) const;
  
  /* get at the value, with const protection */
  std::pair<int,int> get_x() const { return std::make_pair(x.first, x.second); }
  
 private:
  std::pair<int,int> x; 
};


/* implementations */
bool PairIntAbscissa::is_less_than(const Abscissa& rhs) const {  

  if( rhs.abscissa_type() != "pair_int_abscissa" )
    { cerr << "Abscissa:: can't compare rhs=" << rhs.abscissa_type() << " with lhs=pair_int_abscissa, exiting" << endl; exit(1); }

  /* this runtime checking may be a waste of flops */
  
  return ((x.first == (static_cast<const PairIntAbscissa&>(rhs)).get_x().first) && x.second < (static_cast<const PairIntAbscissa&>(rhs)).get_x().second) \
          || (x.first < (static_cast<const PairIntAbscissa&>(rhs)).get_x().first);
}


bool PairIntAbscissa::is_equal_to(const Abscissa& rhs) const {

  if( rhs.abscissa_type() != "pair_int_abscissa" )
    { cerr << "Abscissa:: can't compare rhs=" << rhs.abscissa_type() << " with lhs=pair_int_abscissa, exiting" << endl; exit(1); }
  /* this runtime checking may be a waste of flops */
  
  return x == (static_cast<const PairIntAbscissa&>(rhs)).get_x();
}



#endif
