#pragma warning( disable : 4786 )

#ifndef PAIR_INT_DOUBLE_H
#define PAIR_INT_DOUBLE_H

#include <algorithm>

using namespace std;

class Pair_int_double
{
public:
	Pair_int_double () {}
	Pair_int_double ( const int index, const double value ); 
	int		index () const { return index_; } 
	double	value () const { return value_; } 

//   friend bool         operator < ( const Text & l, const Text & r );

	friend bool     operator >  (  const Pair_int_double & v1,const Pair_int_double & v2 ) 
	{ return ( v1.value() > v2.value() ) ; }
		;
	friend bool     operator <  (  const Pair_int_double & v1,const Pair_int_double & v2 )   
	{ return ( v1.value() < v2.value() ) ; }
	friend bool     operator == (  const Pair_int_double & v1,const Pair_int_double & v2 )   
	{ return ( v1.value() == v2.value() ) ; }


private:
	int		index_;
	double	value_;

};

#endif
