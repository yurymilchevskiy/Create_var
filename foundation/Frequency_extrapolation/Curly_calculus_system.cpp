#include "Curly_calculus_system.h"

Curly_calculus_system::
Curly_calculus_system ( const vector < int > & base ) :
  base_				(base),
  number_of_elements_	(1)
{
	  
	  for (int ii=0;ii<base_.size();ii++ )
		number_of_elements_ *= base_ [ii];
}

	
int	 Curly_calculus_system::
get_cursor_by_array ( const vector <int> array)
{
	int temp = number_of_elements_ ;

	int cursor = 0;

	for (int ii=0;ii<base_.size();ii++ )
	{
		temp /=   base_ [ii];
		cursor += temp * array [ii];
	}

	return cursor ;
}


vector < int >	Curly_calculus_system::
get_array_by_cursor	( int cursor ) 
{
	vector <int> array;
	array.resize( base_.size() ) ;

	int divisor = number_of_elements_ ;

	for (int ii=0;ii<base_.size();ii++ )
	{
		divisor /= base_ [ii];
		array [ii]	= cursor / divisor ;
		cursor		= cursor % divisor ;
	}
	return array ;
}
