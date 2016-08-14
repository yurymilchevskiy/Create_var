#include <cmath>


double exponent_factor ( const int point_number, const int shift ) 
{
	double argument =  6* fabs ( double ( shift )   ) / point_number;
	double  value = exp ( - argument*argument/2    ) ;

	return value ;
}