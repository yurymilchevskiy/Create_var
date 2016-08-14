#ifndef GET_PROPERTY_CONSTANT_H
#define GET_PROPERTY_CONSTANT_H

#include <string>
#include <vector>

using namespace std;

double 
	get_property_constant ( const char one_letter_aa,	const int property_ID );
double 
	get_property_constant ( const int  aa_ID,			const int property_ID );

vector <double> 
	get_property_constant ( const string		 & seq_chain,  const int property_ID );
vector <double> 
	get_property_constant ( const vector < int > & array_aa_ID,	const int property_ID );


#endif