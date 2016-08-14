#pragma warning( disable : 4786 )

#include "Dull_Sum.h"
#include "CowardVariables.h"

#include "../WiseVariables/get_property_ID_by_property_name.h"

using namespace std;

#include <sstream>
#include <cassert>
#include <cmath>
 
Dull_Sum::Dull_Sum (
	CowardVariables & cowa_store, 
	const string	& task_string
) : Property(cowa_store),
    fabs_mode_ ('a')
{

// Dull_Sum 	0 mbb_1 	Mass_and_backbone -1 1  2 
	istringstream ist( task_string );

	ist >> name_ ;
	assert ( name_ == "Dull_Sum" ) ;

	int tmp_int; 
	ist >> tmp_int;
	//is_subsidiary_ = ( tmp_int == 1 ) ? true : false;

	ist >> variable_name_in_list_;

	string current_property_name;
	ist >> current_property_name;
	property_ID_ = get_property_ID_by_property_name ( current_property_name ) ;

	ist >> 	left_border_ ;
	ist >> 	right_border_ ;

	if ( ! (ist >>  power_)  ) 
		power_ = 1;

	ist >>  fabs_mode_ ;
}

Property* Dull_Sum::
clone ( const string & task_string ) const
{
	return new Dull_Sum(cowa_store_,task_string);
}

double    Dull_Sum::calc_value ( 
		const int   position_in_chain/*,  
			  int	var_set_cursor, 
		const		vector < vector < double > >   & Chain_Prime_Constants,
					vector < vector < double > >   & sophisticated_variables */   ) 

{
	int pre_start	= position_in_chain + left_border_;
	int pre_end		= position_in_chain + right_border_;

	int seq_len		= cowa_store_.Chain_Prime_Constants_ [ property_ID_ ].size ();

	int start	= __max ( 0,		pre_start	);
	int end		= __min ( pre_end,	seq_len		);

	double sum = 0;
//	for ( int ii = start; ii < end; ii++ ) 
//		sum += Chain_Prime_Constants  [ property_ID_ ] [ ii ] ;


	if ( fabs_mode_ == 'f' || fabs_mode_ == 'F')
	{
		for ( int ii = start; ii < end; ii++ ) 
			sum += fabs ( cowa_store_.Chain_Prime_Constants_ [ property_ID_ ] [ ii ] ) ;
	}
	else
	{
		for ( int ii = start; ii < end; ii++ ) 
			sum += cowa_store_.Chain_Prime_Constants_[property_ID_][ii];
	}

	//sophisticated_variables_ [position_in_chain][var_set_cursor] = pow ( sum,power_ );

	return pow ( sum,power_ );
}