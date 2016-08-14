#pragma warning( disable : 4786 )
	
#include "c_FourierSmoothed.h"

#include "CowardVariables.h"

#include "../WiseVariables/get_property_ID_by_property_name.h"
#include "../CommonFunc.h"
#include "../WiseVariables/exponent_factor.h" 


using namespace std;

#include <sstream>
#include <cassert>
#include <cmath>
 
c_FourierSmoothed::c_FourierSmoothed ( 
	CowardVariables & cowa_store,
	const string	& task_string
	) : Property(cowa_store)
{

// c_FourierSmoothed 	0 mbb_1 	Mass_and_backbone -1 1  2 
//c_FourierSmoothed	0	fs_NIKI_1,50	NIKI_____3	1,50	1,50	8,00
	istringstream ist( task_string );

	ist >> name_ ;
	assert ( name_ == "c_FourierSmoothed" ) ;

	int tmp_int; 
	ist >> tmp_int;
	//is_subsidiary_ = ( tmp_int == 1 ) ? true : false;

	ist >> variable_name_in_list_;

	string current_property_name;
	ist >> current_property_name;
	property_ID_ = get_property_ID_by_property_name ( current_property_name ) ;

	ist >>	period_lenfth_;
	ist >>	period_number_;

	ist >>  power_ ;


	left_border_	= - period_number_ * period_lenfth_ /2 ;
	right_border_	=	period_number_ * period_lenfth_ /2 ;

	point_number_ =  right_border_ - left_border_;

}

Property* c_FourierSmoothed::
clone ( const string & task_string ) const
{
		//return new c_FourierSmoothed(  task_string,co_task_variable_name_to_index  );
	return new c_FourierSmoothed(cowa_store_, task_string);

}

double    c_FourierSmoothed::calc_value ( 
		const int   position_in_chain ) 
{
	double Coefficient = 2*Pythagorean_Number () / period_lenfth_;

	int pre_start	= position_in_chain + left_border_;
	int pre_end		= position_in_chain + right_border_;

	int seq_len		= cowa_store_.Chain_Prime_Constants_ [ property_ID_ ].size ();

	int start	= __max ( 0,		pre_start	);
	int end		= __min ( pre_end,	seq_len		);
 
	double A_value=0, B_value=0,property_value =0;
	for ( int ii = start; ii < end; ii++ ) 
	{
		double x = ii * Coefficient ;
		double current_value = cowa_store_.Chain_Prime_Constants_ [ property_ID_ ] [ ii ] ;
		double factor = exponent_factor ( point_number_, position_in_chain - ii ) ;

		A_value += current_value * sin (x) * factor ;
		B_value += current_value * cos (x) * factor ;
	}
	double value = sqrt (A_value *A_value + B_value * B_value );

	//sophisticated_variables [position_in_chain][var_set_cursor] = pow ( value_,power_ );

	return  pow(value, power_);
}