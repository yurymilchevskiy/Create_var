#ifndef MEAN_COMPARE_PROBABILITY_INV_H
#define MEAN_COMPARE_PROBABILITY_INV_H

#ifndef PROPERTY_H
#include "Property.h"
#endif

#include <string>
#include <vector>

using namespace std;

#include "../Frequency_extrapolation/Frequency_chain_constants.h"

class Frequency_extrapolation;

class  Mean_compare_probability_inv: public Property
{
public:
	Mean_compare_probability_inv(CowardVariables & cowa_store) :
		Property(cowa_store)
	{}

	explicit Mean_compare_probability_inv(
		CowardVariables & cowa_store,
		const string	& task_string);

	Property* clone(
		const string & task_string) const;

	double    calc_value(const int   position_in_chain);


protected:

	string variable_name_in_list_;

	string frequency_map_name_;
	int index_in_frequency_pull_ ;
	
	int claster_index_ ;

	int left_border_ ;
	int right_border_ ;

	double power_;

	Mean_compare_probability_inv(const Mean_compare_probability_inv&);
	Mean_compare_probability_inv& operator = (const Mean_compare_probability_inv&);
};

#endif