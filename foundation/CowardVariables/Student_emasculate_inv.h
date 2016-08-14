#ifndef STUDENT_EMASCULATE_INV_H
#define STUDENT_EMASCULATE_INV_H

#ifndef PROPERTY_H
#include "Property.h"
#endif

#include <string>
#include <vector>

using namespace std;

#include "../Frequency_extrapolation/Frequency_chain_constants.h"

class Frequency_extrapolation;

class  Student_emasculate_inv: public Property
{
public:
	Student_emasculate_inv(CowardVariables & cowa_store) :
		Property(cowa_store)
	{}

	explicit Student_emasculate_inv(
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

	Student_emasculate_inv(const Student_emasculate_inv&);
	Student_emasculate_inv& operator = (const Student_emasculate_inv&);
};

#endif