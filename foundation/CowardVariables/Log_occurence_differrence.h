#ifndef LOG_OCCURENCE_DIFFERRENCE_H
#define LOG_OCCURENCE_DIFFERRENCE_H

#ifndef PROPERTY_H
#include "Property.h"
#endif

#include <string>
#include <vector>

using namespace std;

#include "../Frequency_extrapolation/Frequency_chain_constants.h"

class Frequency_extrapolation;

class  Log_occurence_differrence: public Property
{
public:
	Log_occurence_differrence(CowardVariables & cowa_store) :
		Property(cowa_store)
	{}

	explicit Log_occurence_differrence(
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

	double power_;

	Log_occurence_differrence(const Log_occurence_differrence&);
	Log_occurence_differrence& operator = (const Log_occurence_differrence&);
};

#endif