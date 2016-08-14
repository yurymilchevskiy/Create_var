#ifndef STUDENT_EMASCULATE_H
#define STUDENT_EMASCULATE_H

#ifndef PROPERTY_H
#include "Property.h"
#endif

#include <string>
#include <vector>

using namespace std;

#include "../Frequency_extrapolation/Frequency_chain_constants.h"

class Frequency_extrapolation;

class  Student_emasculate: public Property
{
public:
	Student_emasculate(CowardVariables & cowa_store) :
		Property(cowa_store)
	{}

	explicit Student_emasculate(
		CowardVariables & cowa_store,
		const string	& task_string);

	//~Student_emasculate() {}

	Property* clone(
		const string & task_string) const;

	double    calc_value(const int   position_in_chain);


	void    calc_value (
				const int   position_in_chain,  
					  int	var_set_cursor, 
					vector < Frequency_chain_constants > & pull_fcc,
					vector < vector < double > >   & sophisticated_variables    );

protected:

	string variable_name_in_list_;

	string frequency_map_name_;
	int index_in_frequency_pull_ ;
	
	int claster_index_ ;

	int left_border_ ;
	int right_border_ ;

	double power_;

	Student_emasculate(const Student_emasculate&);
	Student_emasculate& operator = (const Student_emasculate&);
};

#endif