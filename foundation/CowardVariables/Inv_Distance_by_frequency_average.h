#ifndef INV_DISTANCE_BY_FREQUENCY_AVERAGE_H
#define INV_DISTANCE_BY_FREQUENCY_AVERAGE_H

#ifndef PROPERTY_H
#include "Property.h"
#endif

#include <string>
#include <vector>

using namespace std;

#include "../Frequency_extrapolation/Frequency_chain_constants.h"

class Frequency_extrapolation;

class  Inv_Distance_by_frequency_average: public Property
{
public:
	Inv_Distance_by_frequency_average(CowardVariables & cowa_store) :
		Property(cowa_store)
	{}

	explicit Inv_Distance_by_frequency_average(
		CowardVariables & cowa_store,
		const string	& task_string);

	Property* clone(
		const string & task_string) const;

	double    calc_value(const int   position_in_chain);

/*	
		void    calc_value ( 
				const int   position_in_chain,  
					  int	var_set_cursor, 
				const   vector < vector < vector <double > > > & Chain_Frequency_Constants,
						vector < Frequency_extrapolation * > & frequency_pull_,
				    	vector < vector < double > >   & sophisticated_variables    ); 
*/

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

	Inv_Distance_by_frequency_average(const Inv_Distance_by_frequency_average&);
	Inv_Distance_by_frequency_average& operator = (const Inv_Distance_by_frequency_average&);
};

#endif