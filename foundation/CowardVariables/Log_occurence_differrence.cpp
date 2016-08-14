#pragma warning( disable : 4786 )

#include "Log_occurence_differrence.h"
#include "CowardVariables.h"

//#include "Base_distance_to_claster.h"

#include "../WiseVariables/get_property_ID_by_property_name.h"

#include "../Frequency_extrapolation/Frequency_extrapolation.h"

using namespace std;

#include <sstream>
#include <cassert> 
#include <cmath>
#include <iostream>
#include <fstream>

extern ofstream log_stream;
 
Log_occurence_differrence::Log_occurence_differrence(
	CowardVariables & cowa_store,
	const string	& task_string
	) : Property(cowa_store)
{


// Log_occurence_differrence  0	DUMB	W5_01234_dx2 0     1 
	istringstream ist( task_string );

	ist >> name_ ;
	assert ( name_ == "Log_occurence_differrence" ) ;

	int tmp_int; 
	ist >> tmp_int;
	//is_subsidiary_ = ( tmp_int == 1 ) ? true : false;

	ist >> variable_name_in_list_;

	ist >> frequency_map_name_;
	ist >> 	claster_index_ ;

	if ( ! (ist >>  power_)  ) 
		power_ = 1;

	
	if  (cowa_store_.frequency_name_to_pull_index_.end() != cowa_store_.frequency_name_to_pull_index_.find ( frequency_map_name_)  ) // нашел значит
		index_in_frequency_pull_ = cowa_store_.frequency_name_to_pull_index_[frequency_map_name_];
	else 
	{
		cout << frequency_map_name_ ;
		cout << "Log_occurence_differrence() ERROR: can't associate Frequency_extrapolation object for name " << endl;
	//	log  << "Log_occurence_differrence() ERROR: can't associate Frequency_extrapolation object for name " <<  frequency_map_name_ << endl;
		exit (-1);
	}
}

Property* Log_occurence_differrence::
clone ( 	const string			& task_string)  const
{
		return new Log_occurence_differrence(cowa_store_, task_string);
}
	
double    Log_occurence_differrence::
calc_value ( const int   position_in_chain )  
{
	vector <int> const & occurence = cowa_store_.pull_fcc_[index_in_frequency_pull_ ].respect_occurence ();
	int total_sample_size  = cowa_store_.pull_fcc_[index_in_frequency_pull_ ].total_sample_size();


	vector < double>	const & current_array	= cowa_store_.pull_fcc_[index_in_frequency_pull_ ].distance_to_clusters_sum()[position_in_chain];
	vector < double>    const & global_array    = cowa_store_.pull_fcc_[index_in_frequency_pull_ ].tot_distance_to_clusters_sum();

	double ocu = (double )  occurence [position_in_chain];


	double sum;
	if ( occurence [position_in_chain] ) 
	{
		double test1 = global_array [claster_index_]/total_sample_size   ;
		double test2 = current_array [claster_index_]/occurence [position_in_chain];

		sum = log ( 1 + log ( 1 + ocu ) ) * ( global_array [claster_index_]/total_sample_size   - current_array [claster_index_]/occurence [position_in_chain]) ;
	}
	else 
		sum = 0;

	//sophisticated_variables [position_in_chain][var_set_cursor] = pow ( sum,power_ );

	return pow(sum, power_);
}
