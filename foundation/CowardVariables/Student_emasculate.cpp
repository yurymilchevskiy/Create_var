#pragma warning( disable : 4786 )

#include "Student_emasculate.h"
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
 
Student_emasculate::Student_emasculate(
	CowardVariables & cowa_store,
	const string	& task_string
	) : Property(cowa_store)
{

// Student_emasculate 	0 pl_dist_0 0  	1 


// Student_emasculate  0	DUMB	W5_01234_dx2 0     -2 2    1 
	istringstream ist( task_string );

	ist >> name_ ;
	assert ( name_ == "Student_emasculate" ) ;

	int tmp_int; 
	ist >> tmp_int;
	//is_subsidiary_ = ( tmp_int == 1 ) ? true : false;

	ist >> variable_name_in_list_;

	ist >> frequency_map_name_;
	ist >> 	claster_index_ ;

	ist >> 	left_border_ ;
	ist >> 	right_border_ ;

	if ( ! (ist >>  power_)  ) 
		power_ = 1;

	
	if  (cowa_store_.frequency_name_to_pull_index_.end() != cowa_store_.frequency_name_to_pull_index_.find ( frequency_map_name_)  ) // нашел значит
		index_in_frequency_pull_ = cowa_store_.frequency_name_to_pull_index_[frequency_map_name_];
	else 
	{
		cout << frequency_map_name_ ;
		cout << "Student_emasculate() ERROR: can't associate Frequency_extrapolation object for name " << endl;
	//	log  << "Student_emasculate() ERROR: can't associate Frequency_extrapolation object for name " <<  frequency_map_name_ << endl;
		exit (-1);
	}
}

Property* Student_emasculate::
clone ( const string & task_string) const
{
	return new Student_emasculate(cowa_store_, task_string);
}
	
double Student_emasculate::
calc_value (	const int   position_in_chain )
{
	int pre_start	= position_in_chain + left_border_;
	int pre_end		= position_in_chain + right_border_;

	int seq_len		= cowa_store_.sophisticated_variables_.size ();

	int start	= __max ( 0,		pre_start	);
	int end		= __min ( pre_end,	seq_len		);

	vector <int> const & occurence			= cowa_store_.pull_fcc_ [index_in_frequency_pull_ ].respect_occurence ();
	int total_sample_size			= cowa_store_.pull_fcc_[index_in_frequency_pull_ ].total_sample_size();
	vector < double > const & average_array = cowa_store_.pull_fcc_[index_in_frequency_pull_ ].tot_distance_to_clusters_sum();

	double total_average_value = average_array [claster_index_]/total_sample_size ;

	double sum = 0;	int local_occurence = 0; 

	for ( int ii = start; ii < end; ii++ ) 
	{
		vector < double> const & current_array = cowa_store_.pull_fcc_[index_in_frequency_pull_ ].distance_to_clusters_sum()[ii];
		if ( occurence[ii] ) 
			sum += total_average_value -  current_array [claster_index_]/occurence[ii] ;

		local_occurence += occurence[ii];

	}

	double efficient = sqrt ( (double )  local_occurence ); 
	sum *= efficient;

	int sign =  ( sum >= 0 ) ? 1 : -1;

	//sophisticated_variables [position_in_chain][var_set_cursor] = sign * pow ( fabs(sum), power_ );

	return sign * pow(fabs(sum), power_);
}
