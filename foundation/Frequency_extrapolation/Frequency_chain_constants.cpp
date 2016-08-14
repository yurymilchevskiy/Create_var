#pragma warning( disable : 4786 )

#include "Frequency_chain_constants.h"
#include "Frequency_extrapolation.h" 

#include <fstream>
#include <iostream>
using namespace std;

Frequency_chain_constants::
Frequency_chain_constants( 
	vector < int > & processed_index,
	ifstream & data_stream,
	Frequency_extrapolation * fr_ext )
{
// Attention - dangerous change
//	vector <int>	whole_occurence = 	fr_ext->get_occurence ( data_stream );

	for ( int kk=0;kk < processed_index.size(); kk++ )
	{
		vector <double>  single_distance_to_clusters_sum				 =	fr_ext->get_single_distance_to_clusters_sum					( processed_index[kk] );
		vector <double>  single_squared_distance_to_clusters_sum		 = 	fr_ext->get_single_squared_distance_to_clusters_sum			( processed_index[kk] );
		vector <double>  single_inverse_distance_to_clusters_sum		 =	fr_ext->get_single_inverse_distance_to_clusters_sum			( processed_index[kk] );
		vector <double>  single_inverse_squared_distance_to_clusters_sum =  fr_ext->get_single_inverse_squared_distance_to_clusters_sum ( processed_index[kk] );

		distance_to_clusters_sum_				 .push_back( single_distance_to_clusters_sum				);
		squared_distance_to_clusters_sum_		 .push_back( single_squared_distance_to_clusters_sum		);
		inverse_distance_to_clusters_sum_		 .push_back( single_inverse_distance_to_clusters_sum		);
		inverse_squared_distance_to_clusters_sum_.push_back( single_inverse_squared_distance_to_clusters_sum);

// Attention - dangerous change
//		respect_occurence_.push_back (whole_occurence [processed_index[kk]] );
		int curr_occurence = fr_ext->get_occurence (processed_index[kk] );
		respect_occurence_.push_back (curr_occurence );
// test only 
		if (curr_occurence < 0 )
			cout << "strange occyrence " << endl;
			


	}

// FIX ÌÎÆÅÒ ÁÛÒÜ ÏÎÍÀÄÎÁÈÒÑß
	tot_distance_to_clusters_sum_				  =	fr_ext->get_tot_distance_to_clusters_sum				 ();		 
	tot_squared_distance_to_clusters_sum_		  = fr_ext->get_tot_squared_distance_to_clusters_sum		 ();		 
	tot_inverse_distance_to_clusters_sum_		  =	fr_ext->get_tot_inverse_distance_to_clusters_sum		 ();		 
	tot_inverse_squared_distance_to_clusters_sum_ = fr_ext->get_tot_inverse_squared_distance_to_clusters_sum ();			 
	total_sample_size_							  = fr_ext->get_total_sample_size							 ();

}
