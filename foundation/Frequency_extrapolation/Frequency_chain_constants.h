#ifndef FREQUENCY_CHAIN_CONSTANTS_H
#define FREQUENCY_CHAIN_CONSTANTS_H

#include <vector>
#include <fstream>

using namespace std;

class Frequency_extrapolation ;
class Frequency_extrapolation_ZIP ;

class Frequency_chain_constants
{
public:

	Frequency_chain_constants( 
		vector < int > & processed_index,
		ifstream & data_stream,
		Frequency_extrapolation * fr_ext );
	
	~Frequency_chain_constants() {} ;

	vector <vector <double> > const & distance_to_clusters_sum					() const { return distance_to_clusters_sum_					;}			 
	vector <vector <double> > const & squared_distance_to_clusters_sum			() const { return squared_distance_to_clusters_sum_			;}
	vector <vector <double> > const & inverse_distance_to_clusters_sum			() const { return inverse_distance_to_clusters_sum_			;}
	vector <vector <double> > const & inverse_squared_distance_to_clusters_sum	() const { return inverse_squared_distance_to_clusters_sum_	;}

	vector <double>  const & tot_distance_to_clusters_sum					() const { return tot_distance_to_clusters_sum_				;}
	vector <double>  const & tot_squared_distance_to_clusters_sum			() const { return tot_squared_distance_to_clusters_sum_		;}
	vector <double>  const & tot_inverse_distance_to_clusters_sum			() const { return tot_inverse_distance_to_clusters_sum_		;}
	vector <double>  const & tot_inverse_squared_distance_to_clusters_sum	() const { return tot_inverse_squared_distance_to_clusters_sum_;}
	
	vector < int >	const & respect_occurence() const { return respect_occurence_	;}			 
	int				 total_sample_size() const { return total_sample_size_ ;}			 

private:

	vector <vector <double> > distance_to_clusters_sum_					;
	vector <vector <double> > squared_distance_to_clusters_sum_			;
	vector <vector <double> > inverse_distance_to_clusters_sum_			;
	vector <vector <double> > inverse_squared_distance_to_clusters_sum_	;

	vector <double>  tot_distance_to_clusters_sum_					;
	vector <double>  tot_squared_distance_to_clusters_sum_			;
	vector <double>  tot_inverse_distance_to_clusters_sum_			;
	vector <double>  tot_inverse_squared_distance_to_clusters_sum_	;
	
	vector < int >	 respect_occurence_	;
	int				 total_sample_size_ ;
};


#endif


