#ifndef SINGLE_CLUSTER_RECORD_H 
#define SINGLE_CLUSTER_RECORD_H 

#include <vector>
#include <string>

#include "../Pair_int_double.h"

using namespace std;

class Single_cluster_record
{
public:

	Single_cluster_record () {}

	Single_cluster_record (
		const int		in_data_base_index,
		const int		neighbour_number,
		const double	distance_sum_to_origin,
		const double	distance_sum_to_origin_squared,
		const vector < Pair_int_double >  & pair_neighbour_index_distance);

	int		in_data_base_index				() const  { return in_data_base_index_; } 
	int		neighbour_number				() const  { return neighbour_number_;   }
	double	distance_sum_to_origin			() const  { return distance_sum_to_origin_; }
	double	distance_sum_to_origin_squared	() const  { return distance_sum_to_origin_squared_; }
	vector < Pair_int_double >	pair_neighbour_index_distance () const  { return pair_neighbour_index_distance_;  }

	friend bool	 operator >  (  const Single_cluster_record & v1,const Single_cluster_record & v2 ) {
		return ( v1.neighbour_number () > v2.neighbour_number()) ;									}	


	friend bool	 operator <  (  const Single_cluster_record & v1,const Single_cluster_record & v2 ) {
		return ( v1.neighbour_number () < v2.neighbour_number()) ;									}


	friend bool	 operator ==  (  const Single_cluster_record & v1,const Single_cluster_record & v2 ) { 
		return ( v1.neighbour_number () == v1.neighbour_number()) ;									 }

private:
	int		in_data_base_index_;
	int		neighbour_number_;
	double	distance_sum_to_origin_;
	double	distance_sum_to_origin_squared_;

	vector < Pair_int_double >  pair_neighbour_index_distance_;

};
#endif