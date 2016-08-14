#include "Single_cluster_record.h"

Single_cluster_record ::Single_cluster_record (
	const int		in_data_base_index,
	const int		neighbour_number,
	const double	distance_sum_to_origin,
	const double	distance_sum_to_origin_squared,
	const vector < Pair_int_double >  & pair_neighbour_index_distance) :
		in_data_base_index_				(in_data_base_index), 
		neighbour_number_				(neighbour_number), 
		distance_sum_to_origin_			(distance_sum_to_origin), 
		distance_sum_to_origin_squared_	(distance_sum_to_origin_squared),
		pair_neighbour_index_distance_  ( pair_neighbour_index_distance)

{


}



