#pragma warning( disable : 4786 )

#include "fill_up_appropriate_crude_index_subset.h"
#include < cassert>

void fill_up_appropriate_crude_index_subset (
	vector < vector < int > >   & index_of_degeneration,
	const int position,
	const int lb,
	const int rb, 
	const vector < int > & crude_index_array,
	const int virtual_residue_index, 
	vector < int > & index_subset )
{
	int size = crude_index_array.size();
	
	int window_size = rb -  lb + 1;
	int start = position +  lb ;
	int end   = position +  rb + 1 ;

	assert ( index_subset.size() ==  window_size );

	int kk=0;
	for (int ii=start ; ii<end; ii++,kk++ )
		index_subset[kk] = ( ii >= 0 && ii < size ) ? index_of_degeneration [kk] [ crude_index_array[ii] ] : index_of_degeneration [kk] [ virtual_residue_index ];

}
