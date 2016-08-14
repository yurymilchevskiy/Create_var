#include "aa_sequence_to_index_array.h"

vector < int > aa_sequence_to_index_array (const string & sequence )
{
	int sequence_size = sequence.size ();  
	vector < int > index_array; index_array.resize ( sequence_size );

	for (int ii=0;ii<sequence_size ;ii++)
		index_array[ii]  = aminoacid_to_index (   sequence [ii]   ) ;
		//index_array[ii]  = aminoacid_to_index (  ( const ) sequence [ii]   ) ;

	return index_array;
}