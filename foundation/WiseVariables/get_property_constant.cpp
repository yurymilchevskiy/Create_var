#include "get_property_constant.h"
#include "property_constant.h"

#include "../aminoacid_to_index.h" 
#include "../aa_sequence_to_index_array.h"

using namespace std;


double 	get_property_constant ( const char one_letter_aa,	const int property_ID ) 
{

	int aa_ID = aminoacid_to_index (one_letter_aa);

	/*check only */ double test =  get_property_constant ( aa_ID , property_ID  ) ;

	return get_property_constant ( aa_ID , property_ID  ) ;
}


double 	get_property_constant ( const int  aa_ID,			const int property_ID )
{
	return PhysChemConst [ property_ID  ] [aa_ID ];
}



vector <double>		get_property_constant ( const string		 & seq_chain,  const int property_ID )
{
	vector <int> array_aa_ID =  aa_sequence_to_index_array ( seq_chain );
	return 	get_property_constant ( array_aa_ID,	property_ID );
}

vector <double> 	get_property_constant ( const vector < int > & array_aa_ID,	const int property_ID )
{
	int size = array_aa_ID.size();
	vector < double > array_property ; 	array_property.resize ( size ) ;

	for (int ii=0;ii<size;ii++) 
		array_property[ii] = PhysChemConst [ property_ID  ] [array_aa_ID [ii] ];

	return array_property;
}
