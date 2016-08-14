#pragma warning( disable : 4786 )

#include "get_property_constant.h"
#include "property_constant.h"

#include "../aminoacid_to_index.h" 
#include "../aa_sequence_to_index_array.h"

//#include "Frequency.h"

//vector < Frequency * >  init_Frequency_object_array ();

using namespace std;

 void prime_contants_for_chain ( const string & chain_sequence,vector < vector <double > > & Chain_Prime_Constants  )
{
	Chain_Prime_Constants.clear();
	vector <int> chain_index = 	aa_sequence_to_index_array ( chain_sequence );

//	vector < vector <double > > pr_cnst_set ; 	
	
	for (int ii=0; ii< number_of_Properties ;ii++)

	{
		vector <double> current_chain_property = get_property_constant ( chain_index  ,	ii );
		Chain_Prime_Constants.push_back(current_chain_property );
	}
/*
	static vector < Frequency * >  Frequency_object_array  = init_Frequency_object_array ();

	for ( ii = 0; ii < number_of_Frequency_object; ii++ ) 
	{
		vector < int > processed_index = Frequency_object_array [ii]->translate_sequence_to_degenerate_array( chain_sequence );

		int number_of_classes = Frequency_object_array [ii]->get_number_of_classes();
		for (int kk=0; kk < number_of_classes ;kk++ )
		{
			vector < double > current_frequneces = Frequency_object_array [ii]->get_class_frequenses (processed_index, kk ) ;
			pr_cnst_set.push_back ( current_frequneces ) ;
		}

		vector < double >	occurence_vector	= Frequency_object_array [ii]->get_occurence_vector	( processed_index ) ;
		pr_cnst_set.push_back  (occurence_vector) ;

		vector < double > local_quality_vector  = Frequency_object_array [ii]->get_local_quality_vector  ( processed_index ) ;
		pr_cnst_set.push_back   (local_quality_vector) ;
	}
*/
	//return pr_cnst_set ;
}
