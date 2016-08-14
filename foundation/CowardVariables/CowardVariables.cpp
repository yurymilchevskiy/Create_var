#pragma warning( disable : 4786 )

#include "CowardVariables.h"

#include "../Frequency_extrapolation/Frequency_extrapolation.h"
#include "../Frequency_extrapolation/Frequency_chain_constants.h"


//#include "Sum_of_ready.h"
#include "Dull_Sum.h"
#include "c_FourierSmoothed.h"
/*
#include "Chain_property_log.h"
#include "Chain_normalized_property.h"
#include "c_FourierSmoothed_ml_len.h"
*/
#include "../WiseVariables/Prime_contants_for_chain.h" 

#include "Distance_by_frequency_average.h"
#include "Inv_Distance_by_frequency_average.h"
#include "Log_occurence_differrence.h"
#include "Student_emasculate.h"
#include "Student_emasculate_inv.h"
#include "Mean_compare_probability.h"
#include "Mean_compare_probability_inv.h"


#include <iostream>
#include <fstream>  
#include <cassert> 
#include <sstream>
 

extern ofstream log_stream;

typedef   map < string, Property * > MAP_SEQUENCEMAP_TOVALUE ;
typedef   map < string, int >           MAP_NAME_TO_INDEX ;

CowardVariables:: 
CowardVariables( const string & Sequence_to_values_file_name)
{
	ifstream  in_stream( Sequence_to_values_file_name.c_str() );
	if ( ! in_stream)	
	{	
		cout       << "can't find file "  << Sequence_to_values_file_name << endl;
		log_stream << "can't find file "  << Sequence_to_values_file_name << endl;
		assert (  in_stream);		
		exit (1);	
	}

	vector < string >   all_mapped_definite_task ;

	string current_line;
	int counter = 0;
	while ( getline(in_stream,current_line,'\n' ) )
	{
		if (   current_line[0] == '/'  || 
		       current_line[0] == '#'  || 
		       current_line[0] == ' '  || 
		       current_line[0] == '\n' || 
		       current_line[0] == '\0'  )
		  continue;

		all_mapped_definite_task.push_back( current_line );

		cout << counter << "\t" << current_line << endl;
		counter++;
	} 

	number_of_variables_ = all_mapped_definite_task.size();



	init_known_templates_map ();
	init_frequency_prichindales(all_mapped_definite_task);


	int size = all_mapped_definite_task.size();
///	values_.resize( size );
	number_of_variables_ = size;
 
    MAP_SEQUENCEMAP_TOVALUE ::iterator theIterator_MAP_SEQUENCEMAP_TOVALUE;
	MAP_NAME_TO_INDEX       ::iterator theIterator_MAP_NAME_TO_INDEX;  

	for (int ii=0;ii<size;ii++)
	{
		istringstream ist ( all_mapped_definite_task[ii] );

		string current_procedure_key_word, dummy, current_task_name ;
		ist >> current_procedure_key_word >> dummy >> current_task_name;

		theIterator_MAP_SEQUENCEMAP_TOVALUE = known_templates_map_SequenceMap_toValue_.find ( current_procedure_key_word  );
		if ( theIterator_MAP_SEQUENCEMAP_TOVALUE == known_templates_map_SequenceMap_toValue_.end() )
		{
			log_stream << " WiseSequenceTranslator ERROR: dummy keyword " <<  current_procedure_key_word << endl;
			cout       << " WiseSequenceTranslator ERROR: dummy keyword " <<  current_procedure_key_word << endl;
			exit(1);
		}
		else 
		{
			Property * current_derived_object = known_templates_map_SequenceMap_toValue_[current_procedure_key_word]->clone(all_mapped_definite_task[ii]);
			array_derived_SequenceMap_toValues_.push_back(  current_derived_object );

			if ( current_task_name != "DUMB" && current_task_name != "dumb") 
			{
				if  ( coward_variable_name_to_index_.end() != coward_variable_name_to_index_.find ( current_task_name )  ) 
				{
					log_stream << " WiseSequenceTranslator ERROR: twice assigned variable name: " <<  current_task_name << endl;
					cout       << " WiseSequenceTranslator ERROR: twice assigned variable name: " <<  current_task_name << endl;
					exit(1);

				}
				else 
					coward_variable_name_to_index_ [ current_task_name ] = ii;
			}
		}
	}
}

void  CowardVariables:: 
init_frequency_prichindales(
	const vector < string >   & all_mapped_definite_task)
{
	for (int ii = 0; ii < all_mapped_definite_task.size(); ii++)
	{

		istringstream ist(all_mapped_definite_task[ii]);

		//		Chain_normalized_property 	0 DUMB 	 ZASB820101 1  f
		string dummy_word, word;

		ist >> dummy_word;

		//Dull_Sum		0	DUMB	PTIO830101 - 5 5	1
		//	c_FourierSmoothed	0       DUMB	WERD780101	3.6	8	2

// Тут надо всьтавлять все классы, относящиеся к вычислению предикторов по физ-хим параметрам. 
// Т.е. не использующих Frequency_extrapolation	объекты

		if (dummy_word == "Dull_Sum" ||
			dummy_word == "c_FourierSmoothed")
			continue;

		ist >> dummy_word;
		ist >> dummy_word;

		ist >> word;

		if (frequency_name_to_pull_index_.find(word) == frequency_name_to_pull_index_.end())  // ну не нашлось word
		{
			int current_size = frequency_name_to_pull_index_.size();
			frequency_name_to_pull_index_[word] = current_size;

			Frequency_extrapolation	*curren_frequency = new  Frequency_extrapolation(word, COMMON_USAGE);
			frequency_pull_.push_back(curren_frequency);
		}
	}
}

void   CowardVariables:: 
calc_values (
	  int position ) 
{

	for (int ii=0;ii<number_of_variables_ ;ii++)
	{
		/*double current_value = array_derived_SequenceMap_toValues_[ii]->calc_value (
					position,  
					ii, 
					Chain_Prime_Constants,
					sophisticated_variables    ) ;*/

		double current_value = array_derived_SequenceMap_toValues_[ii]->calc_value(
			position); 

		sophisticated_variables_[position][ii] = current_value;


	}
} 

void   CowardVariables:: 
process_chain ( 
	const	string	&	sequense)
{	
	//pull_fcc_				=	frequency_constant_for_chain(sequense, "together");
	
	frequency_constant_for_chain(sequense, "together");

	//Chain_Prime_Constants_	=	prime_contants_for_chain ( sequense );
	
	prime_contants_for_chain(sequense, Chain_Prime_Constants_);

	int residue_number_in_chain = sequense.size();
	sophisticated_variables_.resize( residue_number_in_chain ); 

	for (int jj=0;jj<residue_number_in_chain; jj++)
	{
		sophisticated_variables_[jj].resize(number_of_variables_); 
		calc_values (jj	) ;
	}
}


void  CowardVariables:: 
init_known_templates_map ()
{
	Dull_Sum*  Dull_Sum_pointer  =									new		Dull_Sum( *this ); ///!!!! like this
	known_templates_map_SequenceMap_toValue_["Dull_Sum"]			=		Dull_Sum_pointer  ;

	c_FourierSmoothed* c_FourierSmoothed_pointer					= new	c_FourierSmoothed (*this);
	known_templates_map_SequenceMap_toValue_["c_FourierSmoothed"]	=		c_FourierSmoothed_pointer; 

	Distance_by_frequency_average*  Distance_by_frequency_average_pointer = new	Distance_by_frequency_average(*this);
	known_templates_map_SequenceMap_toValue_["Distance_by_frequency_average"] = Distance_by_frequency_average_pointer;

	Log_occurence_differrence*  Log_occurence_differrence_pointer = new	Log_occurence_differrence(*this);
	known_templates_map_SequenceMap_toValue_["Log_occurence_differrence"] = Log_occurence_differrence_pointer;

	Inv_Distance_by_frequency_average*  Inv_Distance_by_frequency_average_pointer = new	Inv_Distance_by_frequency_average(*this);
	known_templates_map_SequenceMap_toValue_["Inv_Distance_by_frequency_average"] = Inv_Distance_by_frequency_average_pointer;

	Student_emasculate*  Student_emasculate_pointer = new	Student_emasculate(*this);
	known_templates_map_SequenceMap_toValue_["Student_emasculate"] = Student_emasculate_pointer;

	Student_emasculate_inv*  Student_emasculate_inv_pointer = new	Student_emasculate_inv(*this);
	known_templates_map_SequenceMap_toValue_["Student_emasculate_inv"] = Student_emasculate_inv_pointer;

	Mean_compare_probability*  Mean_compare_probability_pointer = new	Mean_compare_probability(*this);
	known_templates_map_SequenceMap_toValue_["Mean_compare_probability"] = Mean_compare_probability_pointer;

	Mean_compare_probability_inv*  Mean_compare_probability_inv_pointer = new	Mean_compare_probability_inv(*this);
	known_templates_map_SequenceMap_toValue_["Mean_compare_probability_inv"] = Mean_compare_probability_inv_pointer;

}



//vector < Frequency_chain_constants >  CowardVariables::
void CowardVariables::
frequency_constant_for_chain(
	const string & sequence,
	const string & chain_ID)
{

//	vector < Frequency_chain_constants > pull_fcc;

	pull_fcc_.clear();

	int frequency_pull_size = frequency_pull_.size();

	vector < Frequency_chain_constants > xxx;

	for (int ii = 0; ii<frequency_pull_size; ii++)
	{

		ifstream data_stream;
		frequency_pull_[ii]->i_freq_data_stream(
			chain_ID, data_stream);

		Frequency_chain_constants single_fcc =
			frequency_pull_[ii]->get_single_frequency_chain_constants(
				sequence, data_stream);

		pull_fcc_.push_back(single_fcc);

	}

//	return  pull_fcc;
}
