#ifndef COWARDVARIABLES_H
#define COWARDVARIABLES_H

#include "Property.h"

#include	<map>
#include	<vector>

using namespace std;

#include "../Frequency_extrapolation/Frequency_chain_constants.h"

class Frequency_extrapolation;

class CowardVariables
{
public:
CowardVariables() {}
	explicit CowardVariables( const string & Sequence_to_values_file_name);
	~CowardVariables() {};
/*

	vector < vector < vector < double > > >		process_protein (  
		Model & model,
		const  vector < vector < vector < double > > >   Protein_Prime_Constants);
*/

	void process_chain ( 
		const string & sequense );

	int get_number_of_variables () const {return number_of_variables_ ;}  

	map    < string, int >					get_coward_variable_name_to_index () const
												{ return coward_variable_name_to_index_;}  


	vector < vector < double > >	const & get_sophisticated_variables() const { return sophisticated_variables_; }

	vector < vector <double > >     Chain_Prime_Constants_;
	vector < vector < double > >	sophisticated_variables_;

	vector < Frequency_chain_constants >		pull_fcc_;
	vector < Frequency_extrapolation * >		frequency_pull_;
	map    < string, int >						frequency_name_to_pull_index_;

	/*vector < Frequency_chain_constants >
		frequency_constant_for_chain(
			const string & sequence,
			const string & chain_ID);*/

	void 		frequency_constant_for_chain(
			const string & sequence,
			const string & chain_ID); 

private:

	vector < Property *  >				    array_derived_SequenceMap_toValues_;
	map    < string, Property* >			known_templates_map_SequenceMap_toValue_;
	map    < string, int >					coward_variable_name_to_index_;

	int		number_of_variables_ ;


	void init_frequency_prichindales(const vector < string >   & all_mapped_definite_task);


	void  init_known_templates_map ();

	/*void calc_values (
	  int position, 
	  const  vector < vector < double > >   & Chain_Prime_Constants,
	  vector < vector < double > >			& sophisticated_variables) ;*/

	void calc_values(		int position); 

	CowardVariables ( const CowardVariables& );
	void operator = ( const CowardVariables& );
};

#endif
