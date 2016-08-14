// то же что и Frequency_extrapolation nтолько экономит память 
// вместо double	**distance_record_;  будет map <int, double* >
#ifndef FREQUENCY_EXTRAPOLATION_ZIP_H
#define FREQUENCY_EXTRAPOLATION_ZIP_H

#include "Frequency_extrapolation_operating_modes.h"
#include "Frequency_chain_constants.h"

#include <vector>
#include <fstream>
#include <string>
#include <fstream>
#include <map>

using namespace std;

class Curly_calculus_system;
class Sheduler;
class Fragment_base_subtle ;

using namespace std;


class Frequency_extrapolation_ZIP
{
public:
	
	Frequency_extrapolation_ZIP	( const string & name, const Frequency_extrapolation_operating_modes run_mode  );
	~Frequency_extrapolation_ZIP() ; /// FIX 

		void fill_up_frequency_extrapolation ();

    vector < vector <string> >		get_degeneration						() const {return degeneration_; } ;

	vector < int > 					translate_sequence_to_degenerate_array  ( const string &sequence ) ;
	vector < int >					translate_sequence_to_degenerate_array  ( const vector < int > & crude_sequence_index ) ;
  
	vector < vector < string > > 	translate_array_of_cursor_to_sequence	( vector < int >  & degenerate_array );
	vector < string >  				translate_cursor_to_sequence			( int  degenerate_cursor );


	void show_freq_data (const string pdb_chain_ID );
	void show_freq_data_slow ();


	vector <double>   get_single_distance_record (
		const int record_index ,
		ifstream & data_stream);

//	vector <int>   	get_occurence (		ifstream & data_stream);
	int   get_occurence (const int record_index );
		



	vector <double>   	get_single_distance_to_clusters_sum					 ( const int record_index );
	vector <double>   	get_single_squared_distance_to_clusters_sum			 ( const int record_index );
	vector <double>   	get_single_inverse_distance_to_clusters_sum			 ( const int record_index );
	vector <double>   	get_single_inverse_squared_distance_to_clusters_sum  ( const int record_index );



	vector <double>   	get_tot_distance_to_clusters_sum				 ( );
	vector <double>   	get_tot_squared_distance_to_clusters_sum		 ( );
	vector <double>   	get_tot_inverse_distance_to_clusters_sum		 ( );
	vector <double>   	get_tot_inverse_squared_distance_to_clusters_sum ( );
	int					get_total_sample_size							 ( ) ;


	int get_number_of_elements () const  {return number_of_elements_; }


	void prepare_together_freq_data ();

	void o_freq_data_stream ( 	string  base_file_name, ofstream & data_stream );
	void i_freq_data_stream (  	string  base_file_name, ifstream & data_stream );

	Frequency_chain_constants 	get_single_frequency_chain_constants ( 
									const string & sequence, 
									ifstream & data_stream ) ;
private:
	Fragment_base_subtle	*fbs_;
	Sheduler				*sheduler_;
	Curly_calculus_system	*curly_calculus_system_;

	string					name_;

	int shift_;					
	vector < vector <string> >	degeneration_;


	int number_of_elements_ ;							//number of elements in Frequency_extrapolation_base
	int number_of_classes_  ; 
	int record_length_		;
	int number_of_items_in_record_;

	void init_curly_calculus_system						();
	void read_degeneration								( const string & degeneration_file_name );
	void translate_degeneration_to_index				();
	void assign_window_size_left_right_border_value 	();


	int							window_size_;			// number of sequential residues under consideration 
	int							left_border_value_;		// obtainded from position_shift_ (for performance increase only)
	int							right_border_value_;	// obtainded from position_shift_ (for performance increase only)
	
	vector < int >				position_shift_;		// position_shift_ - assigns position for each point of window 
	vector < vector < int > >   index_of_degeneration_; // 


	vector <int > claster_motif_index_ ;
	double		**claster_motif_coordinates_;

	void init_claster_motif(); 

	int fragment_length_;						// number of aminoacids in polypeptide 

	int		*occurence_;	
	double	**distance_record_;

	void 	write_single_chain_data (
		const string & pdb_chain_ID,
		int *occurence ,
		double **distance_record );

	void _write_single_chain_data  (
		const string & pdb_chain_ID);


	void _such_up_single_chain_data (
		const string & pdb_chain_ID );


	map <int,int> index_to_non_zerro_occurence_index_;

	void 	init_map_index_to_nonzerro_index ();

	char *zakroma_;

	void _suck_up_kruto (ifstream & data_stream);
	int *index_to_pos_;

	int non_zerro_occurence_number_;

// запасной поток для файла 
	ifstream standby_datastream_; 

};

#endif
