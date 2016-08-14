//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!


#pragma warning( disable : 4786 )

#include "Frequency_extrapolation.h"
#include "Curly_calculus_system.h"
#include "../aa_sequence_to_index_array.h"
#include "fill_up_appropriate_crude_index_subset.h"

#include "../Fragment_base/Fragment_base_subtle.h"

#include <iostream>

#include "../CommonFunc.h"
#include "../Censorship.h"
#include "../Sheduler.h"

#include "../Fragment_base/Chain_binary.h"
#include "../Cluster_set/Cluster_set.h"
#include "../by_Qmol/kabsch_stolen.h"

#include "../Sane_metrics/Sane_metrics.h"

#include "../make_claster_motif_by_dihedral_set.h"

#include "../Fragment_base/accepted_chain_data.h"

#include <vector>
#include <cassert>
#include <iostream>

/*
void make_claster_motif_by_dihedral_set( 
	const string & path_to_dihedral_store,
	int number_of_classes,
	int fragment_length);
*/

extern ofstream log_stream;
extern Censorship configuration;

double *local_distance_record;

Frequency_extrapolation::
~Frequency_extrapolation	()
{
	if ( fbs_ )			delete fbs_;
	if ( sheduler_ ) 		delete sheduler_;
	if ( claster_motif_coordinates_ ) 		
	{
		for (int kk=0; kk<number_of_classes_; kk++)			delete [] claster_motif_coordinates_[kk];
		delete [] claster_motif_coordinates_;
	}

	if (  occurence_ )		delete [] occurence_;

	if ( zakroma_ )
		delete [] zakroma_;

	if ( index_to_pos_ )
		delete [] index_to_pos_;


	if ( san_me_)
		delete san_me_;
}

Frequency_extrapolation::
Frequency_extrapolation	( 
	const string & name, 
	const Frequency_extrapolation_operating_modes run_mode ):
//	const Cluster_source_operating_modes ):
 name_ (name),
 sheduler_					(0),
 curly_calculus_system_		(0),
 fbs_						(0),
 san_me_					(0),
 claster_motif_coordinates_ (0),
 occurence_					(0),	
// distance_record_			(0),
 zakroma_					(0),
 index_to_pos_				(0),
 cls_						(0)
{
	//string ttt =configuration.option_meaning("Path_to_Frequency_extrapolation")  + name_ + string ("/") + string ("sheduler");
	sheduler_		= new Sheduler  (configuration.option_meaning("Path_to_Frequency_extrapolation")  + name_ + string ("/") + string ("sheduler") ) ;

/// !!! а нужен таки Cluster_set

//	path_to_accepted_chain_list_  = 		sheduler_->option_meaning("LEARNING_SAMPLE_BINARY_FILE");	
	
	denominator_constant_ =   atof( sheduler_->option_meaning("DENOMINATOR_CONSTANT").c_str() );

	string cluster_set_name_ =   sheduler_->option_meaning("CLUSTER_SET_NAME");
	cls_ = new Cluster_set  (	cluster_set_name_,COMMON_USAGE_CLUSTER_SET_MODE);

//тут лежат заданные руками базизные структуры ТОЛЬКО ЕСЛИ УГЛЫ ЗАДАЮТСЯ ВРУЧНУЮ !!!!! 
//	path_to_dihedral_store_ =   sheduler_->option_meaning("PATH_TO_DIHEDRAL_STORE");
	

	string sane_metrics_name   = cls_->get_sane_metrics_mame();

	number_of_classes_ = cls_->number_of_classes();

/*
	Sane_metrics *san_me= new  Sane_metrics ( sane_metrics_name,SANE_METRICS_COMMON_USAGE);	
	san_me_= new  Sane_metrics ( sane_metrics_name,SANE_METRICS_COMMON_USAGE);
	string fbs_file_name = san_me_->get_fragment_base_subtle_name();
	fbs_= new  Fragment_base_subtle (fbs_file_name ,FRAGMENT_BASE_SUBTLE_COMMON_USAGE)  ;
*/
	fragment_length_ = cls_->fragment_length();
	shift_ = fragment_length_/2;


	claster_motif_coordinates_  = cls_->get_claster_motif_coordinates();

//***********************************************************************

	read_degeneration ( 
		configuration.option_meaning("Path_to_Frequency_extrapolation")  + name_ + string ("/") + sheduler_->option_meaning ("DEGENERATION_ASSIGNMENT_FILE") ) ;

	translate_degeneration_to_index ();
	assign_window_size_left_right_border_value ();
	init_curly_calculus_system ();

	number_of_elements_  = curly_calculus_system_	->	get_number_of_elements();

	number_of_items_in_record_ = 4*number_of_classes_ ; 

	//init_claster_motif();

	record_length_		 = number_of_items_in_record_ * sizeof ( double ) ;


	ifstream  data_stream;

	switch ( run_mode  )
	{
		case COMMON_USAGE: 
			
	//		i_freq_data_stream ("together",	data_stream);
			//i_freq_data_stream ("together",	standby_datastream_);
//??			init_claster_motif();

			soak_up_zakroma_from_binary_store();
			init_essential_indexes();

			
			//_suck_up_kruto (data_stream);   // FIX    вместо нее init_essential_indexes();

			break;
		case FILL_UP_FREQUENCY_EXTRAPOLATION_STORE:

			occurence_ = new int [number_of_elements_];
			memset (occurence_ ,0, number_of_elements_*sizeof (int)  );

		// теперь тут, чтобы отличать оот получения вместо кластеров заранее заданных структур 
//			init_claster_motif();
			fill_up_frequency_extrapolation		();
//			prepare_united_chain_data			();
//			prepare_jack_nife_chain_data		();
			break;
/*
		case FILL_UP_FOR_PRELIMINARILY_SETTET_STRUCTURES:
			occurence_ = new int [number_of_elements_];
			memset (occurence_ ,0, number_of_elements_*sizeof (int)  );


		claster_motif_coordinates_ = new double* [number_of_classes_] ;	
		for ( int ii=0; ii < number_of_classes_; ii++  )
			claster_motif_coordinates_[ii]  = new double [fragment_length_*9];


			make_claster_motif_by_dihedral_set( 
				path_to_dihedral_store_,
				number_of_classes_,
				fragment_length_,
				claster_motif_coordinates_); 
			fill_up_frequency_extrapolation		();

		case ANALYSE_CLUSTER_ONLY: 
//			calc_cluster_mutual_distance (distance_mutual);
			break;
*/
		default :
			log_stream	<< "Inadmissible run mode for Frequency class" << endl;
			cout		<< "Inadmissible run mode for Frequency class" << endl;
			exit (-1);
	}
}

void Frequency_extrapolation:: 
init_claster_motif() 
{
//	string claster_origin_structure_file =  sheduler_->option_meaning("CLASTER_ORIGIN_STRUCTURE_LIST");

	string cluster_set_mame = cls_->get_cluster_set_mame();

	string path_to_cluster_set_protocol_file_name = 
		configuration.option_meaning("Path_to_Cluster_set")  +  
		cluster_set_mame + string ("/") + string ("protocol");

	claster_motif_index_ = pull_out_claster_origin_structure_list (
								path_to_cluster_set_protocol_file_name	);
	
	claster_motif_coordinates_ = new double* [number_of_classes_] ;	
	for ( int ii=0; ii < number_of_classes_; ii++  )
		claster_motif_coordinates_[ii]  = new double [fragment_length_*9];

	for (int kk=0;kk<number_of_classes_;kk++)
		fbs_->get_coord( claster_motif_index_ [kk], claster_motif_coordinates_[kk]);
}

double **  Frequency_extrapolation:: 
get_claster_motif_coordinates ()
{
	return claster_motif_coordinates_;
}




void Frequency_extrapolation:: 
fill_up_frequency_extrapolation  ()
{

	
	vector <string> accepted_chain_ID_list;
	vector <int >	accepted_chain_lenth;
	
	string binary_file_name ("accepted_chain_list.bin");

	fill_up_accepted_chain_data ( 
		accepted_chain_ID_list, 
		accepted_chain_lenth,
		binary_file_name);


//	vector < string > accepted_chain_ID_list = fbs_->get_accepted_chain_ID_list ();

	double *cu_cord_set  = new double [fragment_length_*9];

	double *w	= new double [fragment_length_ *9];

	for (int kk=0;kk<(fragment_length_ *9);kk++)
		w[kk] = 1;



	// ******* REFRESH CURRENT ARRAYS ***********
//		vector < int > occurence;							occurence.	    resize ( number_of_elements_  );
//		vector < vector < double > >  distance_record;		distance_record.resize ( number_of_elements_  );

//	int *occurence_global = new int [number_of_elements_];
//	memset (occurence_global ,0, number_of_elements_*sizeof (int)  );

/***
	double **distance_record_global = new double* [number_of_elements_];
	for ( kk=0; kk<number_of_elements_;kk++ )
	{
		distance_record_global [kk] = new double [number_of_items_in_record_ ];
		memset (distance_record_global[kk],0,number_of_items_in_record_ *sizeof (double)) ;
	}
**/

	memset (occurence_ ,0, number_of_elements_*sizeof (int)  );
	for (unsigned ii = 0; ii < accepted_chain_ID_list.size(); ii++ )
	{


		// переделать только lkz создания together database 
		// Тогда это блок не нужен
/*
		{
		string freq_data_file_name = 
			configuration.option_meaning("Path_to_Frequency_extrapolation") + 	
			name_			   + 
			"/base/"		   + 
			accepted_chain_ID_list[ii]	   +  
			".freq_data";


			ifstream is_present_stream_test(freq_data_file_name.c_str());
			if (is_present_stream_test )
			{
				is_present_stream_test.close();
				cout << "Uze bilo "<< ii << "   " << accepted_chain_ID_list[ii] << endl;
				
				continue;   // FIX
			}
		}

		memset (occurence_ ,0, number_of_elements_*sizeof (int)  );
*/

//		for ( int ttt=0; ttt<number_of_elements_;ttt++ )			memset ( distance_record_[ttt],0,4*number_of_classes_ *sizeof (double) );


//		Single_masticated_structure_record * record = 	mss_->read_single_record ( PDB_chain_ID_list_[ii] );

		Chain_binary chain( accepted_chain_ID_list[ii] );
		string sequence =chain.get_sequence();

		bool *is_geometry_admissible	=	chain.get_is_geometry_admissible();
		bool *is_there_coord			=	chain.get_is_there_coord();

//		string sequence = record->get_sequence();
		vector < int > processed_index  = translate_sequence_to_degenerate_array(sequence) ;



		// ***** CREATING CURRENT CONTRTIBUTION to DATABASE DATA *****/
		// int current_number_of_fragments = set_of_coordinate_in_clasters_system.size(); 
// пока предполагаем что рамер окна из current.degeneration и размер фрагмента равны
		int current_number_of_fragments = sequence.length() - window_size_ + 1;
		for (int kk=0; kk< current_number_of_fragments ;kk++ )
		{
			
			chain.extract_fragment(kk,fragment_length_,cu_cord_set );

// check that coorditates  presents and admissible
			bool is_geom_flag = true;
			bool is_coord_flag = true;
			//int right_stop = min (kk + window_size_,current_number_of_fragments ) ;
			int right_stop = kk + window_size_;
			for (int pp=kk; pp< right_stop; pp++)
			{
				is_geom_flag  = is_geom_flag & is_geometry_admissible	[pp];
				is_coord_flag = is_coord_flag & is_there_coord			[pp];
			}
			if (!is_geom_flag || !is_coord_flag) 
				continue ;			// енсли плохие координаты или  их нету то пропускаем

			occurence_					[ processed_index[kk+ shift_] ] ++;
//			occurence_global			[ processed_index[kk+ shift_] ] ++;

			
			if ( distance_record_.find(processed_index[kk+ shift_]) == distance_record_.end()  )  // ну не нашлось word
			{
				distance_record_[ processed_index[kk+ shift_] ] = new double [number_of_items_in_record_ ];
				memset (distance_record_[ processed_index[kk+ shift_] ],0,number_of_items_in_record_*sizeof(double));   /// первый раз
//				local_distance_record = new double [number_of_items_in_record_ ];
//				memset (local_distance_record,0,number_of_items_in_record_*sizeof(double));
				
			}	
			//memset (distance_record_[ processed_index[kk+ shift_] ],0,number_of_items_in_record_*sizeof(double));   /// ТУТ БЫЛА ОЩИбка из-за которой together databse врала

			for (int zz=0;zz<number_of_classes_;zz++)
			{

				bool error_flag= true;
				double distance  = kabsch_rmsd	(
										cu_cord_set , 
										claster_motif_coordinates_[zz], 
										w, 
										(3*fragment_length_),
										error_flag);


				double xx = distance;
//				double yy = 1.0 / (1 + xx );

// АРХИВАЖНО
				double yy = 1.0 / (denominator_constant_ + xx );

				
				distance_record_ [ processed_index[kk+ shift_] ] [zz]  +=  xx;
				distance_record_ [ processed_index[kk+ shift_] ] [zz +   number_of_classes_ ]  +=  xx*xx;

				distance_record_[ processed_index[kk+ shift_] ] [zz  + 2*number_of_classes_ ]  +=  yy;
				distance_record_[ processed_index[kk+ shift_] ] [zz  + 3*number_of_classes_ ]  +=  yy*yy;

/***

				distance_record_global [ processed_index[kk+ shift_] ] [zz]  +=  xx;
				distance_record_global [ processed_index[kk+ shift_] ] [zz +   number_of_classes_ ]  +=  xx*xx;

				distance_record_global[ processed_index[kk+ shift_] ] [zz  + 2*number_of_classes_ ]  +=  yy;
				distance_record_global[ processed_index[kk+ shift_] ] [zz  + 3*number_of_classes_ ]  +=  yy*yy;
****/
			}
		}
///////		
//			write_single_chain_data (
//				accepted_chain_ID_list[ii],
//				occurence_,
//				distance_record_);


//			_write_single_chain_data (	accepted_chain_ID_list[ii]);


		cout << ii << "   " << accepted_chain_ID_list[ii] << endl;
	}

	_write_single_chain_data (	"together");
/***
	memcpy (occurence_,occurence_global, number_of_elements_ * sizeof (int) ) ;

	for (  kk=0; kk<number_of_elements_;kk++ )
		memcpy (distance_record_[kk],distance_record_global[kk],number_of_items_in_record_ *sizeof (double) );

	_write_single_chain_data (	"together");
***/

//elete [] occurence_global ;
	delete [] w;	delete [] cu_cord_set ;

//or ( kk=0; kk<number_of_elements_;kk++ )
//delete [] distance_record_global[kk];
//elete [] distance_record_global;
}


void Frequency_extrapolation::
o_freq_data_stream (  
	string  base_file_name, 
	ofstream & data_stream )
{
	string freq_data_file_name = 
		configuration.option_meaning("Path_to_Frequency_extrapolation") + 	
		name_			   + 
		"/base/"		   + 
		base_file_name	   +  
		".freq_data";

	data_stream.open(freq_data_file_name .c_str() ,ios::binary);
	if ( ! data_stream)	
	{	
		log_stream << "write_single_chain_data(): ERROR -  can't create " << freq_data_file_name << endl;
		cout		<< "write_single_chain_data(): ERROR -  can't create " << freq_data_file_name << endl;
		exit (1);	
	}
}

void Frequency_extrapolation::
i_freq_data_stream (  
	string  base_file_name, 
	ifstream & data_stream )
{
	string freq_data_file_name = 
		configuration.option_meaning("Path_to_Frequency_extrapolation") + 	
		name_			   + 
		"/base/"		   + 
		base_file_name	   +  
		".freq_data";

	data_stream.open(freq_data_file_name .c_str() ,ios::binary);
	if ( ! data_stream)	
	{	
		log_stream <<  "ERROR -  can't create " << freq_data_file_name << endl;
		cout		<< "ERROR -  can't create " << freq_data_file_name << endl;
		exit (1);	
	}
}




// ***********************************************************************************************************************************
void Frequency_extrapolation::
_such_up_single_chain_data (
	const string & pdb_chain_ID )
{
	ifstream data_stream; 
	i_freq_data_stream (pdb_chain_ID,data_stream) ;

	memset (occurence_ ,0, number_of_elements_*sizeof (int)  );
	for ( int kk=0; kk<number_of_elements_;kk++ )
		memset (distance_record_[kk],0,number_of_items_in_record_ *sizeof (double) );


	int non_zerro_occurence_number;
	data_stream.read ( (char* ) & non_zerro_occurence_number ,sizeof (int)  );  /******/
	int kk;
	for ( kk=0; kk<non_zerro_occurence_number; kk++ )
	{
		int ii;
		data_stream.read ( (char* ) & ii,													sizeof (int)    );  /******/
		data_stream.read( (char* ) & occurence_ [ii],										sizeof (int)    );  /******/
		data_stream.read( (char* ) distance_record_ [ii],	number_of_items_in_record_* sizeof (double) );  /******/

		log_stream << kk << " "  << ii << "\t" << occurence_ [ii] <<  "\t" << distance_record_ [ii][0]/occurence_ [ii] << endl;
	}

}
	



void Frequency_extrapolation::
read_degeneration	( const string & degeneration_file_name)
{
	ifstream  tune_stream( degeneration_file_name.c_str() );
	if ( ! tune_stream)	{	
		log_stream << "Frequency: ERROR -  can;t find tune file" << degeneration_file_name << endl;
		cout       << "Frequency: ERROR -  can;t find tune file" << degeneration_file_name << endl;
		exit (1);	
	}

	string current_line;

	while( getline( tune_stream  , current_line, '\n' ) )
	{
		if (   current_line[0] == '/'  || 
			   current_line[0] == '#'  || 
//			   current_line[0] == ' '  || 
			   current_line[0] == '\n' || 
			   current_line[0] == '\0')
			continue;

		int dummy_int;
		string dumy_word;  
		vector < string > single_degeneration;

		istringstream ist (current_line);
		ist >> dummy_int;
		position_shift_.push_back(dummy_int) ;
		while (	ist >> dumy_word )  
			single_degeneration.push_back( dumy_word ) ;

		degeneration_.push_back(single_degeneration) ;
	}
}
vector < int > Frequency_extrapolation::
translate_sequence_to_degenerate_array ( const string &sequence ) 
{
	vector < int > crude_index_array = aa_sequence_to_index_array ( sequence ); 

	return ( translate_sequence_to_degenerate_array	( crude_index_array  ) );
}

vector < int > Frequency_extrapolation::
translate_sequence_to_degenerate_array  ( const vector < int > & crude_index_array ) 
{
	vector < int > index_subset_array;   index_subset_array.resize( window_size_ );

	int sequence_size = crude_index_array.size(); 
	vector < int > degenerate_array ;  degenerate_array.resize(sequence_size); 
	
	for ( int ii=0;ii<sequence_size;ii++ )
	{

		fill_up_appropriate_crude_index_subset (
			index_of_degeneration_,
			ii,
			left_border_value_,	
			right_border_value_, 
			crude_index_array,
			get_virtual_residue_index (),
			index_subset_array); 

			int cursor = curly_calculus_system_->get_cursor_by_array ( index_subset_array );

			degenerate_array [ii] = cursor;
	}

	return degenerate_array ;
}

void Frequency_extrapolation::
init_curly_calculus_system ( )
{
	vector < int >  bases;    bases.resize( window_size_ ) ;
	
	for ( int ii=0; ii< window_size_; ii++ )
		bases[ii] =  degeneration_[ii].size();

	curly_calculus_system_ = new Curly_calculus_system ( bases );
}
void Frequency_extrapolation::
assign_window_size_left_right_border_value ()
{
	window_size_		= position_shift_.size();
	left_border_value_	= position_shift_.front();  
	right_border_value_ = position_shift_.back();
}

void Frequency_extrapolation::
translate_degeneration_to_index ()
{
	map <char, int> temp_map =  set_index_by_aa() ;  // to get size only
	int map_size = temp_map.size(); 

	int degeneration_size = degeneration_.size();
	index_of_degeneration_.resize( degeneration_size );

	int ii;
	for ( ii=0;ii<degeneration_size;ii++ )
		index_of_degeneration_[ii].resize (map_size );

	for ( ii=0;ii<degeneration_size;ii++ )
	{
		int current_degeneration_size = degeneration_[ii].size();
		for ( int jj=0;jj<current_degeneration_size ;jj++)
		{
			int current_word_size = degeneration_[ii][jj].size();
			for (int kk=0;kk< current_word_size  ;kk++)
			{
  				const char aa =  degeneration_[ii][jj][kk];
				int aa_index = aminoacid_to_index ( aa );
				index_of_degeneration_[ii][ aa_index ] = jj;

			}
		}
	}

}


vector < string >  Frequency_extrapolation::
translate_cursor_to_sequence (  int   degenerate_cursor )
{

		vector < string >  redundant_sequence;
		vector < int > current_array = curly_calculus_system_->get_array_by_cursor	( degenerate_cursor  );
		int current_array_size = current_array.size();
		for (int kk=0;kk<current_array_size ;kk++ )
		{
			int current_index = current_array [kk];
			string current_word = degeneration_ [kk][ current_index  ];
			redundant_sequence.push_back ( current_word  ); 
		}

		return redundant_sequence;
}






void Frequency_extrapolation::
init_map_index_to_nonzerro_index ()
{
	ifstream data_stream; 
	i_freq_data_stream ("together",data_stream) ;

	int non_zerro_occurence_number;
	data_stream.read (  (char* ) & non_zerro_occurence_number, sizeof (int) );

	int shift;
	for (int kk = 0; kk < non_zerro_occurence_number;kk++)
	{
		shift = sizeof(int) + (record_length_ + 2* sizeof (int) ) * kk;
		data_stream.seekg(shift,ios::beg);
		int ii;
		data_stream.read (  (char* ) & ii, sizeof (int) );

		index_to_non_zerro_occurence_index_[ii] = kk; 

		log_stream << kk << "  " << ii << endl;
		
	}
	data_stream.clear();
	data_stream.close();


}


void Frequency_extrapolation::
show_freq_data_slow ()
{
	string show_file_name = 
		configuration.option_meaning("Path_to_Frequency_extrapolation") + 	
		name_			+ 
		"/base/"		+ 
		"together"		+  
		".show_freq_data_slow";

	ofstream out(show_file_name.c_str() );
	if ( ! out )	
	{	
		log_stream << "write_single_chain_data(): ERROR -  can't create " << show_file_name<< endl;
		log_stream << "write_single_chain_data(): ERROR -  can't create " << show_file_name<< endl;
		exit (1);	
	}
/*
	ifstream data_stream; 
	i_freq_data_stream ("together",data_stream) ;
*/
/// FIX STRANGE
	int kk;
	for ( kk=0;kk<number_of_elements_;kk++) 
	{
/*			int dd;
			if (kk==74414) 
				cin >> dd;
*/
 		int occurence_kk = get_occurence (kk);

		if (occurence_kk == 0 )
			continue;

		PutVa (kk,out, 8,2,'l');
		out << ": ";
		PutVa (occurence_kk,out, 6,2,'l');

		vector <double> single_distance_to_clusters_sum  
			= get_single_distance_to_clusters_sum(kk);
	
		for (int tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_distance_to_clusters_sum[tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}

		out << " || ";


		vector <double> single_squared_distance_to_clusters_sum  = get_single_squared_distance_to_clusters_sum(kk);
		int tt;
		for (tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_squared_distance_to_clusters_sum[tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}


		vector <double> single_inverse_distance_to_clusters_sum  = get_single_inverse_distance_to_clusters_sum(kk);
		for (tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_inverse_distance_to_clusters_sum[tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}

		vector <double> single_inverse_squared_distance_to_clusters_sum  = get_single_inverse_squared_distance_to_clusters_sum(kk);
		for (tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_inverse_squared_distance_to_clusters_sum [tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}


		vector < string > redundant_words	= translate_cursor_to_sequence	( kk );
		for ( unsigned jj=0;jj< redundant_words.size();jj++ )
			PutVa ( redundant_words[ jj ],out ,redundant_words[ jj ].size()+5,redundant_words[ jj ].size()+2,'l');

		out << endl;

	}

	out << "************************************************************************"<< endl;
	{
		int occurence_kk = get_total_sample_size(); 

		PutVa (kk,out, 8,2,'l');
		out << ": ";
		PutVa (occurence_kk,out, 12,2,'l');

		vector <double> single_distance_to_clusters_sum  
			= get_tot_distance_to_clusters_sum();
			//get_single_distance_to_clusters_sum(kk);
	
		for (int tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_distance_to_clusters_sum[tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}

		out << " || ";


		vector <double> single_squared_distance_to_clusters_sum  = 
			get_tot_squared_distance_to_clusters_sum();
			//get_single_squared_distance_to_clusters_sum(kk);
		int tt;
		for (tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_squared_distance_to_clusters_sum[tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}


		vector <double> single_inverse_distance_to_clusters_sum  = 
			get_tot_inverse_distance_to_clusters_sum();
			//get_single_inverse_distance_to_clusters_sum(kk);

		for (tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_inverse_distance_to_clusters_sum[tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}

		vector <double> single_inverse_squared_distance_to_clusters_sum  = 
			get_tot_inverse_squared_distance_to_clusters_sum();
			//get_single_inverse_squared_distance_to_clusters_sum(kk);

		for (tt=0;tt<number_of_classes_;tt++) 
		{
			double temp = ( occurence_kk > 0 ) ? single_inverse_squared_distance_to_clusters_sum [tt] / occurence_kk  : 0;
			PutVaDouble ( temp, out, 6,2,'l');
		}
	}
}


// ***********************************************************************************************************************************
// Надо написать функцию типа write_single_chain_data () чтобы прописывались только те записи для которых occurence_[ii] != 0;
// ***********************************************************************************************************************************
void Frequency_extrapolation::
_write_single_chain_data (
	const string & pdb_chain_ID )
{
	ofstream data_stream; 
	o_freq_data_stream (pdb_chain_ID,data_stream) ;

	int non_zerro_occurence_number = 0;
	int ii;
	for ( ii=0;ii<number_of_elements_;ii++)
	{
		if ( occurence_ [ii] != 0 )
			non_zerro_occurence_number++;
	}

	int *index_to_pos = new int [number_of_elements_];
	for (ii=0;ii<number_of_elements_;ii++)
		index_to_pos [ii] = -1;


	double *accumulated_distance_record = new double [number_of_items_in_record_];
	memset (accumulated_distance_record ,0,number_of_items_in_record_*sizeof (double));


	data_stream.write ( (char* ) & non_zerro_occurence_number ,sizeof (int)  );  /******/

	int counter = 0;
	for (ii=0;ii<number_of_elements_;ii++)
	{
		if ( occurence_ [ii] == 0 )
			continue;

		data_stream.write ( (char* ) & ii,													sizeof (int)    );  /******/
		data_stream.write ( (char* ) & occurence_ [ii],										sizeof (int)    );  /******/
		data_stream.write ( (char* ) distance_record_ [ii],	number_of_items_in_record_* sizeof (double) );  /******/
	
		index_to_pos [ii] = counter;
		counter++;

		for (int kk=0;kk<number_of_items_in_record_;kk++)
			accumulated_distance_record [kk] +=  distance_record_[ii][kk];
	}

	data_stream.write ( (char* ) index_to_pos,	number_of_elements_ * sizeof (int)    );  /******/
	data_stream.write ( (char* ) accumulated_distance_record ,	number_of_items_in_record_* sizeof (double) );  /******/

	int total_sample_size = 0;
	for (ii=0;ii<number_of_elements_;ii++)		total_sample_size += occurence_[ii];
	data_stream.write ( (char* ) & total_sample_size  ,sizeof (int)  );

	
	log_stream <<  pdb_chain_ID << endl;
	log_stream <<  accumulated_distance_record[0] << endl;
	log_stream <<  accumulated_distance_record[0]*accumulated_distance_record[0] << endl;
	log_stream <<  accumulated_distance_record[30] << endl;
	log_stream <<  total_sample_size << endl;
	log_stream <<  "sigma: " << total_sample_size*accumulated_distance_record[30] - accumulated_distance_record[0]*accumulated_distance_record[0] << endl;
	log_stream <<  "**********************************" << endl;



	
	delete [] accumulated_distance_record ;
	delete [] index_to_pos;
}


void Frequency_extrapolation::
init_essential_indexes()
{
	i_freq_data_stream("together", standby_datastream_);

	standby_datastream_.seekg(0, ios::beg);

	standby_datastream_.read((char*)&non_zerro_occurence_number_, sizeof(int));

	index_to_pos_ = new int[number_of_elements_];
	int shift = non_zerro_occurence_number_* (record_length_ + 2 * sizeof(int)) + sizeof(int);
	standby_datastream_.seekg(shift, ios::beg);
	standby_datastream_.read((char*)index_to_pos_, number_of_elements_ *sizeof(int));

}

void Frequency_extrapolation::
soak_up_zakroma_from_binary_store()
{
	string freq_data_file_name =
		configuration.option_meaning("Path_to_Frequency_extrapolation") +
		name_ +
		"/base/" +
		"together.freq_data";

	ifstream t_stream(freq_data_file_name.c_str(), ios::binary);
	if (!t_stream)
	{
		log_stream << "ERROR -  can't create " << freq_data_file_name << endl;
		cout << "ERROR -  can't create " << freq_data_file_name << endl;
		exit(1);
	}

	t_stream.seekg(0, ios::end);
	long long Length_of_together_file = t_stream.tellg();
	t_stream.close();

	ifstream datastream(freq_data_file_name.c_str(), ios::binary);
	if (!datastream)
	{
		log_stream << "ERROR -  can't create " << freq_data_file_name << endl;
		cout << "ERROR -  can't create " << freq_data_file_name << endl;
		exit(1);
	}


	zakroma_ = new char[Length_of_together_file];


	datastream.read((char*)zakroma_, Length_of_together_file*sizeof(char));


	memcpy(&non_zerro_occurence_number_, zakroma_, sizeof(int));


}


int   Frequency_extrapolation::
get_occurence (	const int record_index	)
{
	

	assert (index_to_pos_); // not defined map yet

	if ( record_index < 0 || record_index >= number_of_elements_ )
		return -1;

	int non_zerro_occurence_record_index = index_to_pos_ [record_index];

	if (non_zerro_occurence_record_index == -1)
		return 0;
	else 
	{
		int shift =  sizeof (int ) +	non_zerro_occurence_record_index*( record_length_ + 2* sizeof (int) )  + 	sizeof (int ) ;

		int phrase_len = sizeof (int);
		int current_occurence;

		if (zakroma_)
			memcpy (&current_occurence,zakroma_ + shift, phrase_len);
		else 
		{
			
//			char *temp	= new char	[ phrase_len ];

			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  &current_occurence,	phrase_len );
			
//			memcpy (&current_occurence,	 temp,	phrase_len );


//			delete [] temp	;

		} 

		return current_occurence;
	}
}

// Ok
vector <double>   Frequency_extrapolation::   //если нету - возвр. нулевой  длины вектор
get_single_distance_to_clusters_sum( const int record_index )
{
	assert (index_to_pos_); // not defined map yet

	vector < double > record ;	record.resize (0); 

	double *array = new double [number_of_classes_];

	if ( record_index < 0 || record_index >= number_of_elements_ ) 		return record ;

	int non_zerro_occurence_record_index = index_to_pos_ [record_index];

	if (non_zerro_occurence_record_index == -1)		return record ;
	else 
	{
		int shift = sizeof (int ) +
					non_zerro_occurence_record_index*( record_length_ + 2* sizeof (int) )  + 
					2*sizeof (int ) ;
		int phrase_len = number_of_classes_*sizeof (double);
		

		if (zakroma_)
			memcpy (array,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  array,	phrase_len );

//			memcpy (array,	 temp,	phrase_len );

//			delete [] temp	;
		}

		record.resize(number_of_classes_); 
		for (int ii=0; ii<number_of_classes_;ii++ ) 		record[ii] = array [ii];

		delete [] array ;

		return record;
	}

}



vector <double>   Frequency_extrapolation::
get_single_squared_distance_to_clusters_sum (
	const int record_index)
{
	assert (index_to_pos_); // not defined map yet

	vector < double > record ;	record.resize (0); 

	double *array = new double [number_of_classes_];

	if ( record_index < 0 || record_index >= number_of_elements_ ) 		return record ;

	int non_zerro_occurence_record_index = index_to_pos_ [record_index];

	if (non_zerro_occurence_record_index == -1)		return record ;
	else 
	{
		int shift = sizeof (int ) +
					non_zerro_occurence_record_index*( record_length_ + 2* sizeof (int) )  + 
					2*sizeof (int ) +
					number_of_classes_*sizeof (double);
		int phrase_len = number_of_classes_*sizeof (double);
	
		if (zakroma_)
			memcpy (array,zakroma_ + shift, phrase_len);
		else 
		{
///			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  array,	phrase_len );
//			memcpy (array,	 temp,	phrase_len );

//			delete [] temp	;

		}

		record.resize(number_of_classes_); 
		for (int ii=0; ii<number_of_classes_;ii++ ) 		record[ii] = array [ii];

		delete [] array ;

		return record;
	}

}

// Ok
vector <double>   Frequency_extrapolation::
get_single_inverse_distance_to_clusters_sum (
	const int record_index )
{
	assert (index_to_pos_); // not defined map yet

	vector < double > record ;	record.resize (0); 

	double *array = new double [number_of_classes_];

	if ( record_index < 0 || record_index >= number_of_elements_ ) 		return record ;

	int non_zerro_occurence_record_index = index_to_pos_ [record_index];

	if (non_zerro_occurence_record_index == -1)		return record ;
	else 
	{
		int shift = sizeof (int ) +
					non_zerro_occurence_record_index*( record_length_ + 2* sizeof (int) )  + 
					2*sizeof (int ) +
					2*number_of_classes_*sizeof (double);
		int phrase_len = number_of_classes_*sizeof (double);
		
//		memcpy (array,zakroma_ + shift, phrase_len);


		if (zakroma_)
			memcpy (array,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  array,	phrase_len );
//			memcpy (array,	 temp,	phrase_len );

//			delete [] temp	;

		}


		record.resize(number_of_classes_); 
		for (int ii=0; ii<number_of_classes_;ii++ ) 		record[ii] = array [ii];

		delete [] array ;

		return record;
	}
}


vector <double>   Frequency_extrapolation::
get_single_inverse_squared_distance_to_clusters_sum (
	const int record_index)
{
	assert (index_to_pos_); // not defined map yet

	vector < double > record ;	record.resize (0); 

	double *array = new double [number_of_classes_];

	if ( record_index < 0 || record_index >= number_of_elements_ ) 		return record ;

	int non_zerro_occurence_record_index = index_to_pos_ [record_index];

	if (non_zerro_occurence_record_index == -1)		return record ;
	else 
	{
		int shift = sizeof (int ) +
					non_zerro_occurence_record_index*( record_length_ + 2* sizeof (int) )  + 
					2*sizeof (int ) +
					3*number_of_classes_*sizeof (double);
		int phrase_len = number_of_classes_*sizeof (double);
		
//		memcpy (array,zakroma_ + shift, phrase_len);

		if (zakroma_)
			memcpy (array,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  array,	phrase_len );
//			memcpy (array,	 temp,	phrase_len );

//			delete [] temp	;

		}


		record.resize(number_of_classes_); 
		for (int ii=0; ii<number_of_classes_;ii++ ) 		record[ii] = array [ii];

		delete [] array ;

		return record;
	}
}

vector <double>   Frequency_extrapolation::	
get_tot_distance_to_clusters_sum ()
{
		double *d_buffer = new double [number_of_classes_];
		assert (index_to_pos_); // not defined map yet  
// Значит что не прочли еще бинарный файл с данными  функцией _suck_up_kruto() или вообще FILL_UP_FREQUENCY_EXTRAPOLATION_STORE 
			        
		int shift = sizeof (int ) +
					non_zerro_occurence_number_* ( record_length_ + 2* sizeof (int) )  + 
				    number_of_elements_ * sizeof (int) ;

		
//		memcpy (			d_buffer,			zakroma_ + shift,			number_of_classes_*sizeof (double) );
		int phrase_len = number_of_classes_*sizeof (double);

		if (zakroma_)
			memcpy (d_buffer,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  d_buffer,	phrase_len );
//			memcpy (d_buffer,	 temp,	phrase_len );

//			delete [] temp	;

		}


		vector <double>   ddd ; 
		ddd .resize (number_of_classes_);

		for (int ii=0; ii<number_of_classes_; ii++)
			ddd [ii] = d_buffer[ii];

		delete [] d_buffer;

		return ddd ;


}
vector <double>   Frequency_extrapolation::	
get_tot_squared_distance_to_clusters_sum()
{
		double *d_buffer = new double [number_of_classes_];
		assert (index_to_pos_); // not defined map yet  
// Значит что не прочли еще бинарный файл с данными  функцией _suck_up_kruto() или вообще FILL_UP_FREQUENCY_EXTRAPOLATION_STORE 
			        
		int shift = sizeof (int ) +
					non_zerro_occurence_number_* ( record_length_ + 2* sizeof (int) )  + 
				    number_of_elements_ * sizeof (int) + 
					+ number_of_classes_*sizeof (double);


//		memcpy (
//			d_buffer,
//			zakroma_ + shift,
//			number_of_classes_*sizeof (double) );

		int phrase_len = number_of_classes_*sizeof (double);

		if (zakroma_)
			memcpy (d_buffer,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  d_buffer,	phrase_len );
//			memcpy (d_buffer,	 temp,	phrase_len );

//			delete [] temp	;

		}



		vector <double>   ddd ; 
		ddd .resize (number_of_classes_);

		for (int ii=0; ii<number_of_classes_; ii++)
			ddd [ii] = d_buffer[ii];

		delete [] d_buffer;

		return ddd ;

}

vector <double>   Frequency_extrapolation::	
get_tot_inverse_distance_to_clusters_sum()
{
			double *d_buffer = new double [number_of_classes_];
		assert (index_to_pos_); // not defined map yet  
// Значит что не прочли еще бинарный файл с данными  функцией _suck_up_kruto() или вообще FILL_UP_FREQUENCY_EXTRAPOLATION_STORE 
			        
		int shift = sizeof (int ) +
					non_zerro_occurence_number_* ( record_length_ + 2* sizeof (int) )  + 
				    number_of_elements_ * sizeof (int) +
					+ 2*number_of_classes_*sizeof (double);


//		memcpy (
//			d_buffer,
//			zakroma_ + shift,
//			number_of_classes_*sizeof (double) );


		int phrase_len = number_of_classes_*sizeof (double);

		if (zakroma_)
			memcpy (d_buffer,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  d_buffer,	phrase_len );
//			memcpy (d_buffer,	 temp,	phrase_len );

//			delete [] temp	;

		}


		vector <double>   ddd ; 
		ddd .resize (number_of_classes_);

		for (int ii=0; ii<number_of_classes_; ii++)
			ddd [ii] = d_buffer[ii];

		delete [] d_buffer;

		return ddd ;


}

vector <double>   Frequency_extrapolation::	
get_tot_inverse_squared_distance_to_clusters_sum ()
{
			double *d_buffer = new double [number_of_classes_];
		assert (index_to_pos_); // not defined map yet  
// Значит что не прочли еще бинарный файл с данными  функцией _suck_up_kruto() или вообще FILL_UP_FREQUENCY_EXTRAPOLATION_STORE 
			        
		int shift = sizeof (int ) +
					non_zerro_occurence_number_* ( record_length_ + 2* sizeof (int) )  + 
				    number_of_elements_ * sizeof (int) +
					+ 3*number_of_classes_*sizeof (double);

//		memcpy (
//			d_buffer,
//			zakroma_ + shift,
//			number_of_classes_*sizeof (double) );

		int phrase_len = number_of_classes_*sizeof (double);

		if (zakroma_)
			memcpy (d_buffer,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  d_buffer,	phrase_len );
//			memcpy (d_buffer,	 temp,	phrase_len );

//			delete [] temp	;
		}

		vector <double>   ddd ; 
		ddd .resize (number_of_classes_);

		for (int ii=0; ii<number_of_classes_; ii++)
			ddd [ii] = d_buffer[ii];

		delete [] d_buffer;
		return ddd ;
}

int		Frequency_extrapolation::	
get_total_sample_size()
{
		double *d_buffer = new double [number_of_classes_];
		assert (index_to_pos_); // not defined map yet  
// Значит что не прочли еще бинарный файл с данными  функцией _suck_up_kruto() или вообще FILL_UP_FREQUENCY_EXTRAPOLATION_STORE 
			        
		int shift = sizeof (int ) +
					non_zerro_occurence_number_* ( record_length_ + 2* sizeof (int) )  + 
				    number_of_elements_ * sizeof (int) +
					+ 4*number_of_classes_*sizeof (double);


		int		total_sample_size; 
//		memcpy (
//			&total_sample_size,
//			zakroma_ + shift,
//			sizeof (int));

		int phrase_len = sizeof (int);

		if (zakroma_)
			memcpy (&total_sample_size,zakroma_ + shift, phrase_len);
		else 
		{
//			char *temp	= new char	[ phrase_len ];
			standby_datastream_.seekg(shift,ios::beg);
			standby_datastream_.read ( (char* )  &total_sample_size,	phrase_len );
//			memcpy (&total_sample_size,	 temp,	phrase_len );

//			delete [] temp	;

		}

		return total_sample_size;

}

Frequency_chain_constants  Frequency_extrapolation::
get_single_frequency_chain_constants ( 
	const string & sequence, 
	ifstream & data_stream ) 
{
		vector < int > processed_index  = translate_sequence_to_degenerate_array(sequence) ;
		Frequency_chain_constants fcc ( processed_index,data_stream, this );
		return fcc;
}

vector <vector <double> >  Frequency_extrapolation::
  calc_cluster_mutual_distance ()
{
	vector <vector <double> > distance_mutual;
	distance_mutual.resize(number_of_classes_);
	for (int ii=0;ii<number_of_classes_;ii++)
		distance_mutual[ii].resize(number_of_classes_);

	double *w	= new double [fragment_length_ *9];

	for (int kk=0;kk<(fragment_length_ *9);kk++)
		w[kk] = 1;
	
	int ii;
	for ( ii=0;ii<number_of_classes_;ii++)
	{
		for (int jj=0;jj<number_of_classes_;jj++)
		{

			bool error_flag= true;
			double distance  = kabsch_rmsd	(
									claster_motif_coordinates_[ii], 
									claster_motif_coordinates_[jj], 
									w, 
									(3*fragment_length_),
									error_flag);


			distance_mutual[ii][jj] = distance ;
		}
	}

	delete [] w;

	return distance_mutual;
}


//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
//  ATTENSION: Frequency_extrapolation_ZIP !!!!!!!!!!
