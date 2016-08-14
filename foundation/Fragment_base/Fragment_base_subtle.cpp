#pragma warning( disable : 4786 )

#include "Fragment_base_subtle.h"
#include "Chain_binary.h"

#include "../CommonFunc.h"
#include "../Censorship.h"

#include "accepted_chain_data.h"

#include "../Sheduler.h"

#include <fstream>

#include "../tri_to_one_and_vice_versa_aa_translation.h"

using namespace std;

extern ofstream log_stream;
extern Censorship configuration;

Fragment_base_subtle::~Fragment_base_subtle()
{ 
	if (!sheduler_ ) 		
		delete sheduler_; 
}

Fragment_base_subtle::
Fragment_base_subtle( 
	const string & fragment_dir_name,
	const Fragment_base_subtle_run_mode run_mode ) :
		 fragment_dir_name_ (fragment_dir_name),
		 sheduler_ (0),	 
		 total_fragment_number_ (0)
{
// надо fragment_length_ всосать откуда-то  *****
//	 fragment_length_ = 5;
// надо fragment_length_ всосать откуда-то  *****

	sheduler_					= new Sheduler  (	
		configuration.option_meaning("Path_to_Fragment_base")  +  
		fragment_dir_name_ + string ("/") + 
		string ("sheduler") ) ;

	fragment_length_ =		atoi ( sheduler_->option_meaning ("FRAGMENT_LENGTH").c_str() );

	string path_to_accepted_chain_list_file =
		//configuration.option_meaning("Path_to_Chain_store") + 
		sheduler_->option_meaning ("ACCEPTED_CHAIN_LIST_FILE_NAME");
	
	fill_up_accepted_chain_data ( 
		accepted_chain_ID_list_, 
		accepted_chain_lenth_,
		path_to_accepted_chain_list_file);

	full_path_to_fragment_base_ = 
		configuration.option_meaning("Path_to_Fragment_base") + 
		fragment_dir_name_ + string ("/fragment_base.bin") ;



	record_size_ = fragment_length_*9*sizeof (double) + 2*sizeof (int);


	if		( run_mode == FRAGMENT_BASE_SUBTLE_FILL_UP )			fill_up ();
	else if ( run_mode == FRAGMENT_BASE_SUBTLE_COMMON_USAGE )		init ();
	else 
	{
		log_stream	<< "Fragment base " <<  fragment_dir_name << "error"  << endl;
		cout		 << "Fragment base " <<  fragment_dir_name << "error"  << endl;
		exit (1);
		
	}
}

void Fragment_base_subtle::
init ()
{
 	in_base_stream_.open ( full_path_to_fragment_base_ .c_str(),ios::binary  );
	if ( ! in_base_stream_  ) 	
	{
		cout		<< " Can't finf file " << full_path_to_fragment_base_  << endl;
		log_stream	<< " Can't find file " << full_path_to_fragment_base_ << endl;
		exit (1);
	}

	
	in_base_stream_.seekg(0,ios::end);
//	int file_length = in_base_stream_.tellg();

	streamoff  file_length = in_base_stream_.tellg();

	total_fragment_number_  = ( (int) file_length  - sizeof (int) ) / record_size_;
}


void Fragment_base_subtle::
fill_up ()
{
	
 	ofstream fragment_base_stream ( full_path_to_fragment_base_ .c_str(),ios::binary   );
	if ( ! fragment_base_stream  ) 	
	{
		cout		<< " Can't create file " << full_path_to_fragment_base_  << endl;
		log_stream	<< " Can't create file " << full_path_to_fragment_base_ << endl;
		exit (1);
	}
// immediatelly write length to ths start position 
//	fragment_base_stream.write( (char* ) & fragment_length_, sizeof (int) );		


	for (unsigned  ii=0;ii<accepted_chain_ID_list_.size();ii++)
	{
		handle_single_chain_for_fragment_base( 
			fragment_base_stream,
			accepted_chain_ID_list_[ii],
			ii);
		if ( ii%100 == 0)
			cout << ii << endl;

	}

}


int Fragment_base_subtle::   // returns number of valid fragments
handle_single_chain_for_fragment_base  ( 
	ofstream & fragment_base_stream ,
	const string & current_chain_ID, const int index_in_accepted_chain_base )
{
	Chain_binary chain( current_chain_ID );

	int				number_of_residues		= chain.get_number_of_residues();
	bool*			is_there_coord			= chain.get_is_there_coord();
	bool*			is_geometry_admissible  = chain.get_is_geometry_admissible();

	double *coord = new double [fragment_length_*9];

	int valid_fragment_number = 0;

	for (int ii=0; ii<number_of_residues - fragment_length_ + 1; ii++)
	{
		bool flag= true;
		for (int kk=0;kk<fragment_length_;kk++)
			flag = flag & is_geometry_admissible [ii+kk];

		if (! flag)
			continue;

		chain.extract_fragment (ii, fragment_length_, coord );
		
		//int pre_pos =fragment_base_stream.tellp();
		streamoff  pre_pos =fragment_base_stream.tellp();


		 fragment_base_stream.write( (char* ) &index_in_accepted_chain_base, sizeof (int) );		
		 fragment_base_stream.write( (char* ) &ii, sizeof (int) );		
	     fragment_base_stream.write( (char*)  coord, fragment_length_*9*sizeof (double) );		

		//int post_pos =fragment_base_stream.tellp();
		 streamoff post_pos =fragment_base_stream.tellp();

		valid_fragment_number++;
	}

	delete [] coord;
	total_fragment_number_ += valid_fragment_number;

	return valid_fragment_number;
}


// 
int Fragment_base_subtle::
get_index_in_accepted_chain_base (
	const int record_base_index)
{

	int index_in_accepted_chain_base;

	if ( record_base_index > total_fragment_number_  ||  record_base_index  < 0 )
		return -1;
	else 
	{
		in_base_stream_.seekg(record_base_index*record_size_ ,ios::beg);
		in_base_stream_.read( (char* ) &index_in_accepted_chain_base, sizeof(int) );		
		return index_in_accepted_chain_base;
	}
}

int Fragment_base_subtle::
get_chain_serial_number(
	const int record_base_index)
{

	int chain_serial_number;

	if ( record_base_index > total_fragment_number_  ||  record_base_index  < 0 )
		return -1;
	else 
	{
		in_base_stream_.seekg(record_base_index*record_size_ + sizeof (int) ,ios::beg);
		in_base_stream_.read( (char* ) &chain_serial_number, sizeof(int) );		
		return chain_serial_number;
	}
}



bool Fragment_base_subtle::
get_coord (		
	const int record_base_index, 
	double *record )
{
	if ( record_base_index > total_fragment_number_  ||  record_base_index  < 0 )
		return false;
	else
	{
		in_base_stream_.seekg(record_base_index*record_size_ + 2*sizeof(int),ios::beg);
		in_base_stream_.read( (char* ) record , record_size_ - 2*sizeof(int) );		
		return true;
	}
}

bool Fragment_base_subtle::
get_chain_index (	const int record_base_index, 
	int *chain_serial_number, 
	int *position_in_chain )
{
	if ( record_base_index > total_fragment_number_  ||  record_base_index  < 0 )
		return false;
	else
	{
		in_base_stream_.seekg(record_base_index*record_size_ ,ios::beg);
		in_base_stream_.read( (char* ) chain_serial_number , sizeof (int) );
		in_base_stream_.read( (char* ) position_in_chain , sizeof (int) );
		return true;
	}
}


bool Fragment_base_subtle::
get_whole_record (
	const int record_base_index,
	char *whole_record)
{
	if ( record_base_index > total_fragment_number_  ||  record_base_index  < 0 )
		return false;
	else
	{
		in_base_stream_.seekg(record_base_index*record_size_ ,ios::beg);
		in_base_stream_.read( (char* ) whole_record, record_size_*sizeof (char) );
		return true;
	}

}




bool Fragment_base_subtle::
fill_up_record_items ( 
	const int record_base_index,
	string	& pdb_chain_ID,
	string	& fragment_sequence,
	int		& serial_number,
	string  & pdb_resudue_number )
{

	int index_in_accepted_chain_base	= get_index_in_accepted_chain_base (record_base_index);
	serial_number						= get_chain_serial_number(record_base_index);

	pdb_chain_ID = accepted_chain_ID_list_ [index_in_accepted_chain_base ];
	int seq_len = accepted_chain_lenth_[index_in_accepted_chain_base ];
	
	Chain_binary cb(pdb_chain_ID);


	vector < string > in_chain_residue_number  = cb.get_in_chain_residue_number();
	pdb_resudue_number = in_chain_residue_number  [serial_number];



	vector < string > residue_names = 	cb.get_residue_names();
	fragment_sequence.resize(0);
	for (int kk=serial_number; kk < serial_number + fragment_length_ ; kk++ )
		fragment_sequence += tri_to_one_letter_aa  ( residue_names [kk]  );


/*

		int postion =  index  *  base_record_length_;
		i_pdbid_pos_stream_.seekg(postion) ;

		int current_position = i_pdbid_pos_stream_.tellg() ;

		if ( current_position == -1 )
			return false;

		char ch_pdb_chain_ID		[100];		memset ( ch_pdb_chain_ID,		0 ,100 );
		char ch_fragment_sequence	[100];		memset ( ch_fragment_sequence,	0 ,100 );
		char ch_serial_number		[100];		memset ( ch_serial_number ,		0 ,100 );
		char ch_pdb_resudue_number	[100];		memset ( ch_pdb_resudue_number,		0 ,100 );
		char record_ch				[1000];		memset ( record_ch,		0 ,1000 );

		i_pdbid_pos_stream_.read( (char* ) record_ch,	base_record_length_* sizeof (char) ); 

		memcpy ( ch_pdb_chain_ID,		record_ch,						5 );
		memcpy ( ch_fragment_sequence,	record_ch + 5,					length_ );
		memcpy ( &serial_number ,		record_ch + 5 + length_  ,		sizeof (int));
		memcpy ( ch_pdb_resudue_number,	record_ch + 5 + length_  +	sizeof(int),	5);


		pdb_chain_ID		= string ( ch_pdb_chain_ID ) ;
		fragment_sequence	= string ( ch_fragment_sequence );
//		serial_number      = atoi	 ( ch_serial_number); 
		pdb_resudue_number = string ( ch_pdb_resudue_number );

		return true;
***/

	return true;
}

