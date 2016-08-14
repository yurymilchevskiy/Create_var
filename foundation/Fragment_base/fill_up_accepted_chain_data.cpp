#pragma warning( disable : 4786 )

#include "../CommonFunc.h"
#include "../Censorship.h"

#include "accepted_chain_data.h"

#include <fstream>

using namespace std;

extern ofstream log_stream;
extern Censorship configuration;

void fill_up_accepted_chain_data ( 
	vector < string >  & accepted_chain_ID_list, 
	vector < int >     & accepted_chain_lenth,
	const string path_to_accepted_chain_list_file)
{
	string binary_file_name = 
		configuration.option_meaning("Path_to_Chain_store") + 	
		path_to_accepted_chain_list_file;

	ifstream binary_stream ( binary_file_name.c_str(),ios::binary );

	if ( ! binary_stream )	{	
		log_stream	 << "ERROR -  can't find binary file" << binary_file_name<< endl;
		cout		 << "ERROR -  can't find binary file" << binary_file_name<< endl;
		exit (1);	
	}

	binary_stream.seekg(0,ios::end);
//	int file_length = binary_stream.tellg();
	streamoff file_length = binary_stream.tellg();


	char * data_from_bin = new char [(int)file_length];
	binary_stream.seekg(0,ios::beg);
	binary_stream.read( (char* ) data_from_bin, (int)file_length*sizeof (char) );		

	int residue_number = (int)file_length / 9;

	accepted_chain_ID_list.	resize(residue_number);
	accepted_chain_lenth.	resize(residue_number);

	char * current_pdb_chain_ID  = new char [6];
	memset (current_pdb_chain_ID,0,6 );

	for ( int ii=0; ii < residue_number; ii++)
	{
		int	  current_length;

		memcpy (current_pdb_chain_ID,	data_from_bin + ii*9,		5*sizeof(char) );
		memcpy (&current_length,		data_from_bin + ii*9 + 5,	sizeof(int) );

		accepted_chain_ID_list[ii]	= string ( current_pdb_chain_ID );
		accepted_chain_lenth[ii]	= current_length;

	}

	delete [] data_from_bin;
	delete [] current_pdb_chain_ID;
}


