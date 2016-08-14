#pragma warning( disable : 4786 )
#pragma warning( disable : 4018 )

#include "Cluster_set_DSSP_test.h"
#include "Cluster_set_DSSP.h"

#include "Single_cluster_record.h"

#include "../Fragment_base/Chain_binary.h"
#include "../Fragment_base/Fragment_base_subtle.h" 

#include "../CommonFunc.h"
#include "../Censorship.h"

#include "../Geometry_util/Geometry_util.h"

#include <fstream>
#include <iostream>

#include "../Pair_int_double.h"

#include "../Chain_store/DSSP_binary.h"

using namespace std;

extern ofstream log_stream;
extern Censorship configuration;

/*
vector <vector <int> >  
analyase_local_structure_presence (
	vector < vector < double > >  & coord_in_cluster_system, 
	const double constrain,
	const int number_of_classes,
	const int fragment_length);

void print_local_result(
	string cluster_set_name,
	string &PDB_chain_ID,
	vector < vector < double > >  & coord_in_cluster_system,
	string & sequence, 
	vector <vector <int> > local_presence,
	const int fragment_length,
	const double constrain );
*/

Cluster_set_DSSP_test::
~Cluster_set_DSSP_test()
{
	cout << "Cluster_set_DSSP_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}

void Cluster_set_DSSP_test::
peretrah_subtle_clasters_show ()
{
	int flank_size = 5;

	string path_to_subtle_clasters_show = string ("D:/Didona/Store/Cluster_set/DSSP_SET_30/subtle_clasters_show");
	ifstream  in_stream ( path_to_subtle_clasters_show .c_str());
	if ( ! in_stream)	
	{	
		cout		<< "Netu " << path_to_subtle_clasters_show << endl;
		log_stream	<< "Netu " << path_to_subtle_clasters_show << endl;
		throw ;
	}

	string path_to_result = string ("D:/Didona/Store/Cluster_set/DSSP_SET_30/subtle_clasters_show.dssp_interdependence");
	ofstream  out_stream ( path_to_result.c_str());
	if ( ! out_stream )	
	{	
		cout		<< "Can't create " << path_to_subtle_clasters_show << endl;
		log_stream	<< "Can't create " << path_to_subtle_clasters_show << endl;
		throw ;
	}




	string current_line;
	while( getline( in_stream , current_line, '\n' ) )
	{
		if (current_line.size() == 0) 
			continue;
		if (   current_line[0] == '/'  || 
			   current_line[0] == '#'  || 
			   current_line[0] == ' '  || 
			   current_line[0] == '\n' || 
			   current_line[0] == '\0')
			continue;
//CLASTER 1        0.00000  388498    2PYXA ALIEW     365       365        |    -42.479   178.983   -68.973   -39.421   174.311   -56.841   -50.497  -179.591   -60.686   -44.436  -178.304   -61.636

		string temp; 
		string cluster_ID ; 
		double rmds;
		int index_in_base;
		string pdb_chain_ID;
		string fragment_sequence;
		int serial_index;

		istringstream ist (current_line);
		ist >> temp; 
		ist >> cluster_ID ; 
		ist >> rmds;
		ist >> index_in_base;
		ist >> pdb_chain_ID;
		ist >> fragment_sequence;
		ist >> serial_index;

		DSSP_binary pdb_ID_dssp (pdb_chain_ID,COMMON_USAGE);

		string extended_DSSP_sequence	= pdb_ID_dssp.get_extended_DSSP_sequence();
		string sequence_by_DSSP			= pdb_ID_dssp.get_sequence();

		string extended_DSSP_sequence_fragment = extended_DSSP_sequence.substr(serial_index,5);
		string sequence_by_DSSP_fragment =		 sequence_by_DSSP.substr(serial_index,5);

		string left_flank_dssp;			left_flank_dssp.resize(0);
		string right_flank_dssp;		right_flank_dssp.resize(0);

		string left_flank_sequence;		left_flank_sequence.resize(0);
		string right_flank_sequence;	right_flank_sequence.resize(0);
/*
		if (serial_index - flank_size > 0 )
		{
			left_flank_dssp = extended_DSSP_sequence.	substr(serial_index - flank_size ,flank_size );
			left_flank_sequence = sequence_by_DSSP.		substr(serial_index- flank_size ,flank_size );
		}
		else 
		{
			for (int ii=serial_index - flank_size ;ii<flank_size ;ii++)
			{
				if ( ii < 0 )
				{
					left_flank_dssp += "*";
					left_flank_sequence +="*";
				}
				else
				{
					left_flank_dssp		+= extended_DSSP_sequence[ii];
					left_flank_sequence +=	sequence_by_DSSP[ii];
				}
			}

		}
*/

		out_stream << cluster_ID << ": ";
		PutVaDouble(rmds,out_stream, 8,5,'l');
		PutVa(index_in_base,out_stream, 8,5,'l');

		out_stream << pdb_chain_ID		<< "\t";
		out_stream << fragment_sequence << endl;




		out_stream << "\t\t\t\t";
		out_stream << sequence_by_DSSP_fragment		<< "\t";
		out_stream << extended_DSSP_sequence_fragment << endl;
	
	}

}


void Cluster_set_DSSP_test::
make_claster_origin_list_and_names_test()
{
	 string path_to_list_protocol_files = "D:/Didona/Store/Cluster_set/DSSP_SET/protocol/LIST_PROTOCOL_FILES"; 

	 ifstream p_stream ( path_to_list_protocol_files.c_str());

		if ( ! p_stream )	{	
			log_stream	 << "ERROR -  can't find binary file" << path_to_list_protocol_files<< endl;
			cout		 << "ERROR -  can't find binary file" << path_to_list_protocol_files<< endl;
			return;	
		}

		string word;
		vector <string> protocol_files_set;
		while ( p_stream	 >> word)
		{
			protocol_files_set.push_back(word);
		}


		string path_to_protocol_store =  string("D:/Didona/Store/Cluster_set/DSSP_SET/");
		map < int, string >  cl_in_base_nu_to_dssp_word;

		int size = protocol_files_set.size();
		for (int ii=0; ii<size; ii++)
		{

			make_claster_origin_list_and_names (
				path_to_protocol_store,
				protocol_files_set[ii],
				cl_in_base_nu_to_dssp_word );

		}

		string path_to_dssp_indexes = "D:/Didona/Store/Cluster_set/DSSP_SET/protocol/DSSP_INDEXES"; 

		ofstream out ( path_to_dssp_indexes.c_str());

		if ( ! out )	{	
			log_stream	 << "ERROR -  create file" << path_to_dssp_indexes<< endl;
			cout		 << "ERROR -  create file" << path_to_dssp_indexes<< endl;
			return;	
		}

		typedef map < int, string > MAP_INT_STRING;
		MAP_INT_STRING::const_iterator theIterator;


		int counter =0;
		for ( theIterator = cl_in_base_nu_to_dssp_word.begin();theIterator != cl_in_base_nu_to_dssp_word.end();theIterator ++ ) 
		{
			int		key		=	(*theIterator).first ;
			string	value	=   (*theIterator).second ;

			cout << key << "\t" <<  value << endl;
			 out << key << "\t" <<  value << endl;
		}
}

void Cluster_set_DSSP_test::
prepare_indexes_and_names_for_final_clasterization_test ()
{
	string cluster_set_name = "DSSP_SET";
// сделать ручками список файлоd типа  *.protocol. Ну лень писать функцию для этого
	 
}
/**
void ListFiles ();
void Cluster_set_DSSP_test::
List_files_test () 
{
	ListFiles (); 
}
*/
void Cluster_set_DSSP_test::
optimize_clasterization_test ()
{

		string cluster_set_name = "DSSP_SET";

		Cluster_set_DSSP ob(
			cluster_set_name,
			string("LARGE_DSSP_CLUSTERS_INDEX"),
			FILL_UP_MODEL_CLUSTER_SET_DSSP_MODE);

		

	//	ob.optimize_clasterization();

}


void Cluster_set_DSSP_test::
handle_all_dssp_words ()
{

		string cluster_set_name = "DSSP_SET";

		string list_file = string ("D:/Didona/Store/Fragment_base_DSSP_bridge/Second/protocol/list");

		ifstream p_stream ( list_file.c_str());

		if ( ! p_stream )	{	
			log_stream	 << "ERROR -  can't find binary file" << list_file<< endl;
			cout		 << "ERROR -  can't find binary file" << list_file<< endl;
			return;	
		}

		string word;
		vector <string> dssp_word;
		while ( p_stream	 >> word)
		{
			dssp_word.push_back(word);
		}

		int set_size = dssp_word.size();
		for (int ii=1190;ii<set_size;ii++)
		{

			Cluster_set_DSSP ob(
				cluster_set_name,
				dssp_word[ii],
				FILL_UP_MODEL_CLUSTER_SET_DSSP_MODE);

		}

	//	ob.optimize_clasterization();

}
