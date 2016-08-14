#pragma warning( disable : 4786 )

#include "../Cluster_set/Cluster_set.h"

#include "../Censorship.h" 
#include "../Sheduler.h"
#include "../CommonFunc.h"
#include <iostream>

#include <cassert>

extern Censorship configuration;   
extern ofstream log_stream;

using namespace std;



vector <int > pull_out_claster_origin_structure_list (
	const string & exact_path_to)
{

	ifstream source_stream ( exact_path_to .c_str() );
	if ( ! source_stream ) 	
	{
		cout		<<	exact_path_to << " not found " << endl;
		log_stream	<<	exact_path_to << " not found " << endl;
		exit (-1);
	}


	vector < int > choosen_indexes;


	int  counter =0 ;
	string current_line;
	while ( getline(source_stream,current_line,'\n' ) )
	{

		counter ++ ;
		string word;
		istringstream ist( current_line );

		ist >> word ;
		ist >> word ;
		ist >> word ;
		ist >> word ;

		int number_of_classes ;

		ist >> number_of_classes  ;

	//	assert ( check_number == waited_cluster_number );
	//	if ( check_number != waited_cluster_number )
	//	{
	//		cout		<<	"wrong waited cluster  in clusterization protocol file" << endl;
	//		log_stream	<<	"wrong waited cluster  in clusterization protocol file" << endl;
	//		exit (-2);
	//	}

  //		choosen_indexes.resize (waited_cluster_number );

		choosen_indexes.resize (number_of_classes);

		for (int ii=0; ii< number_of_classes; ii++ ) 
			ist >> choosen_indexes[ii];

	}


	return choosen_indexes;
}