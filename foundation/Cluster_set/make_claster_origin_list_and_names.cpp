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

// Имеем наборы базисных структур, полученных для разных dssp words. Эти наборы в файлах  *.protocol. 
// Записаны в последней строке.
// cl_in_base_nu_to_dssp_word связывает номер в базе данных фрагментов с dssp word (типа HHHTT)
// К dssp word добавляется еще номер в списке из файла *.protocol. Станет  HHHTT_3. _3 это четвер
void  make_claster_origin_list_and_names (
	const string & path_to_protocol_store,
	const string & protocol_file,
//	const string & name,
//	const string & extension,
	map < int, string > & cl_in_base_nu_to_dssp_word )
{

	string name = protocol_file.substr(0,5);

	string path_to_file = path_to_protocol_store + protocol_file;
	ifstream source_stream ( path_to_file.c_str() );
	if ( ! source_stream ) 	
	{
		cout		<<	path_to_file << " not found " << endl;
		log_stream	<<	path_to_file << " not found " << endl;
		exit (-1);
	}

	vector < int > choosen_indexes;

	int  counter =0 ;
	string current_line;

	int number_of_classes ;

	while ( getline(source_stream,current_line,'\n' ) )
	{

		counter ++ ;
		string word;
		istringstream ist( current_line );

		ist >> word ;
		ist >> word ;
		ist >> word ;
		ist >> word ;

		ist >> number_of_classes  ;
		choosen_indexes.resize (number_of_classes);

		for (int ii=0; ii< number_of_classes; ii++ ) 
			ist >> choosen_indexes[ii];
	}
	
	for (int ii=0; ii <  number_of_classes ; ii++)
	{
		ostringstream ost ;
		ost << name << "_" << ii;
		string current_dssp_word = ost.str();

		if ( cl_in_base_nu_to_dssp_word .find( choosen_indexes[ii] ) == cl_in_base_nu_to_dssp_word .end()  ) 
			cl_in_base_nu_to_dssp_word [ choosen_indexes[ii] ] = current_dssp_word;
		else 
		{
			log_stream	<< choosen_indexes[ii] << " repeated " << endl;
			cout		<< choosen_indexes[ii] << " repeated " << endl;
		}
	}
}