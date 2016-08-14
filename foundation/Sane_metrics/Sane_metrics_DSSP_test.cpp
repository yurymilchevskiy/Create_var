#pragma warning( disable : 4786 )

#include "Sane_metrics_DSSP_test.h"
#include "Sane_metrics_DSSP.h"

#include "../Fragment_base/Fragment_base_subtle.h"
#include "../by_Qmol/kabsch_stolen.h"
#include "../Fragment_base/Chain_binary.h" 
#include "../Cluster_set/Cluster_set.h" 
#include "../CommonFunc.h"

#include "../BioPolymerMechanics/foundation/Model.h"
#include "../BioPolymerMechanics/foundation/Core_iterator.h"
#include "../BioPolymerMechanics/foundation/Atom.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <cassert>

using namespace std;

extern ofstream log_stream;

Sane_metrics_DSSP_test::
~Sane_metrics_DSSP_test()
{
	cout << "Sane_metrics_DSSP_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}

void Sane_metrics_DSSP_test::
first_constructor_test ()
{
		Sane_metrics_DSSP   sa_me(  
				 string("DSSP_SET"),
				 string("HHHHH"),
		         SANE_METRICS_DSSP_FILL_UP) ;

}

void Sane_metrics_DSSP_test::
fill_up_subset_by_whole_base_test ()
{
	cout << "fill_up_subset_by_whole_base_test ()" << endl;
		Sane_metrics_DSSP   sa_me(  
				 string("DSSP_SET"),
				 string("LARGE_DSSP_CLUSTERS_INDEX"),
		         SANE_METRICS_DSSP_FILL_UP_SUBSET_BY_WHOLE_BASE) ;

}




void Sane_metrics_DSSP_test::
create_all_existing_DSSP_words_metrics ()
{

	//vector <string> dssp_words = get_all_existing_DSSP_words();

	string path_to_list_file = "D:/Didona/Store/Fragment_base_DSSP_bridge/Second/protocol/list";


	ifstream in ( path_to_list_file .c_str() );
	if ( ! in ) 	
	{
		cout		<<	path_to_list_file << " can't find " << endl;
		log_stream	<<	path_to_list_file << " can't find " << endl;
		throw;
	}

	string  current_word ;
	vector <string> dssp_words;
	while (in >> current_word ) 
		dssp_words.push_back(current_word);
	
	
	int size= dssp_words.size();

// fix 
	for (int ii=3331; ii<size; ii++)
	{
		Sane_metrics_DSSP   sa_me(  
				 string("DSSP_SET"),
				 dssp_words[ii],
		         SANE_METRICS_DSSP_FILL_UP) ;
	}
}

void Sane_metrics_DSSP_test::
show_values_for_dssp_words ()
{

	Sane_metrics_DSSP   sa_me(  
			string("DSSP_SET"),
			string("LARGE_DSSP_CLUSTERS_INDEX"),
		    SANE_METRICS_COMMON_DSSP_USAGE) ;


	int number_of_records = sa_me.get_number_of_records();
	
	double *metrics_sqrt = new double [ number_of_records *sizeof (double) ];
	sa_me.fill_metrics( 0,metrics_sqrt );


	double *metrics_1= new double [ number_of_records *sizeof (double) ];
	sa_me.fill_metrics( 1,metrics_1);

	double *metrics_2= new double [ number_of_records *sizeof (double) ];
	sa_me.fill_metrics( 2,metrics_2);

	double *metrics_3= new double [ number_of_records *sizeof (double) ];
	sa_me.fill_metrics( 3,metrics_3);

	double *metrics_4= new double [ number_of_records *sizeof (double) ];
	sa_me.fill_metrics( 4,metrics_4);

	vector <int> base_indexes_for_dssp_words  = sa_me.get_base_indexes_for_dssp_words () ;

	string path_to_list_file = "D:/Didona/Store/Fragment_base_DSSP_bridge/Second/protocol/LARGE_DSSP_CLUSTERS_INDEX.and_names";
	ifstream in ( path_to_list_file .c_str() );
	if ( ! in ) 	
	{
		cout		<<	path_to_list_file << " can't find " << endl;
		log_stream	<<	path_to_list_file << " can't find " << endl;
		throw;
	}


	string path_to_result_file = "D:/Didona/Store/Fragment_base_DSSP_bridge/Second/protocol/LARGE_DSSP_CLUSTERS_INDEX_compare_metrics.txt";
	ofstream out ( path_to_result_file .c_str() );
	if ( ! out ) 	
	{
		cout		<<	path_to_result_file << " can't find " << endl;
		log_stream	<<	path_to_result_file << " can't find " << endl;
		throw;
	}


	string word;
	int index_in_base;

	map <int, string >  index_in_base_to_dssp_word;

	while (in >> index_in_base >> word)
	{
		index_in_base_to_dssp_word [index_in_base] = word;
	}

	
	assert (number_of_records == index_in_base_to_dssp_word.size() );

	for (int ii=0;ii<number_of_records;ii++)
	{
		int current_index = base_indexes_for_dssp_words [ii];
		string dssp_word = index_in_base_to_dssp_word [current_index];

		PutVa (dssp_word,out,8,1,'l');				out << '\t';
		PutVa (current_index,out,10,1,'l');			out << '\t';

		PutVaDouble (metrics_sqrt[ii],out,10,3,'l');out << '\t';
		PutVaDouble (metrics_1[ii],out,10,3,'l');	out << '\t';
		PutVaDouble (metrics_2[ii],out,10,3,'l');	out << '\t';
		PutVaDouble (metrics_3[ii],out,10,3,'l');	out << '\t';
		PutVaDouble (metrics_4[ii],out,10,3,'l');	
		out << endl;

	}

}


