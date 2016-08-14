#pragma warning( disable : 4786 )

#include "Fill_up_accepted_chain_data_test.h"

#include "accepted_chain_data.h"

#include "CommonFunc.h"

#include "Censorship.h"

#include <fstream>

using namespace std;

extern ofstream log_stream;
extern Censorship configuration;

Fill_up_accepted_chain_data_test::
~Fill_up_accepted_chain_data_test()
{
	cout << "Fill_up_accepted_chain_data_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}

void Fill_up_accepted_chain_data_test::
first_simple_test ()
{

	vector < string >   accepted_chain_ID_list;
	vector < int >      accepted_chain_lenth ;

	string binary_file_name = 
		configuration.option_meaning("Path_to_Chain_store") + 	
		string("accepted_chain_list.bin");



	fill_up_accepted_chain_data ( 
		accepted_chain_ID_list, 
		accepted_chain_lenth,
		binary_file_name );


	string full_path_to_test = "TEST/Fill_up_accepted_chain_data_test";

	ofstream out_stream ( full_path_to_test.c_str() );
	if ( ! out_stream ) 	
	{
		cout		<< " Can't create file " << full_path_to_test	<< endl;
		log_stream	<< " Can't create file " << full_path_to_test<< endl;
		exit (1);
	}

	for (int ii=0;ii<accepted_chain_lenth.size();ii++)
		out_stream << accepted_chain_ID_list[ii] << "  " << accepted_chain_lenth[ii] << endl;

}