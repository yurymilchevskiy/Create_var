#ifndef FILL_UP_ACCEPTED_CHAIN_DATA_H 
#define FILL_UP_ACCEPTED_CHAIN_DATA_H 

#include <vector>
#include <string>

using namespace std; 

void fill_up_accepted_chain_data ( 
	vector < string >  & accepted_chain_ID_list, 
	vector < int >     & accepted_chain_lenth,
	const string path_to_accepted_chain_list_file);

#endif 