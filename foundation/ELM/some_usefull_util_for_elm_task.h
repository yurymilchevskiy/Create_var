#ifndef SOME_USEFULL_UTIL_FOR_ELM_TASK_H
#define SOME_USEFULL_UTIL_FOR_ELM_TASK_H

#include <sstream>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

#include "../CowardVariables/CowardVariables.h"

using namespace std;

// READ SEQUENCE ID (ANE SETTER BY USER) AND SET OF POSITION IN SEQUENCE CONTAINING REGULAR EXPRESSION STARTS
int fill_up_StartPointsMap(
	const string & path_to_store,
	map < string, vector<int> > & StartPointsMap);

// Return ready variables set for futher DA calculations. Variables setted by CowardVariables.task
vector < vector <double> > single_set_for_ELM_data_for_DA(
	const string & path_to_motif_store,
	map < string, vector<int> > & StartPointsMap,
	CowardVariables *cowa_creator_);


#endif
