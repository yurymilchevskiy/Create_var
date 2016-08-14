#include "ELM_Test.h"

#include <fstream>

#include "save_as_fasta.h"

#include "some_usefull_util_for_elm_task.h"

#include "../CowardVariables/CowardVariables.h"
#include "../CommonFunc.h"

using namespace std;

extern ofstream log_stream;
//extern Censorship configuration;

ELM_Test::
~ELM_Test()
{
	cout << "ELM_Test PASSED: " << "  failed: " << get_failed() << "	passed: " << get_passed() << endl;
}

void ELM_Test::dummy_preparing_ELM_data_for_DA()
{
	string path_to_motif_store = "D:/Didona/Store/ELM/MOD_N-GLC_1/motif/";
	string path_to_control_store = "D:/Didona/Store/ELM/MOD_N-GLC_1/control/";



	string path_to_motif_StartPointsMap		= "D:/Didona/Store/ELM/MOD_N-GLC_1/motif/" + string("StartPointsMap");
	string path_to_control_StartPointsMap	= "D:/Didona/Store/ELM/MOD_N-GLC_1/control/" + string("StartPointsMap");


//	string file_name = path_to_motif_store;
//	ifstream  in_motif(file_name.c_str());

	map < string, vector<int> >motif_StartPointsMap; 
	map < string, vector<int> >control_StartPointsMap;

	int total_motiff_number =
		fill_up_StartPointsMap(
			path_to_motif_StartPointsMap,
			motif_StartPointsMap);
	

	int total_control_number =
		fill_up_StartPointsMap(
			path_to_control_StartPointsMap,
			control_StartPointsMap);

	string path_to_cowardvariables_task_file = "D:/Didona/Store/ELM/MOD_N-GLC_1/model4/CowardVariables.task";
	CowardVariables *cowa_creator_ = new CowardVariables(path_to_cowardvariables_task_file);

	typedef   map < string, vector <int> > MAP_SEQUENCEMAP_TO_VECTOR_INT;
	MAP_SEQUENCEMAP_TO_VECTOR_INT::iterator theMAP_SEQUENCEMAP_TO_VECTOR_INT;

	vector < vector <double> > sophisticated_variables_motif =
		single_set_for_ELM_data_for_DA(
			path_to_motif_store,
			motif_StartPointsMap,
			cowa_creator_);

	
	vector < vector <double> > sophisticated_variables_control =
		single_set_for_ELM_data_for_DA(
			path_to_control_store,
			control_StartPointsMap,
			cowa_creator_);

	int number_of_variables = cowa_creator_->get_number_of_variables();

	delete cowa_creator_; // больше не нужен генератор переменных


	string path_to_ready_DA_data = "D:/Didona/Store/ELM/MOD_N-GLC_1/model4/DA_data.txt";

	ofstream out(path_to_ready_DA_data.c_str() );
	if (!out) {
		log_stream	<< "ERROR -  can't create file" << path_to_ready_DA_data << endl;
		cout		<< "ERROR -  can't create file" << path_to_ready_DA_data << endl;
		throw "ERROR - can't create file for DA_data";
	}


	int number_of_classes = 2;
	int number_of_motif_cases	= sophisticated_variables_motif.size();
	int number_of_control_cases = sophisticated_variables_control.size();
	int number_of_cases = number_of_motif_cases + number_of_control_cases;

	PutVa(number_of_classes,	out, 5,		3, 'r');
	PutVa(number_of_cases,		out, 10,	3, 'r');
	PutVa(number_of_variables,	out, 5,		3, 'r');
	out << endl;
	
	
	int current_class_index = 0;
	for (int ii = 0; ii < number_of_motif_cases; ii++)
	{
		PutVa(current_class_index, out, 5, 3, 'r');
		for (int kk = 0; kk < number_of_variables;kk++)
			PutVaDouble(sophisticated_variables_motif[ii][kk], out, 12, 6, 'r');

		out << endl;
	}
	 current_class_index = 1;
	for (int ii = 0; ii < number_of_control_cases; ii++)
	{
		PutVa(current_class_index, out, 5, 3, 'r');
		for (int kk = 0; kk < number_of_variables; kk++)
			PutVaDouble(sophisticated_variables_control[ii][kk], out, 12, 6, 'r');

		out << endl;
	}
	

}
