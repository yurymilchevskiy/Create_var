#include "some_usefull_util_for_elm_task.h"
#include "save_as_fasta.h"

vector < vector <double> > single_set_for_ELM_data_for_DA(
	const string & path_to_motif_store	,
	map < string, vector<int> > & StartPointsMap,
	CowardVariables *cowa_creator_)
{

	typedef   map < string, vector <int> > MAP_SEQUENCEMAP_TO_VECTOR_INT;
	MAP_SEQUENCEMAP_TO_VECTOR_INT::iterator theIterator;


	vector <vector <double> > da_set;
	int counter;
	for (theIterator = StartPointsMap.begin(); theIterator != StartPointsMap.end(); theIterator++)
	{
	string key = (*theIterator).first;
		vector <int> indexes = (*theIterator).second;

		string path_to_fasta_file = path_to_motif_store + key + string(".fasta");

		string fasta_head, sequence;
		read_fasta(
			path_to_fasta_file,
			fasta_head,
			sequence);

		vector < vector < double > > sophisticated_variables;
		cowa_creator_->process_chain(sequence);
		sophisticated_variables = cowa_creator_->get_sophisticated_variables();

	//	int shift = 2; // AttENSION. 
		for (int ii = 0; ii < indexes.size(); ii++)
		{
			da_set.push_back(sophisticated_variables[indexes[ii]]);
		}
		counter++;
		if (counter % 100 == 0)
			cout << counter << endl;
	}
	return da_set;
};


