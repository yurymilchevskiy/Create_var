#pragma warning( disable : 4786 )

#include <sstream>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;


extern ofstream log_stream;
//extern Censorship configuration;

int fill_up_StartPointsMap(
	const string & path_to_store,
	map < string, vector<int> > & StartPointsMap)
{
	ifstream  in(path_to_store.c_str());
	if (!in)	
	{
		log_stream << "CowardVariables_test: ERROR -  can't create file" << endl;
		cout       << "CowardVariables_test: ERROR -  can't create file" << endl;
		throw "CowardVariables_test: ERROR -  can't create file";
		return 0;
	}


	string current_line;
	while (getline(in,current_line, '\n'))
	{
		if (current_line.size() == 0)
			continue;
		if (current_line[0] == '/' ||
			current_line[0] == '#' || 
			current_line[0] == ' ' ||
			current_line[0] == '\n' ||
			current_line[0] == '\0')
			continue;

		vector <int> index_set;
		string key, meaning;
		
		istringstream ist(current_line);
		{
			ist >> key;
			
			int current_index;
			while (ist >> current_index)
				index_set.push_back(current_index);
		
			StartPointsMap[key] = index_set;
		}
		index_set.resize(0);
	}

	return StartPointsMap.size();


}