#include < vector >

using namespace std;

void fill_up_appropriate_crude_index_subset (
	vector < vector < int > >   & index_of_degeneration,
	const int position,
	const int lb,
	const int rb, 
	const vector < int > & crude_index_array,
	const int virtual_residue_index, 
	vector < int > & index_subset );