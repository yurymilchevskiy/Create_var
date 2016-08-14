#include <fstream>
#include <iostream>
#include <string>

#include "save_as_fasta.h"

using namespace std;
extern ofstream log_stream;
 
void save_as_fasta(
	const string & Path_to_file,
	const string & header,
	const string & sequence)
{
	ofstream out(Path_to_file.c_str());
	if (!out) {
		cout << "ERROR -  can't create " << Path_to_file << endl;
		throw "save_as_fasta ERROR";
	} 

	out << header << endl;
	int seq_len = sequence.size();
	int counter = 0;
	for (int ii = 0; ii < seq_len; ii++)
	{
		counter++;
		out << sequence[ii];
		if (counter % 60 == 0)
		{
			out << endl;
			counter = 0;
		}
	}
}