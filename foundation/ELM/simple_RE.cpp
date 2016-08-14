#include  "../Censorship.h"
#include  <iostream>
#include <sstream>

#include "simple_RE.h"

extern Censorship configuration;
extern ofstream log_stream;

simple_RE::simple_RE( const string & source_file )
{
	ifstream in( source_file.c_str() );

	if (!in) {

		log_stream	<< "simple_RE ERROR: can't file" << source_file << endl;
		cout		<< "simple_RE ERROR: can't file" << source_file << endl;
		throw;
	}

	string current_line;
	while (getline(in, current_line, '\n'))
	{
		if (current_line.size() == 0)
			continue;
		if (current_line[0] == '/' ||
			current_line[0] == '#' ||
			current_line[0] == ' ' ||
			current_line[0] == '\n'||
			current_line[0] == '\0')
			continue;

		int postion;
		bool This_or_All_but;
		string word;
		string buffer;

		{
			istringstream ist(current_line);
			ist >> postion >> buffer >> word;
			if (buffer == "true")		This_or_All_but = true;
			else if (buffer == "false")	This_or_All_but = false;
			else throw "error: WRONG WORD CHECK false or true";
		}
		single_RE_record sre_o(postion, This_or_All_but, word);
		record_.push_back(sre_o);
	}

}

vector <int> simple_RE::find_motifs_in_seq( const string &sequence )
{

	vector <int> motiff_start_positions;
	int seglen = sequence.size();
	int size_SR = 3;
	size_SR = record_.size();
	
	bool coinsidence = false;
	int ii=0, jjj=0, kk=0;
	for (ii = 0; ii < seglen - size_SR + 1; ii++)
	{
		
		for (jjj = 0; jjj < size_SR; jjj++)
		{

			for (kk = 0; kk < record_[jjj].word_.size();kk++ )
			{
				char ch1 = sequence[ii + jjj];
				char ch2 = record_[jjj].word_[kk];

				if (record_[jjj].This_or_All_but_)
				{
					if ( ch1 != ch2 )
					{
						coinsidence = false;
						continue;
					}
					else
					{
						coinsidence = true;
						break;
					}
				}
				else
				{
					if (ch1 == ch2 )
					{
						coinsidence = false;
						continue;
					}
					else
					{
						coinsidence = true;
						break;
					}
				}
			}
			if (!coinsidence)				
				break;
		}
		if (coinsidence)
			motiff_start_positions.push_back(ii);
	
	}
	return motiff_start_positions;
}