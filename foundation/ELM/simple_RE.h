#ifndef SIMPLE_RE_H
#define SIMPLE_RE_H

#include <string>
#include <vector>

// SIMPLEST REGULAR EXPRESSION
using namespace std;

class single_RE_record
{
public:
	single_RE_record::single_RE_record(
		int postion,
		bool This_or_All_but,
		string word) :
			postion_(postion),
			This_or_All_but_(This_or_All_but),
			word_(word)
	{

	}

//private:
	int postion_;
	bool This_or_All_but_;
	string word_;
};

class simple_RE
{
public:
	simple_RE( const string & source_file );

	int get_size() const {
		return record_.size();
	}

	vector <int> find_motifs_in_seq(const string &sequence);
	//vector <int> simple_RE:find_motifs_in_seq(const string &sequence)

private:
	vector <single_RE_record> record_;

};
#endif
