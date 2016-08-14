#ifndef SAVE_AS_FASTA_H
#define SAVE_AS_FASTA_H

#include <string>

using namespace std;

void save_as_fasta(
	const string & Path_to_file,
	const string & header,
	const string & sequence);

void read_fasta(
	const string & path_to_fasta_file,
	string & fasta_head,
	string & sequence);

#endif