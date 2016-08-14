#include "save_as_fasta.h"
#include <iostream>
#include <fstream>

void read_fasta(
	const string & path_to_fasta_file,
	string & fasta_head,
	string & sequence)
{
	ifstream in(path_to_fasta_file.c_str());
	if (!in) {
		cout << "ERROR -  can't find fasta file" << path_to_fasta_file << endl;
		throw "ERROR - can't create file for DA_data";
	}

	getline(in, fasta_head, '\n');

	char ch;
	while (in >> ch)
		sequence += ch;
}
