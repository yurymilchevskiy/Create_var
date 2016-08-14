#ifndef	AMINOACID_TO_INDEX_H
#define	AMINOACID_TO_INDEX_H

#pragma warning( disable : 4786 )

#include <map> 

using namespace std;

	map < char, int >  set_index_by_aa() ;
	map < int, char >  set_aa_by_index () ;

	int aminoacid_to_index		( const char aa );
	char index_to_aminoacid		( const int  index );

	int get_virtual_residue_index	(); // returns virtual residue index: for example for residue inside chain
	int get_size_aminoacid_set		() ;
	
	bool is_standard_aa ( const char aa );

#endif