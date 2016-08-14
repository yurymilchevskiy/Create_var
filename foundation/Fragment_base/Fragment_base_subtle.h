#ifndef FRAGMENT_BASE_SUBTLE_H 
#define FRAGMENT_BASE_SUBTLE_H 

enum Fragment_base_subtle_run_mode
{
	FRAGMENT_BASE_SUBTLE_FILL_UP,
	FRAGMENT_BASE_SUBTLE_COMMON_USAGE
};

#include <vector>
#include <string>
#include <fstream>

class	Sheduler;

using namespace std; 

class Fragment_base_subtle
{	
public: 
	Fragment_base_subtle	( 
		const string & fragment_dir_name,   
		const Fragment_base_subtle_run_mode run_mode);

	Fragment_base_subtle	() {}; 


	~Fragment_base_subtle	();

	int handle_single_chain_for_fragment_base  ( 
		ofstream & fragment_base_stream ,
		const string & current_chain_ID, const int index_in_accepted_chain_base );

	int get_total_fragment_number () const { return total_fragment_number_;} 

	bool get_coord (	const int record_base_index, 
						double *record );

	bool get_chain_index (	const int record_base_index, 
		int *chain_serial_number, 
		int *position_in_chain );

	int get_fragment_length	() const { return fragment_length_;} 

	int get_record_size		() const { return record_size_; }

	bool get_whole_record (
		const int record_base_index,
		char *whole_record);
	
	string get_accepted_chain_ID_by_index(const int index) { 
		return accepted_chain_ID_list_[index];
	}

	vector < string >  get_accepted_chain_ID_list () const { return  accepted_chain_ID_list_ ;} 


	bool fill_up_record_items ( 
		const int index,
		string	& pdb_chain_ID,
		string	& fragment_sequence,
		int		& serial_number,
		string  & pdb_resudue_number);


int get_index_in_accepted_chain_base	(	const int record_base_index);
int get_chain_serial_number				(	const int record_base_index);

protected:
	string fragment_dir_name_;
	int fragment_length_;  // gauged as residue number; That  is length of polypeptide fragment

	int total_fragment_number_;

	Sheduler *sheduler_;

	vector < string >  accepted_chain_ID_list_;  //
	vector <int>       accepted_chain_lenth_;  

	string full_path_to_fragment_base_;

	ifstream in_base_stream_;

	void fill_up ();
	void init ();

	int record_size_;


	string path_to_dihedral_store_;

};

#endif
