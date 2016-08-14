#ifndef CHAIN_BINARY_H 
#define CHAIN_BINARY_H 

#include <vector>
#include <string>

using namespace std; 

class Chain_binary
{	
public: 
	Chain_binary ( const string & pdb_chain_ID );


	Chain_binary ( // Только для случая когда координаты главной цепи берутся извне. 
		const string & pdb_chain_ID,
		int		number_of_residues_,
		double	*coord_set_);

	~Chain_binary ();

	int		get_number_of_residues	() const { return number_of_residues_; }
	string	get_pdb_chain_ID		() const { return pdb_chain_ID_; }

	vector <string> get_residue_names () ;
	vector <string> get_in_chain_residue_number (); 

	string	get_sequence () ;


	int		*get_serial_index			() ;
	bool    *get_is_there_coord			() ;
	bool	*get_is_geometry_admissible () ;
	double  *get_coord_set				() ;

	void	print_protocol ();

	bool	extract_fragment (const int start_pos, const int length, double *coord );

	void	save_pdb_fragment (const int start_pos,const int length, const string & path_to_pdb);
//	void	save_pdb_fragment (const double *coord,const int length, const string & path_to_pdb);

	void	save_pdb_fragment (
		const double *coord,
		const int length,
		const string & path_to_pdb,
		vector <string> & residue_names,
		vector <string> & in_chain_residue_number);

	void save_pdb_fragment (
		const double *coord,
		const int length, 
		const string & path_to_pdb,
		vector <string> & residue_names,
		vector <string> & in_chain_residue_number,
		char chain_ID);

	void	fragment_to_principal_axes (
		const int start_pos, 
		const int length, 
		double *coord );

	vector < vector < double > >  const & get_set_of_coordinate_in_clasters_system() const
	{
		return set_of_coordinate_in_clasters_system_; 
	}

	void Chain_binary::
		positioning_chain_by_clasters_set(
			double **claster_motif_coordinates,
			const int fragment_length,
			const int number_of_classes,
			vector < vector < double > >  & set_of_coordinate_in_clasters_system);

    // похоже, ВОЗВРАЩАет НЕПРАВИЛЬНОЕ ЗНАЧЕНИЕ!!!!!! Правильное внутри функции
	vector < vector < double > >   & positioning_chain_by_clasters_set (
		double **claster_motif_coordinates,
		const int fragment_length,
		const int number_of_classes ) ;





	vector < vector < double > >  positioning_chain_by_clasters_set ( 
		double **claster_motif_coordinates,
		const int fragment_length,
		const int number_of_classes,
		const int sequence_length) ;  // атавизм - есть тут number_of_residues_ & get_sequence();

	vector < double > get_torsion_angles ();

	int get_number_of_residues ()  { return number_of_residues_;}

	
	void center_of_mass(vector < double > & c_m_coord);

	double* get_coord_set() const { return coord_set_; }

private:
	string	pdb_chain_ID_ ;
	int		number_of_residues_;
//	char	**residue_names_;
//	char	**in_chain_residue_number_;

	int		*serial_index_;
	bool	*is_there_coord_;
	bool	*is_geometry_admissible_;
	double	*coord_set_;

	char	*data_from_bin_;
	void	suck_in_file ();

	vector < vector < double > > set_of_coordinate_in_clasters_system_;
};



vector < vector < double > > 
two_chain_distance_set (
	const string & pdb_chain_ID_1,
	const string & pdb_chain_ID_2,
	const int fragment_length );

bool align_two_chains (
		double *cord_set_1,
		double *cord_set_2,
		int fragment_length);



#endif

