#pragma warning( disable : 4786 )

#include "Chain_binary.h"

#include "../Censorship.h"

#include "../CommonFunc.h"
#include <iostream>
#include <fstream>

#include "../Geometry_util/Fragmrent_tranformation.h"
#include "../by_Qmol/EigenValues3D.h"
#include "../by_Qmol/kabsch_stolen.h"

#include "../tri_to_one_and_vice_versa_aa_translation.h"

#include "../Fragment_base/fill_up_fragment_torsion_angles.h"

#include <cassert>

extern Censorship configuration;
extern ofstream log_stream;

Chain_binary::
~Chain_binary ()
{
	 if (	serial_index_	) 				delete [] serial_index_;

	 if (	is_there_coord_	)				delete [] 
	     is_there_coord_;

	 if (	is_geometry_admissible_	)				 delete [] is_geometry_admissible_;

	 if (	coord_set_		)						 delete [] coord_set_;

	 if (	data_from_bin_	)				 delete [] data_from_bin_;

}
Chain_binary::
Chain_binary ( // Только для случая когда координаты главной цепи берутся извне. 
		const string & pdb_chain_ID,
		int		number_of_residues,
		double	*coord_set): 
 pdb_chain_ID_ (pdb_chain_ID ), 
 number_of_residues_		(number_of_residues),
 serial_index_				(0),
 is_there_coord_			(0),
 is_geometry_admissible_	(0),
 coord_set_					(0),
 data_from_bin_				(0)

{

	is_geometry_admissible_ = (bool *) new bool [number_of_residues_];
	is_there_coord_			= (bool *) new bool [number_of_residues_];
	coord_set_ = (double *) new double [9*number_of_residues_];

	for (int ii=0; ii<number_of_residues_;ii++)
	{
		is_geometry_admissible_[ii] = true;
		is_there_coord_	[ii]		= true;
	}

	memcpy (coord_set_,coord_set, 9*number_of_residues_*sizeof (double) );
}


Chain_binary::
Chain_binary ( const string & pdb_chain_ID ) :
 pdb_chain_ID_ (pdb_chain_ID ), 
 number_of_residues_		(0),
 serial_index_				(0),
 is_there_coord_			(0),
 is_geometry_admissible_	(0),
 coord_set_					(0),
 data_from_bin_				(0)
{

	suck_in_file ();

}
void Chain_binary::
suck_in_file ()
{
	string binary_file_name = 
		configuration.option_meaning("Path_to_Chain_store") + 	
		string ("/binary/") + string (pdb_chain_ID_)  + string (".bin");

	ifstream binary_stream ( binary_file_name.c_str(),ios::binary );

	if ( ! binary_stream )	{	
		log_stream	 << "ERROR -  can't find binary file" << binary_file_name<< endl;
		cout		 << "ERROR -  can't find binary file" << binary_file_name<< endl;
		exit (1);	
	}

	binary_stream.seekg(0,ios::end);
	//int file_length = binary_stream.tellg();
	streamoff file_length = binary_stream.tellg();
	

	number_of_residues_ = ( (int) file_length - 9 ) / 86;

	data_from_bin_ = new char [(int) file_length];
	binary_stream.seekg(0,ios::beg);
	binary_stream.read( (char* ) data_from_bin_, file_length*sizeof (char) );		

	char *name_test = new char [10];
	memset (name_test,0,10);
	memcpy (name_test,data_from_bin_, 5*sizeof (char) );
	string name_test_string (name_test);

	for (int ii=0;ii<5;ii++)
		name_test_string[ii] = toupper (name_test_string[ii]);

	log_stream << " suck_in_file () stage: strange pdb_chain_ID_  :"  << name_test_string  << "  in spite of " << pdb_chain_ID_ << endl;

	assert ( name_test_string == pdb_chain_ID_ );

	delete [] name_test;

}


vector <string> Chain_binary::
get_residue_names () 
{
	int start_pos  = 9;
	vector < string > residue_names;
	residue_names.resize(number_of_residues_);

	char buff[10]; memset (buff,0,10);

	for (int ii=0; ii<number_of_residues_; ii++)
	{
		memcpy (buff,data_from_bin_+start_pos+ii*3,3 );
		residue_names[ii] = ( (string) buff ) ;
	}

	return residue_names;
}

string Chain_binary::
get_sequence () 
{
	vector < string > residue_names = get_residue_names () ;
	string sequence;
	for (int ii=0; ii<residue_names.size(); ii++ ) 
		sequence += tri_to_one_letter_aa (residue_names[ii]);

	return sequence;
}



vector <string> Chain_binary::
get_in_chain_residue_number ()
{
	int start_pos  = 9 + number_of_residues_*3;

	vector < string > in_chain_residue_number;
	in_chain_residue_number.resize(number_of_residues_);

	char buff[6]; memset (buff,0,6);

	for (int ii=0; ii<number_of_residues_; ii++)
	{
		memcpy (buff,data_from_bin_+start_pos+ii*5,4 );
		in_chain_residue_number[ii] = ( (string) buff ) ;

	}

	return in_chain_residue_number;

}

int* Chain_binary::
get_serial_index () 
{
	int start_pos  = 9 + number_of_residues_*3 + number_of_residues_*5;
	
	serial_index_ = (int *) new int [number_of_residues_];

	memcpy (serial_index_,data_from_bin_ +  start_pos,  number_of_residues_*sizeof (int) );

	return serial_index_;
}

bool* Chain_binary::
get_is_there_coord () 
{
	int start_pos  = 9 + number_of_residues_*3 + number_of_residues_*5 + number_of_residues_*4;

	is_there_coord_ = (bool *) new bool [number_of_residues_];

	memcpy (is_there_coord_, data_from_bin_ +  start_pos,  number_of_residues_*sizeof (bool) );

	return  is_there_coord_;
}

bool* Chain_binary::
get_is_geometry_admissible () 
{
	int start_pos  = 9 + number_of_residues_*3 + number_of_residues_*5 + number_of_residues_*4 + number_of_residues_*1;

	is_geometry_admissible_ = (bool *) new bool [number_of_residues_];

	memcpy (is_geometry_admissible_, data_from_bin_ +  start_pos,  number_of_residues_*sizeof (bool) );

	return  is_geometry_admissible_;
}

double* Chain_binary::
get_coord_set ()
{
	int start_pos  = 9 + number_of_residues_*3 + number_of_residues_*5 + number_of_residues_*4 +number_of_residues_*1 + number_of_residues_*1 ;

	coord_set_ = (double *) new double [9*number_of_residues_];

	memcpy (coord_set_, data_from_bin_ +  start_pos,  9*number_of_residues_*sizeof (double) );

	return coord_set_;
}


void Chain_binary::
print_protocol () 
{
	 
	string full_path_to_protocol_file = 
		configuration.option_meaning("Path_to_Chain_store") + 	
		 string ("/protocol/") + string (pdb_chain_ID_)  + string (".binprot");

	ofstream out ( full_path_to_protocol_file.c_str() );
	if ( ! out ) 	
	{
		out			<<	full_path_to_protocol_file  << " can't create " << endl;
		log_stream	<<	full_path_to_protocol_file  << " can't create " << endl;
		return; 
	}


	vector <string> residue_names			= get_residue_names () ;
	vector <string> in_chain_residue_number = get_in_chain_residue_number (); 
	bool*			is_there_coord			= get_is_there_coord() ;
	bool*			is_geometry_admissible  = get_is_geometry_admissible() ;

	for ( int ii=0; ii < number_of_residues_; ii++ )
	{

		PutVa ( residue_names			[ii],	out, 5,4,'l');
		PutVa ( residue_names			[ii],	out, 5,4,'l');
		PutVa ( in_chain_residue_number	[ii],	out, 6,5,'l');

		if ( is_there_coord[ii] ) 	out << "C";
		else						out << "-";

		if ( is_geometry_admissible[ii] ) 	out << "G";
		else									out << "-";

		
		out << endl;
	}
}
	
bool	Chain_binary::
extract_fragment (
	const int start_pos, 
	const int length, 
	double *coord )
{
	if (start_pos < 0 || (start_pos+length-1) > number_of_residues_ )
		return false;

//	if ( ! coord_set_ ) 				get_coord_set (); 

	if ( ! is_geometry_admissible_ ) 		
		get_is_geometry_admissible (); 


	bool g_flag = true;
	for (int ii=start_pos; ii<start_pos+length; ii++)
		g_flag &= is_geometry_admissible_[ii];

	if (!g_flag)
		return false;
	else
	{
		if (! coord_set_)
			get_coord_set();

		memcpy (coord,coord_set_+start_pos*9,9*length*sizeof(double) );
		return true;
	}
}

void	Chain_binary::
save_pdb_fragment (
	const int start_pos,
	const int length, 
	const string & path_to_pdb)
{

	ofstream out( path_to_pdb.c_str() );

	if ( ! out) 	
	{
		cout		<<	" can't create file " << path_to_pdb  << endl;
		log_stream	<<	" can't create file " << path_to_pdb  << endl;
		exit (1);
	}

	double *coord = new double [length *9];
	
	extract_fragment (start_pos,length,coord );

	vector <string> residue_names			= get_residue_names () ;
	vector <string> in_chain_residue_number = get_in_chain_residue_number (); 

	int counter = 1;
	for (int ii = start_pos; ii < start_pos + length; ii++)
	{
		for (int kk=0;kk<3;kk++) 
		{
			PutVa ("ATOM  ",out,6, 6 ,'r');				//   1 -  6        Record name     "ATOM  "         
			PutVa (counter,	out,3, 3 ,'r');	   			//   7 - 11        Atom serial number.

			counter++;

			out  << " ";						//** 12  space
			
			string atom_name;
			switch (kk) 
			{
				case 0: atom_name = " N ";break;
				case 1: atom_name = " CA";break;
				case 2: atom_name = " C ";break;
			}

			PutVa (atom_name ,out,5, 5 ,'l');				
			out << "  ";

			PutVa (residue_names[ii],				out,3, 3,'l');				//   18 - 20        Residue name    resName       Residue name.                         
			out << "  ";

			PutVa (in_chain_residue_number[ii],out,4, 4,'r');	//   23 - 26        Integer         resSeq        Residue sequence number.              

			PutVa ("    ",			out,3, 3,'l');	//** 28-30 пропуск 

			PutVaDouble (coord[9*(ii-start_pos) +3*kk      ],	out,8, 3,'r');	//31 - 38        Real(8.3)       x             Orthogonal coordinates for X in   Angstroms.                           
			PutVaDouble (coord[9*(ii-start_pos) +3*kk  + 1 ],	out,8, 3,'r');	//39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in       
			PutVaDouble (coord[9*(ii-start_pos) +3*kk  + 2 ],	out,8, 3,'r');	//47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in       

			out << endl;
		}
	}
	 delete [] coord;
}

void	Chain_binary::
save_pdb_fragment (
	const double *coord,
	const int length, 
	const string & path_to_pdb,
	vector <string> & residue_names,
	vector <string> & in_chain_residue_number)
{

	ofstream out( path_to_pdb.c_str() );

	if ( ! out) 	
	{
		cout		<<	" can't create file " << path_to_pdb  << endl;
		log_stream	<<	" can't create file " << path_to_pdb  << endl;
		exit (1);
	}


//	vector <string> residue_names			= get_residue_names () ;
//	vector <string> in_chain_residue_number = get_in_chain_residue_number (); 

	int counter = 1;

	
	for (int ii = 0; ii < length; ii++)
	{
		for (int kk=0;kk<3;kk++) 
		{
			PutVa ("ATOM  ",out,6, 6 ,'r');				//   1 -  6        Record name     "ATOM  "         
			PutVa (counter,	out,3, 3 ,'r');	   			//   7 - 11        Atom serial number.

			counter++;

			out  << " ";						//** 12  space
			
			string atom_name;
			switch (kk) 
			{
				case 0: atom_name = " N ";break;
				case 1: atom_name = " CA";break;
				case 2: atom_name = " C ";break;
			}

			PutVa (atom_name ,out,5, 5 ,'l');				
			out << "  ";

			PutVa (residue_names[ii],				out,3, 3,'l');				//   18 - 20        Residue name    resName       Residue name.                         
			out << "  ";

			PutVa (in_chain_residue_number[ii],out,4, 4,'r');	//   23 - 26        Integer         resSeq        Residue sequence number.              

			PutVa ("    ",			out,3, 3,'l');	//** 28-30 пропуск 

			PutVaDouble (coord[9*ii +3*kk      ],	out,8, 3,'r');	//31 - 38        Real(8.3)       x             Orthogonal coordinates for X in   Angstroms.                           
			PutVaDouble (coord[9*ii +3*kk  + 1 ],	out,8, 3,'r');	//39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in       
			PutVaDouble (coord[9*ii +3*kk  + 2 ],	out,8, 3,'r');	//47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in       

			out << endl;
		}
	}
}

// This version contain chain ID
void	Chain_binary::
save_pdb_fragment (
	const double *coord,
	const int length, 
	const string & path_to_pdb,
	vector <string> & residue_names,
	vector <string> & in_chain_residue_number,
	char chain_ID)
{

	ofstream out( path_to_pdb.c_str() );

	if ( ! out) 	
	{
		cout		<<	" can't create file " << path_to_pdb  << endl;
		log_stream	<<	" can't create file " << path_to_pdb  << endl;
		exit (1);
	}


//	vector <string> residue_names			= get_residue_names () ;
//	vector <string> in_chain_residue_number = get_in_chain_residue_number (); 

	int counter = 1;

	
	for (int ii = 0; ii < length; ii++)
	{
		for (int kk=0;kk<3;kk++) 
		{
			PutVa ("ATOM  ",out,6, 6 ,'r');				//   1 -  6        Record name     "ATOM  "         
			PutVa (counter,	out,3, 3 ,'r');	   			//   7 - 11        Atom serial number.

			counter++;

			out  << " ";						//** 12  space
			
			string atom_name;
			switch (kk) 
			{
				case 0: atom_name = " N ";break;
				case 1: atom_name = " CA";break;
				case 2: atom_name = " C ";break;
			}

			PutVa (atom_name ,out,5, 5 ,'l');				
			out << "  ";

			PutVa (residue_names[ii],				out,3, 3,'l');				//   18 - 20        Residue name    resName       Residue name.                         
			out << " ";
			out << chain_ID;

			PutVa (in_chain_residue_number[ii],out,4, 4,'r');	//   23 - 26        Integer         resSeq        Residue sequence number.              

			PutVa ("    ",			out,3, 3,'l');	//** 28-30 пропуск 

			PutVaDouble (coord[9*ii +3*kk      ],	out,8, 3,'r');	//31 - 38        Real(8.3)       x             Orthogonal coordinates for X in   Angstroms.                           
			PutVaDouble (coord[9*ii +3*kk  + 1 ],	out,8, 3,'r');	//39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in       
			PutVaDouble (coord[9*ii +3*kk  + 2 ],	out,8, 3,'r');	//47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in       

			out << endl;
		}
	}
}


void	Chain_binary::
fragment_to_principal_axes (
	const int start_pos, 
	const int length, 
	double *coord )
{
	extract_fragment ( start_pos, length, coord );

	double val[3],vect[3][3],a[3][3];   // val - eigen val; vect - eigen vect; a - inertia tensor
	
	centroid ( coord,  length);
	inertia_tensor ( a, coord, length);

	bool true_flag = EigenSystem3D( a, vect, val );

	double det= det_3x3 (vect);
	if (det < 0) 
	{
		vect[0][2] *= -1;
		vect[1][2] *= -1;
		vect[2][2] *= -1;
	}
	det= det_3x3 (vect);

	assert ( det > 0.99 );

	rotate_fragment (length,coord,vect);
}


vector < vector < double > >  Chain_binary::
positioning_chain_by_clasters_set ( 
	double **claster_motif_coordinates,
	const int fragment_length,
	const int number_of_classes,
	const int sequence_length) 
{

//		string sequence = get_sequence ();
		int current_number_of_fragments = sequence_length - fragment_length + 1;

		vector < vector < double > >  set_of_coordinate_in_clasters_system;
		set_of_coordinate_in_clasters_system.resize( current_number_of_fragments );
		for (int ii=0;ii<current_number_of_fragments;ii++) 
		{
			set_of_coordinate_in_clasters_system[ii].resize( number_of_classes );
			for (int ttt=0;ttt<number_of_classes ;ttt++)
				set_of_coordinate_in_clasters_system[ii][ttt] = -1;
		}

		double *cu_cord_set  = new double [fragment_length*9];
		double *w	= new double [fragment_length *9];
		for (int kk=0;kk<(fragment_length *9);kk++)
			w[kk] = 1;


	//	is_geometry_admissible_	=	get_is_geometry_admissible();
//		is_there_coord_			=	get_is_there_coord();

		int kk;

		for (kk=0; kk< current_number_of_fragments ;kk++ )
		{
			
			extract_fragment(kk,fragment_length,cu_cord_set );

// check that coorditates  presents and admissible
	/*		bool is_geom_flag = true;
			bool is_coord_flag = true;
			int right_stop = kk + fragment_length;
			for (int pp=kk; pp< right_stop; pp++)
			{
				is_geom_flag  = is_geom_flag & is_geometry_admissible_	[pp];
				is_coord_flag = is_coord_flag & is_there_coord_			[pp];
			}
			if (!is_geom_flag || !is_coord_flag) 
				continue ;			// енсли плохие координаты или  их нету то пропускаем
*/

			for (int zz=0;zz<number_of_classes;zz++)
			{

				bool error_flag= true;
				double distance  = kabsch_rmsd	(
										cu_cord_set , 
										claster_motif_coordinates[zz], 
										w, 
										(3*fragment_length),
										error_flag);

				set_of_coordinate_in_clasters_system[kk][zz] = distance ;
			}
		}

		delete [] cu_cord_set;		delete [] w;

		return set_of_coordinate_in_clasters_system;
}

vector < vector < double > >   & Chain_binary::
positioning_chain_by_clasters_set ( 

	double **claster_motif_coordinates,
	const int fragment_length,
	const int number_of_classes	) 
{

		string sequence = get_sequence ();
		int current_number_of_fragments = sequence.length() - fragment_length + 1;

		vector < vector < double > >  set_of_coordinate_in_clasters_system;
		set_of_coordinate_in_clasters_system.resize( current_number_of_fragments );
		for (int ii=0;ii<current_number_of_fragments;ii++) 
		{
			set_of_coordinate_in_clasters_system[ii].resize( number_of_classes );
			for (int ttt=0;ttt<number_of_classes ;ttt++)
				set_of_coordinate_in_clasters_system[ii][ttt] = -1;
		}

		double *cu_cord_set  = new double [fragment_length*9];
		double *w	= new double [fragment_length *9];
		for (int kk=0;kk<(fragment_length *9);kk++)
			w[kk] = 1;


		is_geometry_admissible_	=	get_is_geometry_admissible();
		is_there_coord_			=	get_is_there_coord();

		int kk;

		for (kk=0; kk< current_number_of_fragments ;kk++ )
		{
			
			extract_fragment(kk,fragment_length,cu_cord_set );

// check that coorditates  presents and admissible
			bool is_geom_flag = true;
			bool is_coord_flag = true;
			//int right_stop = min (kk + window_size,current_number_of_fragments ) ;
			int right_stop = kk + fragment_length;
			for (int pp=kk; pp< right_stop; pp++)
			{
				is_geom_flag  = is_geom_flag & is_geometry_admissible_	[pp];
				is_coord_flag = is_coord_flag & is_there_coord_			[pp];
			}
			if (!is_geom_flag || !is_coord_flag) 
				continue ;			// енсли плохие координаты или  их нету то пропускаем


			for (int zz=0;zz<number_of_classes;zz++)
			{

				bool error_flag= true;
				double distance  = kabsch_rmsd	(
										cu_cord_set , 
										claster_motif_coordinates[zz], 
										w, 
										(3*fragment_length),
										error_flag);

				set_of_coordinate_in_clasters_system[kk][zz] = distance ;
			}
		}

		delete [] cu_cord_set;		delete [] w;

		return set_of_coordinate_in_clasters_system;
}

void Chain_binary::
positioning_chain_by_clasters_set(
	double **claster_motif_coordinates,
	const int fragment_length,
	const int number_of_classes,
	vector < vector < double > >  & set_of_coordinate_in_clasters_system)
{

	string sequence = get_sequence();
	int current_number_of_fragments = sequence.length() - fragment_length + 1;

//	vector < vector < double > >  set_of_coordinate_in_clasters_system;
	set_of_coordinate_in_clasters_system.resize(current_number_of_fragments);
	for (int ii = 0; ii<current_number_of_fragments; ii++)
	{
		set_of_coordinate_in_clasters_system[ii].resize(number_of_classes);
		for (int ttt = 0; ttt<number_of_classes; ttt++)
			set_of_coordinate_in_clasters_system[ii][ttt] = -1;
	}

	double *cu_cord_set = new double[fragment_length * 9];
	double *w = new double[fragment_length * 9];
	for (int kk = 0; kk<(fragment_length * 9); kk++)
		w[kk] = 1;


	is_geometry_admissible_ = get_is_geometry_admissible();
	is_there_coord_ = get_is_there_coord();

	int kk;

	for (kk = 0; kk< current_number_of_fragments; kk++)
	{

		extract_fragment(kk, fragment_length, cu_cord_set);

		// check that coorditates  presents and admissible
		bool is_geom_flag = true;
		bool is_coord_flag = true;
		//int right_stop = min (kk + window_size,current_number_of_fragments ) ;
		int right_stop = kk + fragment_length;
		for (int pp = kk; pp< right_stop; pp++)
		{
			is_geom_flag = is_geom_flag & is_geometry_admissible_[pp];
			is_coord_flag = is_coord_flag & is_there_coord_[pp];
		}
		if (!is_geom_flag || !is_coord_flag)
			continue;			// енсли плохие координаты или  их нету то пропускаем


		for (int zz = 0; zz<number_of_classes; zz++)
		{

			bool error_flag = true;
			double distance = kabsch_rmsd(
				cu_cord_set,
				claster_motif_coordinates[zz],
				w,
				(3 * fragment_length),
				error_flag);

			set_of_coordinate_in_clasters_system[kk][zz] = distance;
		}
	}

	delete[] cu_cord_set;		delete[] w;

}


vector < vector < double > > 
two_chain_distance_set (
	const string & pdb_chain_ID_1,
	const string & pdb_chain_ID_2,
	const int fragment_length )
{

	Chain_binary chain_1( pdb_chain_ID_1 );
	Chain_binary chain_2( pdb_chain_ID_2 );

	string sequence_1			= chain_1.get_sequence ();
	int number_of_fragments_1	= sequence_1.length() - fragment_length + 1;

	string sequence_2			= chain_2.get_sequence ();
	int number_of_fragments_2	= sequence_2.length() - fragment_length + 1;


	vector < vector < double > >  distance_set;
	distance_set.resize( number_of_fragments_1 );
	for (int ii=0;ii<number_of_fragments_1;ii++) 
	{
		distance_set[ii].resize( number_of_fragments_2 );
		for (int ttt=0;ttt<number_of_fragments_2 ;ttt++)
			distance_set[ii][ttt] = -1;
	}

	double *w	= new double [fragment_length *9];
	for (int kk=0;kk<(fragment_length *9);kk++)
		w[kk] = 1;


	double *cord_set_1  = new double [fragment_length*9];
	double *cord_set_2  = new double [fragment_length*9];

	


	bool* is_geometry_admissible_1 = (bool *) new bool [sequence_1.length()];
	bool* is_geometry_admissible_2 = (bool *) new bool [sequence_2.length()];

	memcpy (is_geometry_admissible_1, chain_1.get_is_geometry_admissible(),  sequence_1.length() *sizeof (bool) );
	memcpy (is_geometry_admissible_2, chain_2.get_is_geometry_admissible(),  sequence_2.length()*sizeof (bool) );


//	is_geometry_admissible_1	=	chain_1.get_is_geometry_admissible();
//	is_geometry_admissible_2	=	chain_2.get_is_geometry_admissible();

	bool* is_there_coord_1 = (bool *) new bool [sequence_1.length()];
	bool* is_there_coord_2 = (bool *) new bool [sequence_2.length()];

	memcpy (is_there_coord_1 , chain_1.get_is_there_coord(),  sequence_1.length()*sizeof (bool) );
	memcpy (is_there_coord_2 , chain_2.get_is_there_coord(),  sequence_2.length()*sizeof (bool) );

//	is_there_coord_1			=	chain_1.get_is_there_coord();
//	is_there_coord_2			=	chain_2.get_is_there_coord();

	int kk;
	for (kk=0; kk< number_of_fragments_1;kk++ )
	{
		
		chain_1.extract_fragment(kk,fragment_length,cord_set_1 );

		bool is_geom_flag = true;
		bool is_coord_flag = true;
		int right_stop = kk + fragment_length;
		for (int pp=kk; pp< right_stop; pp++)
		{
			is_geom_flag  = is_geom_flag & is_geometry_admissible_1	[pp];
			is_coord_flag = is_coord_flag & is_there_coord_1		[pp];
		}
		if (!is_geom_flag || !is_coord_flag) 			
			continue ;			// если плохие координаты или  их нету то пропускаем


			for (int zz=0;zz<number_of_fragments_2;zz++)
			{

			bool is_geom_flag = true;
			bool is_coord_flag = true;
			int right_stop = zz + fragment_length;
			for (int pp=zz; pp< right_stop; pp++)
			{
				is_geom_flag  = is_geom_flag & is_geometry_admissible_2	[pp];
				is_coord_flag = is_coord_flag & is_there_coord_2		[pp];
			}
			if (!is_geom_flag || !is_coord_flag) 			
				continue ;			// если плохие координаты или  их нету то пропускаем


				chain_2.extract_fragment(zz,fragment_length,cord_set_2 );

				bool error_flag= true;
				double distance  = kabsch_rmsd	(
										cord_set_1, 
										cord_set_2, 
										w, 
										(3*fragment_length),
										error_flag);

				distance_set[kk][zz] = distance ;

//				cout << kk << " " << zz << endl;
			}

	}

	delete [] is_geometry_admissible_1;	delete [] is_geometry_admissible_2;
	delete [] is_there_coord_1;			delete [] is_there_coord_2;
	delete [] cord_set_1;				delete [] cord_set_2;



	return distance_set;
}

/*
vector < vector < int > >  Chain_binary::
analyse_setted_regular_structure_presence_in_chain(
	const int fragment_length )
{


}
*/

bool align_two_chains (
	double *cord_set_1,
	double *cord_set_2,
	int fragment_length)
{

	double *w	= new double [fragment_length *9];
	for (int kk=0;kk<(fragment_length *9);kk++)
		w[kk] = 1;

	bool error_flag= true;
	double distance  = kabsch_rmsd	(
							cord_set_1, 
							cord_set_2, 
							w, 
							(3*fragment_length),
							error_flag);

	return error_flag;

}

vector < double > Chain_binary::
get_torsion_angles ()
{

	vector < double > torsion_angles;
		fill_up_fragment_torsion_angles ( 
			coord_set_,		
			number_of_residues_,
			torsion_angles, 
			'r');

	return torsion_angles;
}

void Chain_binary::center_of_mass(vector < double > & c_m_coord)
{
	c_m_coord.resize(3);

	double x, y, z;
	double x_ave=0, y_ave=0, z_ave=0;

	int true_coord_counter = 0; 
	for (int ii = 0; ii < number_of_residues_; ii++)
	{
		if (is_geometry_admissible_[ii])
		{
			x_ave += coord_set_[9 * ii];			// That's coordinate for N atom!
			y_ave += coord_set_[9 * ii + 1];
			z_ave += coord_set_[9 * ii + 2];
			true_coord_counter++;
		}
//		else
//			cout << "aga!";
	}

	c_m_coord[0] = x_ave / true_coord_counter;
	c_m_coord[1] = y_ave / true_coord_counter;
	c_m_coord[2] = z_ave / true_coord_counter;
	
}