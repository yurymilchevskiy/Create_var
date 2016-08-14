#pragma warning( disable : 4786 )

#include "Chain_binary_test.h"
#include "Chain_binary.h"

#include "../CommonFunc.h"

#include <fstream>
#include <cassert>

#include "../Geometry_util/Fragmrent_tranformation.h"

//	void kabsch_stolen(double *ref, double *vec, double *w, const int num);
//	double rmsd(double *ref, double *vec, double *w, const int num);

#include "../by_Qmol/kabsch_stolen.h"

using namespace std;

extern ofstream log_stream;

Chain_binary_test::
~Chain_binary_test()
{
	cout << "Chain_binary_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}


void Chain_binary_test::
test1 ()
{
	{ Chain_binary cb( "2VB1A");cb.print_protocol();}
/*	{ Chain_binary cb( "1R6JA");cb.print_protocol();}
	{ Chain_binary cb( "2DSXA");cb.print_protocol();}
	{ Chain_binary cb( "1UCSA");cb.print_protocol();}
	{ Chain_binary cb( "1P9GA");cb.print_protocol();}
	{ Chain_binary cb( "2I17A");cb.print_protocol();}
	{ Chain_binary cb( "1GCIA");cb.print_protocol();}
*/
}

void Chain_binary_test::
remind_fragmnet_test ()
{
	ofstream out( "TEST/SobKamni/2VB1A.check" );

	if ( ! out) 	
	{
		cout		<<	" can't open test file " << endl;
		log_stream	<<	" can't open test file " << endl;
		exit (1);
	}

	 Chain_binary cb( "2VB1A");cb.print_protocol();
	 cb.get_coord_set();

 	 int number_of_residues = cb.get_number_of_residues();
	 int length = 10;
	 double *coord = new double [length *9*5];

	 
	vector <string> total_residue_names				= cb.get_residue_names () ;
	vector <string> total_in_chain_residue_number	= cb.get_in_chain_residue_number (); 


	 for ( int ii=0; ii<number_of_residues-length+1; ii++ ) 
	 {
		bool flag = cb.extract_fragment(ii,length,coord);
		string is_OK;
		if (flag)	is_OK = "OK!";
		else		is_OK = "---";

		PutVa (is_OK,out,5,1,'l');

		for (int kk=0;kk<length;kk++ )
		{
			PutVa (total_residue_names[ii+kk],out,3,0,'l');
		}

		out << endl;
	 }

	out.close();
}

void Chain_binary_test::
get_fragment_test1 ()
{

	ofstream out( "TEST/Chain_binary_test_get_fragment_test1" );

	if ( ! out) 	
	{
		cout		<<	" can't open test file " << endl;
		log_stream	<<	" can't open test file " << endl;
		exit (1);
	}

	 Chain_binary cb( "2VB1A");cb.print_protocol();
	 cb.get_coord_set();
	
	 int number_of_residues = cb.get_number_of_residues();
	 int length = 5;
	 double *coord = new double [length *9*5];

	 for ( int ii=0; ii<number_of_residues; ii++ ) 
	 {
		bool flag = cb.extract_fragment(ii,length,coord);
		string is_OK;
		if (flag)	is_OK = "OK!";
		else		is_OK = "---";

		PutVa (is_OK,out,5,1,'l');
		out << endl;
	 }

	out.close();

		int start_pos	= 5;
		length	= 5;
 
		string path_to_pdb = "TEST/Chain_binary_test_pdb.pdb";
		cb.save_pdb_fragment (start_pos,length,path_to_pdb);


		vector <string> total_residue_names				= cb.get_residue_names () ;
		vector <string> total_in_chain_residue_number	= cb.get_in_chain_residue_number (); 

		vector <string> residue_names;				
		vector <string> in_chain_residue_number;
		int ii;
		for ( ii=0;ii<length;ii++)
		{
			residue_names.push_back(total_residue_names[start_pos+ii]); 
			in_chain_residue_number.push_back(total_in_chain_residue_number[start_pos+ii]);
		}

		cb.extract_fragment (start_pos,length,coord );
		cb.save_pdb_fragment (	coord, length, path_to_pdb, residue_names, in_chain_residue_number);


		centroid ( coord,  length);
		path_to_pdb = "TEST/Chain_binary_test_centroid_pdb.pdb";
//		cb.extract_fragment (start_pos,length,coord );
		cb.save_pdb_fragment (	coord, length, path_to_pdb, residue_names, in_chain_residue_number);


		ofstream out1( "TEST/Chain_binary_test_inertia_tensor" );

		if ( ! out1) 	
		{
			cout		<<	" can't open test file " << endl;
			log_stream	<<	" can't open test file " << endl;
			exit (1);
		}

		double ti[3][3];

		inertia_tensor ( ti, coord, length);
    //    print_3x3_matris (ti,length, out1 );


	delete [] coord;
}

#include "../by_Qmol/EigenValues3D.h"

void Chain_binary_test::
EigenValues3D_test ()
{
//	int length = 5;

	 Chain_binary cb( "2VB1A");cb.print_protocol();
	 cb.get_coord_set();

	
	 int number_of_residues = cb.get_number_of_residues();
	 int length = 5;
//	 double *coord = new double [length *9*5];
	
	 double a[3][3];
	 double *coord = new double [length *9*5];

	int start_pos	= 5;
	length	= 5;
 
	cb.extract_fragment (start_pos,length,coord );
	centroid ( coord,  length);
	inertia_tensor ( a, coord, length);

	double val[3];
	double vect[3][3];

	bool true_flag = EigenSystem3D( a, vect, val );

	rotate_fragment (length,coord,vect);

	double b[3][3];

	inertia_tensor ( b, coord, length);

	double Error=0;
	Error += fabs( b[0][0] - val[0]);
	Error += fabs( b[1][1] - val[1]);
	Error += fabs( b[2][2] - val[2]);

	Error += fabs( b[0][1]);
	Error += fabs( b[0][2]);
	Error += fabs( b[1][0]);
	Error += fabs( b[2][1]);
	Error += fabs( b[1][2]);
	Error += fabs( b[2][0]);

	test_ ("CHecking eigen rotation", Error < 0.000001);

	double determinant = det_3x3 (vect);

		vector <string> total_residue_names				= cb.get_residue_names () ;
		vector <string> total_in_chain_residue_number	= cb.get_in_chain_residue_number (); 

		vector <string> residue_names;				
		vector <string> in_chain_residue_number;

		for (int ii=0;ii<length;ii++)
		{
			residue_names.push_back(total_residue_names[start_pos+ii]); 
			in_chain_residue_number.push_back(total_in_chain_residue_number[start_pos+ii]);
		}

		string path_to_pdb = "TEST/Chain_binary_test_ready.pdb";

		cb.save_pdb_fragment (	coord, length, path_to_pdb, residue_names, in_chain_residue_number);


	delete [] coord;
}
#define LENGTH 5
//11 я тверск ямская 5 
void Chain_binary_test::fragment_to_principal_axes_test ()
{
	 Chain_binary cb( "2VB1A");
	 cb.get_coord_set();
	 cb.get_is_geometry_admissible();
	 int number = cb.get_number_of_residues();
 
			 double *coord= new double [(LENGTH *9)];
			 cb.extract_fragment (90,5,coord );

 	 double *coord1= new double [(LENGTH *9)];
 	 double *coord2= new double [(LENGTH *9)];
 	 double *coord2a= new double [(LENGTH *9)];
 	 double *coord2b= new double [(LENGTH *9)];
 	 double *coord2c= new double [(LENGTH *9)];

	int start_pos1 = 4;
	int start_pos2 = 35;

	  cb.fragment_to_principal_axes (
			start_pos1, 
			LENGTH, 
			coord1 );

	  cb.fragment_to_principal_axes (
			start_pos2, 
			LENGTH, 
			coord2 );


	  double dX1 = coord1[42]  - coord1[0];
	  double dY1 = coord1[43]  - coord1[1];
	  double dZ1 = coord1[44]  - coord1[2];


  	  double dX2 = coord2[42]  - coord2[0];
	  double dY2 = coord2[43]  - coord2[1];
	  double dZ2 = coord2[44]  - coord2[2];



	memcpy (coord2a,coord2,(LENGTH *9)*sizeof(double));
	memcpy (coord2b,coord2,(LENGTH *9)*sizeof(double));
	memcpy (coord2c,coord2,(LENGTH *9)*sizeof(double));

	double matr2a [3][3] = {
		{-1,0,0},
		{0,-1,0},
		{0,0,1}};

	double matr2b [3][3] = {
		{1,0,0},
		{0,-1,0},
		{0,0,-1}};

	double matr2c [3][3] = {
		{-1,0,0},
		{0, 1,0},
		{0, 0,-1}};

/*	double matr2a [3][3] = {
		{0,0,-1},
		{0,1,0},
		{1,0,0}};

	double matr2b [3][3] = {
		{0,0,1},
		{0,-1,0},
		{1,0,0}};

	double matr2c [3][3] = {
		{0,0,1},
		{0,1,0},
		{-1,0,0}};*****/
	rotate_fragment (LENGTH,coord2a,matr2a);
	rotate_fragment (LENGTH,coord2b,matr2b);
	rotate_fragment (LENGTH,coord2c,matr2c);


  	  double dX2a = coord2a[42]  - coord2a[0];
	  double dY2a = coord2a[43]  - coord2a[1];
	  double dZ2a = coord2a[44]  - coord2a[2];


  	  double dX2b = coord2b[42]  - coord2b[0];
	  double dY2b = coord2b[43]  - coord2b[1];
	  double dZ2b = coord2b[44]  - coord2b[2];

		vector <string> total_residue_names				= cb.get_residue_names () ;
		vector <string> total_in_chain_residue_number	= cb.get_in_chain_residue_number (); 

	  	vector <string> residue_names_1;				
		vector <string> in_chain_residue_number_1;
	  	vector <string> residue_names_2;				
		vector <string> in_chain_residue_number_2;

		for ( int ii=0;ii<LENGTH;ii++)
		{
			residue_names_1.push_back(total_residue_names[start_pos1+ii]); 
			in_chain_residue_number_1.push_back(total_in_chain_residue_number[start_pos1+ii]);
			residue_names_2.push_back(total_residue_names[start_pos2+ii]); 
			in_chain_residue_number_2.push_back(total_in_chain_residue_number[start_pos2+ii]);
		}
		string path_to_pdb1 = "TEST/a1.pdb";
		string path_to_pdb2 = "TEST/a2.pdb";
		string path_to_pdb2a = "TEST/a2a.pdb";
		string path_to_pdb2b = "TEST/a2b.pdb";
		string path_to_pdb2c = "TEST/a2c.pdb";

		cb.save_pdb_fragment (	coord1, LENGTH, path_to_pdb1, residue_names_1, in_chain_residue_number_1);

		cb.save_pdb_fragment (	coord2, LENGTH, path_to_pdb2, residue_names_2, in_chain_residue_number_2);
		cb.save_pdb_fragment (	coord2a, LENGTH, path_to_pdb2a, residue_names_2, in_chain_residue_number_2);
		cb.save_pdb_fragment (	coord2b, LENGTH, path_to_pdb2b, residue_names_2, in_chain_residue_number_2);
		cb.save_pdb_fragment (	coord2c, LENGTH, path_to_pdb2c, residue_names_2, in_chain_residue_number_2);


	delete []	coord;

	delete []	coord1;
	delete []	coord2;
}

void Chain_binary_test::
kabsch_stolen_test()

{

	 Chain_binary cb( "2VB1A");
	 cb.get_coord_set();
	 cb.get_is_geometry_admissible();
	 int num = cb.get_number_of_residues();
 
		//	 double *coord= new double [(LENGTH *9)];
		//	 cb.extract_fragment (90,5,coord );

 	 double *ref= new double [(LENGTH *9)];
 	 double *vec= new double [(LENGTH *9)];
	 double *w	= new double [(LENGTH *9)];

	 for (int kk=0;kk<(LENGTH *9);kk++)
		w[kk] = 1;
 
	int start_pos1 = 4;
	int start_pos2 = 101;

/*	  cb.fragment_to_principal_axes (
			start_pos1, 
			LENGTH, 
			ref );

	  cb.fragment_to_principal_axes (
			start_pos2, 
			LENGTH, 
			vec );
*/

	bool is_exist1 = cb.extract_fragment ( start_pos1, LENGTH, ref );
	bool is_exist2 = cb.extract_fragment ( start_pos2, LENGTH, vec );
	assert (is_exist1 & is_exist2 );
	
//	centroid ( double *coord, const int length)
	centroid ( ref, LENGTH);
	centroid ( vec, LENGTH);


	bool error_flag=true;;

	kabsch_stolen	(ref, vec, w, (3*LENGTH),error_flag );
	double _RMSD = rmsd			(ref, vec, w, (3*LENGTH) );
	cout << _RMSD << endl;



  	string path_to_pdb1 = "TEST/a1.pdb";
	string path_to_pdb2 = "TEST/a2.pdb";

	vector <string> total_residue_names				= cb.get_residue_names () ;
	vector <string> total_in_chain_residue_number	= cb.get_in_chain_residue_number (); 


	vector <string> residue_names_1;				
	vector <string> in_chain_residue_number_1;
	vector <string> residue_names_2;				
	vector <string> in_chain_residue_number_2;

	for ( int ii=0;ii<LENGTH;ii++)
	{
		residue_names_1.push_back(total_residue_names[start_pos1+ii]); 
		in_chain_residue_number_1.push_back(total_in_chain_residue_number[start_pos1+ii]);
		residue_names_2.push_back(total_residue_names[start_pos2+ii]); 
		in_chain_residue_number_2.push_back(total_in_chain_residue_number[start_pos2+ii]);
	}


	cb.save_pdb_fragment (	ref, LENGTH, path_to_pdb1, residue_names_1, in_chain_residue_number_1);
	cb.save_pdb_fragment (	vec, LENGTH, path_to_pdb2, residue_names_2, in_chain_residue_number_2);


	delete [] ref;
	delete [] vec;
	delete [] w;

}

void Chain_binary_test::
kabsch_rmsd_test()
{

	 Chain_binary cb( "2VB1A");
	 cb.get_coord_set();
	 cb.get_is_geometry_admissible();
	 int num = cb.get_number_of_residues();

 	 int start_pos1 = 4;
	 int start_pos2 = 101;


 	 double *ref= new double [(LENGTH *9)];
 	 double *vec= new double [(LENGTH *9)];
	 double *w	= new double [(LENGTH *9)];


 	 bool is_exist1 = cb.extract_fragment ( start_pos1, LENGTH, ref );
	 bool is_exist2 = cb.extract_fragment ( start_pos2, LENGTH, vec );
	 assert (is_exist1 & is_exist2 );

 	 for (int kk=0;kk<(LENGTH *9);kk++)		w[kk] = 1;
	 bool error_flag=true;;
	 double RMSD = kabsch_rmsd (ref, vec, w, (3*LENGTH), error_flag);

	 cout << "RMSD by kabsch_rmsd " << RMSD << endl;

	delete [] ref;
	delete [] vec;
	delete [] w;
}

void Chain_binary_test::
two_chain_distance_set_test()
{
	ofstream out( "TEST/two_chain_distance_set_test_length10" );

	if ( ! out) 	
	{
		cout		<<	" can't open test file " << endl;
		log_stream	<<	" can't open test file " << endl;
		exit (1);
	}


	vector < vector < double > >  distance_set =
		two_chain_distance_set (
			"119LA","2OB5A",	10);

	for (int ii=0;ii<distance_set.size();ii++)
	{
		for (int jj=0;jj<distance_set[ii].size();jj++)
			PutVaDouble (distance_set[ii][jj],out,7,3,'l'); 

		out << endl;
	}
}


void Chain_binary_test::
center_of_mass_test()
{
	Chain_binary cb("2VB1A"); 

	double *coord_set_ = cb.get_coord_set();
	cb.get_is_geometry_admissible();
	cb.get_is_there_coord();

	vector <double> c_m_coord; 
	cb.center_of_mass(c_m_coord);

	cout << c_m_coord[0] << endl;
	cout << c_m_coord[1] << endl;
	cout << c_m_coord[2] << endl;

}