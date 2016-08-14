#pragma warning( disable : 4786 )

#include "Frequency_extrapolation_test.h"
#include "Frequency_extrapolation.h"
 //#include "Frequency_extrapolation_ZIP.h"

#include "../Censorship.h"

#include "../CommonFunc.h"

#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

extern Censorship configuration;
extern ofstream log_stream;

Frequency_extrapolation_test::~Frequency_extrapolation_test()
{
	cout << "Frequency_extrapolation_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}

void Frequency_extrapolation_test::
check_integrity ()
{
	
}


void Frequency_extrapolation_test::
fill_up_setted_basis_structure()
{

	{
	//	Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W3",	COMMON_USAGE);
//		fre.show_freq_data_slow(); 
	}

	{
	//	Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W4N",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
	//	fre.show_freq_data_slow(); 
	}

	{
	//	Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W7_trivial_PG",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
	//	fre.show_freq_data_slow(); 
	}




	{
//		Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_w9_tail",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
	//	fre.show_freq_data_slow(); 
	}


	{
//	Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W5_noDEG_ZIP",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
//	Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W5_noDEG_ZIP",		COMMON_USAGE );
//	fre.show_freq_data_slow(); 
	}


	{
//		Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W6_3_trivial_PG",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
//		fre.show_freq_data_slow(); 
	}

	{
		Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W6_4_trivial_PG",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
	//	fre.show_freq_data_slow(); 
	}

	{
		Frequency_extrapolation	fre   ( "_SETTED_BASIS_STRUCTURES_3/c30_W4C",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
		//fre.show_freq_data_slow(); 
	}





}
/*
void Frequency_extrapolation_test::
check_integrity ()
{
	Frequency_extrapolation	fre   ( "c30_W5_noDEG_ZIP",		COMMON_USAGE );

	int number_of_elements = fre.get_number_of_elements () ;
	for (int ii =0; ii< number_of_elements;ii++)
	{
		int occurence = fre.get_occurence (ii);
		vector <double> single_distance_to_clusters_sum =  	fre.get_single_distance_to_clusters_sum	 ( ii);

		if (occurence!=0)
		{
			log_stream << ii << " " << occurence  << " " << single_distance_to_clusters_sum[0] << endl;
			cout << ii << " " << occurence  << " " << single_distance_to_clusters_sum[0] << endl;

		}
	}
	
}  
*/
// для 
void Frequency_extrapolation_test::analyse_cluster_only_test()
{
	Frequency_extrapolation	fre   ( "c30_W5_noDEG_ZIP",		ANALYSE_CLUSTER_ONLY );
	vector <vector <double> > distance_mutual = fre.calc_cluster_mutual_distance ();


	ofstream  out( "D:/Agony/TEST/analyse_cluster_only_tes");

	if ( ! out)	{	
		log_stream << "CowardVariables_test: ERROR -  can't create file" << endl;
		cout       << "CowardVariables_test: ERROR -  can't create file" << endl;
		exit (1);	
	}

	int number_of_classes  = distance_mutual.size();
	
	for ( int ii=0;ii<number_of_classes;ii++)
	{
		for (int jj=0;jj<number_of_classes;jj++)
			PutVaDouble (distance_mutual[ii][jj],out,10,5,'l');

		out << endl;
	
	}





}

void Frequency_extrapolation_test::   // creates databases
fill_up()
{
	{
	 Frequency_extrapolation fre   ( "c30_w11_tail_GP",		COMMON_USAGE);
		fre.show_freq_data_slow(); 
	} 

	{
  	cout << "c30_W7_tail_GP_copy" << endl;
	 Frequency_extrapolation fre   ( "c30_W7_tail_GP_copy",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	}
	
	{
//	 Frequency_extrapolation fre   ( "c30_W7_tail_GP",		COMMON_USAGE);
//	 fre.show_freq_data_slow(); 
	}
/*
	{
	 Frequency_extrapolation fre   ( "c30_W4N",		COMMON_USAGE);
	 fre.show_freq_data_slow(); 
	}

	{
	 Frequency_extrapolation fre   ( "c30_W4C",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	}

	{
	 Frequency_extrapolation fre   ( "c30_W4C",		COMMON_USAGE);
	 fre.show_freq_data_slow(); 
	}
*/
}

void Frequency_extrapolation_test::   // creates databases
first_test_extended()
{

	{
	//	Frequency_extrapolation	fre   ( "DSSP_30/c30_W5_noDEG_ZIP",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);

	//	Frequency_extrapolation	fre   ( "DSSP_30/c30_W5_noDEG_ZIP",	COMMON_USAGE);
	//	fre.show_freq_data_slow(); 
	}


	{
	//	Frequency_extrapolation	fre   ( "DSSP_30/c30_W5_noDEG_ZIP_denom",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
		//fre.show_freq_data_slow(); 
	}



	{
//	 Frequency_extrapolation fre   ( "DSSP_30/c30_w11_tail_GP",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	//	fre.show_freq_data_slow(); 
	} 

	{
//	 Frequency_extrapolation fre   ( "DSSP_30/c30_w11_tail",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	//	fre.show_freq_data_slow(); 
	} 

	{
 // 	cout << "c30_W7_tail_GP_copy" << endl;	 
		Frequency_extrapolation fre   ( "DSSP_30/c30_W7_tail_GP",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	}


	{
//	 Frequency_extrapolation fre   ( "DSSP_30/c30_W4C",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	}

	{
	 Frequency_extrapolation fre   ( "DSSP_30/c30_W4N",		FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	// fre.show_freq_data_slow(); 
	}


	{
	//	Frequency_extrapolation	fre   ( "DSSP_30/c30_W6_3_trivial_PG",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
		//fre.show_freq_data_slow(); 
	}

	{
		//Frequency_extrapolation	fre   ( "DSSP_30/c30_W6_4_trivial_PG",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
		//fre.show_freq_data_slow(); 
	}

	{
	//	Frequency_extrapolation	fre   ( "DSSP_30/c30_W7_trivial_PG",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
		//fre.show_freq_data_slow(); 
	}

	{
		//Frequency_extrapolation	fre   ( "DSSP_30/c30_w9_tail",	FILL_UP_FREQUENCY_EXTRAPOLATION_STORE);
	//fre.show_freq_data_slow(); 
	}

}


void Frequency_extrapolation_test::
second_test()
{
	Frequency_extrapolation	fre ( "c30_w5_complex", COMMON_USAGE);

	int  occurence =  fre.get_occurence (3);

	vector <double>  const & single_distance_to_clusters_sum =  	fre.get_single_distance_to_clusters_sum					 ( 3);
	vector <double>  const & squared_distance_to_clusters_sum=  	fre.get_single_squared_distance_to_clusters_sum			 ( 3);
	vector <double>  const & single_inverse_distance_to_clusters_sum = 	fre.get_single_inverse_distance_to_clusters_sum			 ( 3);
	vector <double>  const & single_inverse_squared_distance_to_clusters_sum =  	fre.get_single_inverse_squared_distance_to_clusters_sum  ( 3);

	vector <double> const & tot_distance_to_clusters_sum   	= fre.get_tot_distance_to_clusters_sum				 ( );
	vector <double> const & tot_squared_distance_to_clusters_sum  =	fre.get_tot_squared_distance_to_clusters_sum		 ( );
	vector <double> const & tot_inverse_distance_to_clusters_sum  =	fre.get_tot_inverse_distance_to_clusters_sum		 ( );
	vector <double>  const & tot_inverse_squared_distance_to_clusters_sum = 	fre.get_tot_inverse_squared_distance_to_clusters_sum ( );
	int		total_sample_size = 			fre.get_total_sample_size							 ( ) ;

	int number_of_elements = fre.get_number_of_elements ();


}
/*
void Frequency_extrapolation_test::
show_test()
{
	Frequency_extrapolation	fre ( "W5_01234_dx2", COMMON_USAGE);

	fre.show_freq_data ("1A2PA");
	fre.show_freq_data ("1A6M_");

	fre.show_freq_data ("together");


}
void Frequency_extrapolation_test::
get_single_distance_record_test()
{

	Frequency_extrapolation	fre ( "W5_01234_dx2", COMMON_USAGE);

		string freq_data_file_name = 
			configuration.option_meaning("Path_to_Frequency_extrapolation_store") + 	
			"W5_01234_dx2" + 
			"/base/"			+ 
			"1A2PA" +  
			".freq_data";


		ifstream data_stream(freq_dat
			a_file_name .c_str() ,ios::binary);
		if ( ! data_stream)	
		{	
			log_stream << "write_single_chain_data(): ERROR -  can't create " << freq_data_file_name << endl;
			log_stream << "write_single_chain_data(): ERROR -  can't create " << freq_data_file_name << endl;
			exit (1);	
		}

		vector <double>   single_distance_record ;

		single_distance_record = 	fre.get_single_distance_record (2098 ,data_stream);
		cout << "2098 " << single_distance_record [0] << endl;

		single_distance_record = 	fre.get_single_distance_record (2114 ,data_stream);
		cout << "2114  " << single_distance_record [0] << endl;

		single_distance_record = 	fre.get_single_distance_record (2120 ,data_stream);
		cout << "2120   " << single_distance_record [0] << endl;


		cout << "******" << endl;
		vector < int > occurence = fre.get_occurence (data_stream);

		cout << "ocuurence 2098:  " <<  occurence [2098] << endl;
		cout << "ocuurence 2114:  " <<  occurence [2114] << endl;
		cout << "ocuurence 2120:  " <<  occurence [2120] << endl;
		cout << "ocuurence 2151 :  " <<  occurence [2151] << endl;

}




void Frequency_extrapolation_test::
prepare_jack_nife_chain_data_test ()
{
//	Frequency_extrapolation	fre ( "W5_01234_d_pro_gly", FILL_UP_FREQUENCY_EXTRAPOLATION_STORE );
	Frequency_extrapolation	fre ( "W5_01234_dx2", COMMON_USAGE );

//	fre.prepare_united_chain_data	();
	fre.prepare_jack_nife_chain_data();

	fre.show_freq_data ("1A6M_");
	fre.show_freq_data ("_1A6M_");
	fre.show_freq_data ("together");

}
*********/

void Frequency_extrapolation_test::
show_freq_data_test ()
{
	Frequency_extrapolation	fre   ( "c30_W5_noDEG_ZIP",	COMMON_USAGE);
//	fre.prepare_together_freq_data();
//	Frequency_extrapolation	fre   ( "W5_01234_trivial_P_AL_VI",	COMMON_USAGE);

	//fre.show_freq_data("together");
}


void Frequency_extrapolation_test::
prepare_together_freq_data_test ()
{
//	Frequency_extrapolation	fre   ( "W5_01234_trivial_P_AL_VI",	COMMON_USAGE);	fre.show_freq_data_slow ();
//	fre.show_freq_data ("together");	fre.show_freq_data_short ("together");
	Frequency_extrapolation	fre  ( "c30_w9_tail",	COMMON_USAGE);	
	//fre.prepare_together_freq_data();

	fre.show_freq_data_slow ();

}

 
