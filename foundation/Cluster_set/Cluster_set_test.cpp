#pragma warning( disable : 4786 )
#pragma warning( disable : 4018 )

#include "Cluster_set_test.h"
#include "Cluster_set.h"

#include "Single_cluster_record.h"

#include "../Fragment_base/Chain_binary.h"
#include "../Fragment_base/Fragment_base_subtle.h" 

#include "../CommonFunc.h"
#include "../Censorship.h"

#include "../Geometry_util/Geometry_util.h"

#include <fstream>
#include <iostream>

#include "../Pair_int_double.h"

using namespace std;

extern ofstream log_stream;
extern Censorship configuration;


vector <vector <int> >  
analyase_local_structure_presence (
	vector < vector < double > >  & coord_in_cluster_system, 
	const double constrain,
	const int number_of_classes,
	const int fragment_length);

void print_local_result(
	string cluster_set_name,
	string &PDB_chain_ID,
	vector < vector < double > >  & coord_in_cluster_system,
	string & sequence, 
	vector <vector <int> > local_presence,
	const int fragment_length,
	const double constrain );

Cluster_set_test::
~Cluster_set_test()
{
	cout << "Cluster_set_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}
void Cluster_set_test::
prepare_cluster_picture_and_description ()
{
	string  cluster_set_name = string ("30_5a_20");

	Cluster_set cla_ob(
		cluster_set_name,
		COMMON_USAGE_CLUSTER_SET_MODE);

	string host_dir = cla_ob.get_host_dir() + string("analyze_cluster/");
	string analyse_protocol = host_dir + string("analise_protocol");

	ofstream out( analyse_protocol .c_str() );
	if ( ! out	)
	{	
		log_stream << "can't create " << analyse_protocol<< endl;
		cout       << "can't create " << analyse_protocol<< endl;
		exit (1);	
	}


	vector <int > choosen_indexes = cla_ob.get_claster_motif_index();

	for (int ii=0; ii<choosen_indexes.size(); ii++ )
	{
		cout << choosen_indexes [ii] << endl;
	}

	Fragment_base_subtle *fragment_base = cla_ob.get_fragment_base();


	int length = cla_ob.fragment_length();
//	double *cord_set_1 = new double [length *9*10];


	double *cord_set_buffer_1 = new double [length *9*100];
	double *cord_set_buffer_2 = new double [length *9*100];

	double *cord_set_2 = new double [length *9*10];
	vector <string> residue_names_2;				
	vector <string> in_chain_residue_number_2;


	{
			string pdb_chain_ID;
			string fragment_sequence;
			int serial_number;
			string pdb_resudue_number;

			fragment_base->fill_up_record_items ( 
				choosen_indexes [0] ,
				pdb_chain_ID,
				fragment_sequence,
				serial_number,
				pdb_resudue_number	);


			out << choosen_indexes [0] << "\t" ;
			out << pdb_chain_ID			<< "\t" ;
			out << fragment_sequence	<< "\t" ;
			out << serial_number		<< "\t" ;
			out << pdb_resudue_number	<< endl;

			
			
			Chain_binary *cb = new Chain_binary ( pdb_chain_ID);
			int start_pos = serial_number;
			
			
//			cb->extract_fragment (start_pos,length,cord_set_1 );
			vector <string> total_residue_names				= cb->get_residue_names () ;
			vector <string> total_in_chain_residue_number	= cb->get_in_chain_residue_number (); 

			
			for ( int kk=0;kk<length;kk++)
			{
				residue_names_2.push_back(total_residue_names[start_pos+kk]); 
				in_chain_residue_number_2.push_back(total_in_chain_residue_number[start_pos+kk]);
			}

			cb->extract_fragment (start_pos,length,cord_set_2 );

			origin_to_zero (cord_set_2,length*3);

			ostringstream ost ;
			int ii=0;
			ost << ii ;
			string s_number = ost.str();

//			string path_to_pdb = host_dir + s_number + string (".ent");
//			cb->save_pdb_fragment (	cord_set_2, length, path_to_pdb, residue_names_2, in_chain_residue_number_2);

	}

	
	string chain_ID_set = "!BCDEFGHIJKLMNOPQRSTUVWXYZ12345";

	for (int ii=1; ii<choosen_indexes.size(); ii++ )
	{
			string pdb_chain_ID;
			string fragment_sequence;
			int serial_number;
			string pdb_resudue_number;

			fragment_base->fill_up_record_items ( 
				choosen_indexes [ii] ,
				pdb_chain_ID,
				fragment_sequence,
				serial_number,
				pdb_resudue_number	);


			out << choosen_indexes [ii] << "\t" ;
			out << pdb_chain_ID			<< "\t" ;
			out << fragment_sequence	<< "\t" ;
			out << serial_number		<< "\t" ;
			out << pdb_resudue_number	<< endl;

			
			
			Chain_binary *cb = new Chain_binary ( pdb_chain_ID);
			int start_pos = serial_number;
			
			
//			cb->extract_fragment (start_pos,length,cord_set_1 );
			vector <string> total_residue_names				= cb->get_residue_names () ;
			vector <string> total_in_chain_residue_number	= cb->get_in_chain_residue_number (); 

			vector <string> residue_names;				
			vector <string> in_chain_residue_number;
			
			for ( int kk=0;kk<length;kk++)
			{
				residue_names.push_back(total_residue_names[start_pos+kk]); 
				in_chain_residue_number.push_back(total_in_chain_residue_number[start_pos+kk]);
			}

			double *cord_set_1 = new double [length *9*10];

			cb->extract_fragment (start_pos,length,cord_set_1 );

			ostringstream ost ;
			ost << ii ;
			string s_number = ost.str();

	//		string path_to_pdb = host_dir + s_number + string (".ent");
	//		cb->save_pdb_fragment (	cord_set_1, length, path_to_pdb, residue_names, in_chain_residue_number);


		bool is_ok_alignment = align_two_chains (
			cord_set_1,
			cord_set_2,
			length);


		double rotation_matrix[9];
		memset (rotation_matrix,0,9*sizeof(double));

		make_rotation_matrix_by_cartesian_triad (
			cord_set_2,	//double *c_1,
			cord_set_2+3,	//double *cur,
			cord_set_2+6,	//double *c_2,
			rotation_matrix);	//double *rotation_matrix_nb);

		memset(cord_set_buffer_1,0,9*length*sizeof(double));
		memset(cord_set_buffer_2,0,9*length*sizeof(double));

		for (int kk=0;kk<length*3;kk++)
		{
			Geometry_util::multiplication_row_by_3x3_by_transposed ( 
				cord_set_2+3*kk,//const double	a[3],
				rotation_matrix,//const double	b[9], 
				cord_set_buffer_2+3*kk); // double			r[3] );

			Geometry_util::multiplication_row_by_3x3_by_transposed ( 
				cord_set_1+3*kk,//const double	a[3],
				rotation_matrix,//const double	b[9], 
				cord_set_buffer_1+3*kk); // double			r[3] );


		}

// Withous rotation 
/*		string path_to_pdb_no_rotation_1 = host_dir +   s_number + string ("_no_rotation") + string (".ent");
		cb->save_pdb_fragment (	cord_set_1, length, path_to_pdb_no_rotation_1, residue_names, in_chain_residue_number);
		string path_to_pdb_0_no_rotation_2 = host_dir +   string ("0_no_rotation") + string (".ent");
		cb->save_pdb_fragment (	cord_set_2, length, path_to_pdb_0_no_rotation_2, residue_names_2, in_chain_residue_number_2);
*/
//  rotation 
		string path_to_pdb_rotation_1 = host_dir +   s_number + string (".pdb");
		cb->save_pdb_fragment (	cord_set_buffer_1, length, path_to_pdb_rotation_1, residue_names, in_chain_residue_number,chain_ID_set[ii]);
		if (ii==1)
		{
			string path_to_pdb_rotation_2  = host_dir +   string ("0") + string (".pdb");
			cb->save_pdb_fragment (	cord_set_buffer_2, length, path_to_pdb_rotation_2, residue_names_2, in_chain_residue_number_2,'A');
		}

	
		delete [] cord_set_1;

	}
	
	delete [] cord_set_buffer_1 ;
	delete [] cord_set_buffer_2 ;
	delete [] cord_set_2;

	/*
		fragment_base_->fill_up_record_items ( 
			global_index  ,
			pdb_chain_ID,
			fragment_sequence,
			serial_number,
			pdb_resudue_number	);
*/




}


void Cluster_set_test::
optimize_clasterization_test ()
{
		string cluster_set_name = "30_5a_20";

		Cluster_set ob(
			cluster_set_name,
			FILL_UP_MODEL_CLUSTER_SET_MODE);


		ob.optimize_clasterization();

}
void Cluster_set_test::
check_manually_settin_mode()
{
	string cluster_set_name = "Manulally_3";
	Cluster_set ob(
		cluster_set_name,
		COMMON_USAGE_CLUSTER_SET_MODE);


}

#include "../Fragment_base/accepted_chain_data.h"
#include "../Fragment_base/Chain_binary.h"

#define  MAX_STURCTURE_LENGTH 100

void Cluster_set_test::
analyse_analyse_setted_regular_structure_presence()
{

	double constrain = 0.8;
	string cluster_set_name = "Manulally_3";
	Cluster_set ob(
		cluster_set_name,
		COMMON_USAGE_CLUSTER_SET_MODE);

	int number_of_classes = ob.number_of_classes();
	int fragment_length = ob.fragment_length();


	double **claster_motif_coordinates = ob.get_claster_motif_coordinates();



	string binary_file_name ("accepted_chain_list.bin");
	vector <string>		accepted_chain_ID_list;
	vector <int>		accepted_chain_lenth;


	fill_up_accepted_chain_data ( 
		accepted_chain_ID_list, 
		accepted_chain_lenth,
		binary_file_name );


//accepted_chain_ID_list.resize(0);
//accepted_chain_ID_list.push_back("1ET1A");


	vector < vector < int > > global_presence;
	global_presence.resize(number_of_classes);
	for (int kk=0;kk<number_of_classes;kk++)
		global_presence[kk].resize(MAX_STURCTURE_LENGTH);

	

	for (int ii =0; ii < accepted_chain_ID_list.size(); ii++ )
	{
		cout << ii << " " << accepted_chain_ID_list[ii] << endl; 
		if ( accepted_chain_ID_list[ii]  == "3IOXA" ) 
			continue;

		Chain_binary cb( accepted_chain_ID_list[ii]);

		
		vector < vector < double > >  coord_in_cluster_system = cb.positioning_chain_by_clasters_set ( 
			claster_motif_coordinates,
			fragment_length,
			number_of_classes ) ;


		vector <vector <int> > local_presence =  analyase_local_structure_presence (
			coord_in_cluster_system, 
			constrain,
			number_of_classes,
			fragment_length);
		

		string sequence = cb.get_sequence();
		print_local_result(
			cluster_set_name,
			accepted_chain_ID_list[ii],
			coord_in_cluster_system,
			sequence, 
			local_presence,
			fragment_length,
			constrain );


		for (int tt=0;tt<number_of_classes ;tt++ )
		{
			for (int jj=0;jj<MAX_STURCTURE_LENGTH  ;jj++ )
				global_presence [tt][jj] += local_presence[tt][jj];
		}
	}

		string result_file_name = 		
			configuration.option_meaning("Path_to_Cluster_set") + 	cluster_set_name +
			string ("/occurece_protocol/") + 
			"FINAL_OCCURENCE" + string (".TXT") ;

		ofstream out( result_file_name .c_str() );
		if ( ! out	)
		{	
			log_stream << "result_file_name   can't create " << result_file_name<< endl;
			cout       << "result_file_name   can't create " << result_file_name<< endl;
			exit (1);	
		}

		for (int ii=0;ii<number_of_classes ;ii++)
		{
			for (int jj=0;jj<MAX_STURCTURE_LENGTH  ;jj++)
			{
				PutVa(global_presence [ii][jj],out,8,3,'l');
			}
			out << endl;
		}


}

void print_local_result(
	string cluster_set_name,
	string &PDB_chain_ID,
	vector < vector < double > >  & coord_in_cluster_system,
	string & sequence, 
	vector <vector <int> > local_presence,
	const int fragment_length,
	const double constrain )
{


	Chain_binary cb( PDB_chain_ID);  // нужно только в том месте где стыдно

	int shift = fragment_length/2;

		string result_file_name = 		
			configuration.option_meaning("Path_to_Cluster_set") + 	cluster_set_name +
			string ("/occurece_protocol/") + 
			PDB_chain_ID + string (".struct_occurence") ;

		ofstream out( result_file_name .c_str() );
		if ( ! out	)
		{	
			log_stream << "result_file_name   can't create " << result_file_name<< endl;
			cout       << "result_file_name   can't create " << result_file_name<< endl;
			exit (1);	
		}


		for ( int kk=0;kk<coord_in_cluster_system.size();kk++)
		{
			PutVa(sequence[kk+shift],out,3,1,'l');
			for (int ii=0;ii<local_presence.size();ii++)  // это то же что number_of_classesэ
			{
				if (coord_in_cluster_system [kk][ii] < constrain && coord_in_cluster_system [kk][ii] != -1 ) 
				{
					out << "Yes "  ;
// Вот тут начинается через жопу. Стыдно должно быть
					if (ii==1)
					{


					}

				}
				else 
					out << "No  "  ;

			}
			out << "   ";
			int test = local_presence.size();
			for (int ii=0;ii<test ;ii++)  // это то же что number_of_classesэ
			{
				PutVaDouble (coord_in_cluster_system[kk][ii],out,8,3,'l');
			}
			out << endl;

		}

		for (int ii=0;ii<local_presence.size();ii++)
		{
			for (int jj=0;jj<local_presence[0].size();jj++)
			{
				PutVa(local_presence[ii][jj],out,8,3,'l');
			}
			out << endl;
		}
}

// используется только здесь 

vector <vector <int> >  
analyase_local_structure_presence (
	vector < vector < double > >  & coord_in_cluster_system, 
	const double constrain,
	const int number_of_classes,
	const int fragment_length)
{
	vector < vector < int > > local_presence;
	local_presence.resize(number_of_classes);
	for (int kk=0;kk<number_of_classes;kk++)
		local_presence[kk].resize(MAX_STURCTURE_LENGTH);

	int current_len ;
	bool catch_flag;

	for (int ii=0;ii<number_of_classes;ii++)
	{
		catch_flag=false;

		for (int kk=0;kk<coord_in_cluster_system.size();kk++)
		{
			if (coord_in_cluster_system[kk][ii] == -1) 
			{
				catch_flag=false;
				continue;
			}

			if ( catch_flag == false )
			{
				current_len = 0;
				if (coord_in_cluster_system[kk][ii] < constrain )
				{
					catch_flag = true;
					current_len = fragment_length;
				}
				else 
					continue;
			}
			else 
			{
				if (coord_in_cluster_system[kk][ii] < constrain )
				{
					current_len ++;
				}
				else
				{
					local_presence [ii][current_len] ++;
					catch_flag = false ;
		
					current_len = 0;

				}
			}			

		}
		if (catch_flag) 
			local_presence [ii][current_len] ++;

	}
	
	return local_presence ;
}


void Cluster_set_test::
prepare_cluster_distance_matrix_test()
{
	//string cluster_set_name = "Manulally_3";
	string cluster_set_name = "30_5a_20";
	Cluster_set ob(
		cluster_set_name,
		COMMON_USAGE_CLUSTER_SET_MODE);

	vector < vector <double> >  cluster_distance_matrix = 
		ob.prepare_cluster_distance_matrix() ;

	string path_to_cluster_distance_matrix_file = configuration.option_meaning("Path_to_Cluster_set")  +  ob.get_cluster_set_mame () + string ("/") + 		string ("cluster_distance_matrix") ;


	ofstream output ( path_to_cluster_distance_matrix_file.c_str());
	if ( ! output )	
	{	
		log_stream << "path_to_cluster_distance_matrix_file: ERROR -  can't cluster_distance_matrix_file" << endl;
		cout       << "path_to_cluster_distance_matrix_file: ERROR -  can't cluster_distance_matrix_file" << endl;
		exit (1);	
	}


	output << cluster_distance_matrix .size() << endl;
	for (int tt=0;tt<cluster_distance_matrix .size();tt++)
	{
		for (int pp=0;pp<cluster_distance_matrix .size();pp++)
			PutVaDouble (cluster_distance_matrix[tt][pp],output ,8,3,'l');
		output << endl ;
	}

/// !!!!!!!!!!!!!!!!!!!!
	vector < Pair_int_double > row_indexed_metrics; 

	vector < vector < Pair_int_double > > matrix_indexed_metrics; 
	int ii;
	for (  int kk=0; kk<cluster_distance_matrix .size();kk++)
	{
		for (  ii=0;ii<cluster_distance_matrix .size();ii++)
			row_indexed_metrics.push_back ( Pair_int_double (ii,cluster_distance_matrix[kk][ii]) );
		sort (row_indexed_metrics.begin(),row_indexed_metrics.end() );
		matrix_indexed_metrics.push_back(row_indexed_metrics);
		row_indexed_metrics.resize(0);
	}
	output << "Indexed array (sorted by distance)" << endl;


	for (  int kk=0; kk<cluster_distance_matrix .size();kk++)
	{
		for (  ii=0;ii<cluster_distance_matrix .size();ii++)
		{
			PutVa (matrix_indexed_metrics[kk][ii].index(),output ,3,0,'l');
			output << '(';
			PutVaDouble (matrix_indexed_metrics[kk][ii].value(),output ,5,3,'l');
			output << ')';
		}
		output << endl;
	}

// SCORE MATRIX 

	int nu_clu = cluster_distance_matrix .size();

	vector < vector <double> >  score_matrix_1;
	score_matrix_1.resize(cluster_distance_matrix .size());
	for (int tt=0;tt<cluster_distance_matrix .size();tt++)
		score_matrix_1[tt].resize(cluster_distance_matrix .size() );

	for (int kk=0;kk<cluster_distance_matrix .size();kk++)
	{
		for (int ii=0;ii<cluster_distance_matrix .size();ii++)
		{
			int index = matrix_indexed_metrics[kk][ii].index();
			score_matrix_1[kk][index]= nu_clu - ii;
		}
	}


	for (  int kk=0; kk<cluster_distance_matrix .size();kk++)
	{
		for (  ii=0;ii<cluster_distance_matrix .size();ii++)
		{
			PutVa (score_matrix_1[kk][ii],output ,8,0,'l');

/*			output << '(';
			PutVaDouble (matrix_indexed_metrics[kk][ii].value(),output ,5,3,'l');
			output << ')';
*/
		}
		output << endl;
	}

	

}

void Cluster_set_test::
check_COMMON_USAGE_CLUSTER_SET_MODE_test ()
{
		string cluster_set_name = "30_5a_20";

		Cluster_set ob(
			cluster_set_name,
			COMMON_USAGE_CLUSTER_SET_MODE);


		vector <int> claster_motif_index = ob.get_claster_motif_index();

		int number_of_classes = ob.number_of_classes();
		int fragment_length  = ob.fragment_length();

		double **claster_motif_coordinates  = ob.get_claster_motif_coordinates();
		for (int ii=0;ii<number_of_classes;ii++)
		{
			for (int kk=0;kk<fragment_length*9;kk++)
				cout << claster_motif_coordinates[ii][kk] << '\t';
			cout << endl;
		}

		test_( "number_of_classes ",	number_of_classes  	== 30	);
		test_( "fragment_length ",		fragment_length  	== 5	);
}


void Cluster_set_test::
newlife_constructor_test()
{
		string cluster_set_name = "30_5a_debug";

		Cluster_set ob(
			cluster_set_name,
			FILL_UP_MODEL_CLUSTER_SET_MODE);

		//ob.optimize_clasterization();

		int waited_claster_number = 30;
		double upper_dist=		1.7;	
		double lower_dist=		0.8;
		double step_dist=      0.1;

		

		string metrics_mode("SQUARE_ROOT_INVERSE_DISTANCE");


		vector <int> choosen_indexes = ob.algorithm_2(
			upper_dist,				
			lower_dist,				
			step_dist,				
			waited_claster_number,
			metrics_mode);


			vector < Single_cluster_record >  claster_diversity;
			double  distance_scattering, quadrate_distance_scattering; 

			ob.regulate_choseen_index_subsample_version( 
				choosen_indexes,
				claster_diversity,	
				distance_scattering,
				quadrate_distance_scattering ) ;


			const string output_file_name = string("D:/Didona/Test/plain_cluster_show_test");
				ob.plain_claster_show_subsample_version ( 
					output_file_name,
					claster_diversity,
					distance_scattering,
					quadrate_distance_scattering);



	ofstream output ( "D:/Didona/Test/newlife_constructor_test");
	if ( ! output )	
	{	
		log_stream	<< "can't create " << endl;
		cout		 << "can't create " << endl;
		exit (1);	
	}

	for (int ii=0;ii<choosen_indexes.size();ii++)
		output << choosen_indexes[ii] << endl;

}
/*
void Cluster_set_test::
constructor_test ()
{

	string  fragment_base_subtle_name = string ("5a");
	int zip_factor = 50; 
	int shift = 3;

	Cluster_set cla_ob(
		fragment_base_subtle_name,
		zip_factor,
		shift,
		CLUSTER_SET_FILL_UP);

	double upper_dist = 30;
		
	double	lower_dist				= 1.5;			
	double	step_dist				=  .1;			
	int		waited_claster_number	= 50;

//	string metrics_mode = "SQUARED_INVERSE_DISTANCE";

//	string metrics_mode = "FOURTH_POWER_INVERSE_DISTANCE";

//	string metrics_mode = "PLAIN_INVERSE_DISTANCE";
	string metrics_mode = "SQUARED_INVERSE_DISTANCE";

//	string metrics_mode = "MULTIPLIED_PLAIN_SQUARED_INVERSE_DISTANCE";
//	string metrics_mode = "SQUARE_ROOT_INVERSE_DISTANCE";

//SQUARED_INVERSE_DISTANCE

	vector < int >  choosen_indexes = cla_ob.algorithm_2(
			upper_dist,
			lower_dist,
			step_dist,
			waited_claster_number,
			metrics_mode);

	vector < Single_cluster_record >  claster_diversity;
	double  distance_scattering, quadrate_distance_scattering; 

	cla_ob.regulate_choseen_index( 
		choosen_indexes,
		claster_diversity,
		distance_scattering,
		quadrate_distance_scattering) ;


	cla_ob.plain_claster_show (
		"TEST/Cluster_set_test",
		claster_diversity,
		distance_scattering,
		quadrate_distance_scattering) ;

}
*/
/*
void Cluster_set_test::
optimize_clasterization_test ()
{

	string  fragment_base_subtle_name = string ("5a");
	int zip_factor = 50; 
	int shift = 3;


	Cluster_set cla_ob(
		fragment_base_subtle_name,
		zip_factor,
		shift,
		CLUSTER_SET_FILL_UP);


	int waited_claser_number= 8;
	cla_ob.optimize_clasterization(waited_claser_number);

}





void Cluster_set_test::
mutual_distance_for_BS_show_test ()
{

	vector <int > choosen_indexes = pull_out_claster_origin_structure_list ( 
		"D:/Agony/Store/Sane_metrics/5a/50_3/30/protocol", 
		30 );

	for (int ii=0; ii<choosen_indexes.size(); ii++ )
	{
		cout << choosen_indexes [ii] << endl;
	}

	string  fragment_base_subtle_name = string ("5a");
	int zip_factor = 50; 
	int shift = 3;

	Cluster_set cla_ob(
		fragment_base_subtle_name,
		zip_factor,
		shift,
		CLUSTER_SET_COMMON_USAGE);



	ofstream output ( "D:/Agony/Store/Supplement/BS_mutual_distance.txt");
	if ( ! output )	
	{	
		log_stream << "fill_up_database(): ERROR -  can't create rejected_files.protocol" << endl;
		cout       << "fill_up_database(): ERROR -  can't create rejected_files.protocol" << endl;
		exit (1);	
	}


	vector < vector <double> > Cluster_distance_matrix;

	Cluster_distance_matrix.resize( choosen_indexes.size());
	for (int tt=0;tt<choosen_indexes.size();tt++)
		Cluster_distance_matrix[tt].resize( choosen_indexes.size() );
	int tt;
	for (tt=0;tt<choosen_indexes.size();tt++)
		Cluster_distance_matrix[tt][tt]=0;


	cla_ob.mutual_distance_for_BS_show( 
		choosen_indexes,
		Cluster_distance_matrix);

	output << "             ";
	for (int ii=0;ii<Cluster_distance_matrix.size();ii++)
		PutVa (ii,output,6,3 ,'l');	
	output << endl << "______________________________________________________________________________________________________________________________________________________________________________________________"<< endl;

	for (int ii=0;ii<Cluster_distance_matrix.size();ii++)
	{
		PutVa (ii,output,6,3 ,'l');		PutVa ("|",output,6,3 ,'l');	
		for (int jj=0;jj<Cluster_distance_matrix[0].size();jj++)
		{
			PutVaDouble (Cluster_distance_matrix[ii][jj],output,6,3 ,'l');	
		}
		output << endl;
	}
}
*/

void Cluster_set_test::
pull_out_claster_origin_structure_list_test ( )
{
	vector <int > choosen_indexes = pull_out_claster_origin_structure_list ( 
		"D:/Didona/Store/Cluster_set/DSSP_SET_30/protocol");

	for (int ii=0; ii<choosen_indexes.size(); ii++ )
	{
		cout << choosen_indexes [ii] << endl;
	}

	string  cluster_set_name = string ("DSSP_SET_30");
	Cluster_set cla_ob(
		cluster_set_name,
		COMMON_USAGE_CLUSTER_SET_MODE);

	vector < Single_cluster_record >  claster_diversity;
	double  distance_scattering, quadrate_distance_scattering; 

	cla_ob.regulate_choseen_index_for_whole_base( 
		claster_diversity,
		distance_scattering,
		quadrate_distance_scattering) ;

	cla_ob.subtle_claster_show_for_whole_base(
		"D:/Didona/Store/Cluster_set/DSSP_SET_30/subtle_clasters_show",
		claster_diversity);

}