

#pragma warning( disable : 4786 )

#include "Sane_metrics_test.h"
#include "Sane_metrics.h"


#include "../Fragment_base/Fragment_base_subtle.h"

#include "../by_Qmol/kabsch_stolen.h"

#include "../Fragment_base/Chain_binary.h" 

#include "../Cluster_set/Cluster_set.h" 
#include "../CommonFunc.h"


#include "../BioPolymerMechanics/foundation/Model.h"
#include "../BioPolymerMechanics/foundation/Core_iterator.h"
#include "../BioPolymerMechanics/foundation/Atom.h"


#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

extern ofstream log_stream;

void local_analyser (string &angles_pull, string &angles_pull_name,string output_path );

Sane_metrics_test::
~Sane_metrics_test()
{
	cout << "Sane_metrics_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}

void Sane_metrics_test::
test_subbtle_error()
{
	//	string  name = string ("5a");

// for debug only
	string  name = string ("5a_debug");
	


	Sane_metrics *sa_me = new  Sane_metrics   (  
		name,
		SANE_METRICS_COMMON_USAGE);




	int zip_factor	= sa_me->get_zip_factor();
	int shift		= sa_me->get_shift();




}
void Sane_metrics_test::
test_of_zip_integrity()
{

//	string  fragment_base_subtle_name = string ("5a");

	// for debug only
	string  name = string ("3a_20");

	
	Sane_metrics   ob3(  
		name,
		SANE_METRICS_FILL_UP) ;

}

/// Этот тест надог бы к готовым кластерам

void Sane_metrics_test::
distance_to_canonical_structure_from_BS()
{
	string output_path = "D:/Didona/Store/Supplement/Analyse_standard_confortmations/"; 

	string standard_alpha_helix_pull = "-60 -45 180 -60 -45 180 -60 -45 180 -60 -45 180 -60 -45 180";
	string alpha_name = "standard_alpha_helix";

	string standard_antiparallel_betta	= "-140 135 180 -140 135 180 -140 135 180 -140 135 180 -140 135 180 ";

	string standard_parallel_betta		= "-120 115 180 -120 115 180 -120 115 180 -120 115 180 -120 115 180 ";

	string helix_tri_dec = "-49 -26 180 -49 -26 180 -49 -26 180 -49 -26 180 -49 -26 180 ";


	local_analyser (standard_alpha_helix_pull, string("standard_alpha_helix"),output_path );
	local_analyser (standard_antiparallel_betta, string("antiparallel_betta"),output_path );
	local_analyser (standard_parallel_betta, string("parallel_betta"),output_path );
	local_analyser (helix_tri_dec , string("helix_tri_dec"),output_path );


}

//используется только для 

void local_analyser (string &angles_pull, string &angles_pull_name,string output_path )
{
	

	Model* model = new Model;

	int len_seq = 5;

	model->join_aminoacid("Gly",1);	
	for (int ii=0;ii<len_seq-1;ii++)
		model->join_aminoacid("Gly",0);



	vector < Atom * > core_atoms = model->get_core_atoms();
	Atom * initial_atom = core_atoms.front();

	vector < double > matrix; matrix.resize(9);
	matrix[0]=1;
	matrix[4]=1;
	matrix[8]=1;
	initial_atom->set_rotation_matrix(matrix);
	initial_atom->set_x(0);
	initial_atom->set_y(0);
	initial_atom->set_z(0);


	vector <double> phi, psi, ome;
	phi.resize(len_seq);
	psi.resize(len_seq);
	ome.resize(len_seq);;
	
	istringstream ist (angles_pull);

	for (int kk=0;kk<len_seq;kk++)
	{
		ist >>phi[kk]; 
		phi[kk] *= Pythagorean_Number () /180;

		ist >>psi[kk];	
		psi[kk]*= Pythagorean_Number () /180;

		ist >>ome[kk];	
		ome[kk]*= Pythagorean_Number () /180;

	}

	model->init_phi_psi_owega_backbone_set ();
	for (int ii=0; ii < len_seq; ii++ )
	{
		model->set_phi(ii,phi[ii]);
		model->set_psi(ii,psi[ii]);
		model->set_ome(ii,ome[ii]);
	}

	
	model->calc_cartesain_coordinates ( initial_atom );


	string path_to_ent_store = output_path + angles_pull_name + string(".ent");
	model->save_as(path_to_ent_store);


	double *coord = new double [ len_seq * 9 * sizeof (double) ];
	for (int ii=0; ii < len_seq; ii++ )
	{
		Atom *  ome_atom = model->ome_atom_set (ii) ; 
		Atom *  phi_atom = model->phi_atom_set (ii) ;
		Atom *  psi_atom = model->psi_atom_set (ii) ;
	
		coord [9*ii+0] = ome_atom->x();
		coord [9*ii+1] = ome_atom->y();
		coord [9*ii+2] = ome_atom->z();

		coord [9*ii+3] = phi_atom->x();
		coord [9*ii+4] = phi_atom->y();
		coord [9*ii+5] = phi_atom->z();

		coord [9*ii+6] = psi_atom->x();
		coord [9*ii+7] = phi_atom->y();
		coord [9*ii+8] = phi_atom->z();

	}

// достать из базы BS 
	int number_of_clasters = 30;
	//int fragment_length_ = len_seq;
	vector <int > claster_motif_index= pull_out_claster_origin_structure_list (
		"D:/Agony/Store/Sane_metrics/5a/50_3/30/protocol");


	//string fbs_file_name =  sheduler_->option_meaning("FRAGMENT_BASE_SUBTLE_SOURCE");

//	number_of_classes_ =  atoi ( sheduler_->option_meaning("NUMBER_OF_CLASSES").c_str() ) ;
	//number_of_clasters_  =  atoi ( sheduler_->option_meaning("NUMBER_OF_CLASSES").c_str() ) ;

	Fragment_base_subtle *fbs= new  Fragment_base_subtle ("5a" ,FRAGMENT_BASE_SUBTLE_COMMON_USAGE)  ;

	int fragment_length = fbs->get_fragment_length();

	int shift = fbs->get_fragment_length()/2;


	double **claster_motif_coordinates = new double* [number_of_clasters] ;	
	for ( int ii=0; ii < number_of_clasters; ii++  )
		claster_motif_coordinates[ii]  = new double [fragment_length*9];

	for (int kk=0;kk<number_of_clasters;kk++)
		fbs->get_coord( claster_motif_index [kk], claster_motif_coordinates[kk]);

	double *cu_cord_set  = new double [fragment_length*9];
	double *w	= new double [fragment_length *9];
	for (int kk=0;kk<(fragment_length *9);kk++)
		w[kk] = 1;

	vector <double> set_of_coordinate_in_clasters_system ; 
	set_of_coordinate_in_clasters_system.resize(number_of_clasters);

			for (int zz=0;zz<number_of_clasters;zz++)
			{

				bool error_flag= true;
				double distance  = kabsch_rmsd	(
										coord , 
										claster_motif_coordinates[zz], 
										w, 
										(3*fragment_length),
										error_flag);

				set_of_coordinate_in_clasters_system[zz] = distance ;
			}

	string path_to_explanation_store = output_path + angles_pull_name + string(".explanation");

	ofstream out ( path_to_explanation_store.c_str() );
	if ( ! out ) 	
	{
		out			<<	path_to_explanation_store  << " can't create " << endl;
		log_stream	<<	path_to_explanation_store  << " can't create " << endl;
		exit(-1);
	}

	cout << angles_pull_name << endl;
	cout << angles_pull << endl;
	for (unsigned kk=0;kk<set_of_coordinate_in_clasters_system.size(); kk++)
	{
		PutVa (kk, out, 6,2,'l');
		PutVaDouble(set_of_coordinate_in_clasters_system[kk],out,8,3,'l');
		out << endl;
	}


	delete [] cu_cord_set;		delete [] w;

	delete  fbs;
	delete [] coord;
// надоело 
// fillup backbone cordinated 


}


/* antiparallel
 string angles_pull = "–140° 135° 180 –140° 135° 180 –140° 135° 180 –140° 135° 180 –140° 135° 180";
parallel
string angles_pull = "–120° 115° 180 –120° 115° 180 –120° 115° 180 –120° 115° 180 –120° 115° 180 ) parallel


********/


/*
void Sane_metrics_test::
calibration_single_extenal_structure()
{
		int fragment_length_ = 5;

		double current_Phi,current_Psi,current_omega;

		current_Phi = -60;
		current_Psi=  -45; 	
		current_omega = 180;

		Model* model = new Model;

		//int len_seq = 5;
		model->join_aminoacid("Gly",1);	
		for (int ii=0;ii<fragment_length_-1;ii++)
			model->join_aminoacid("Gly",0);

		model->init_phi_psi_owega_backbone_set ();   // только раз использовать

		for (int ii=0; ii < fragment_length_; ii++ )
		{

			model->set_phi(ii,current_Phi*Pythagorean_Number ()/180);
			model->set_psi(ii,current_Psi*Pythagorean_Number ()/180);
			model->set_ome(ii,current_omega*Pythagorean_Number()/180);
		}


		model->calc_cartesain_coordinates ( initial_atom );

		double *probe_fragment  = new double [fragment_length_*9];

		for (int tt=0;tt<core_atoms.size(); tt++)
		{
			string pdb_atom_name = core_atoms[tt]->get_pdb_atom_name();
			if ( pdb_atom_name == "N" || pdb_atom_name == "CA" || pdb_atom_name == "C") {
				probe_fragment  [3*tt]=core_atoms[tt]->x();
				probe_fragment  [3*tt+1]=core_atoms[tt]->y();
				probe_fragment  [3*tt+2]=core_atoms[tt]->z();
				cout << probe_fragment  [3*tt] << " ";
				cout << probe_fragment  [3*tt+1] << " ";
				cout << probe_fragment  [3*tt+2] << endl;
			}
		}

	string  name = string ("5a_20");   // вообще не нужно здезь Sane_metrics. Только Fragment_base
	
	Sane_metrics *sa_me = new  Sane_metrics   (  
		name,
		SANE_METRICS_COMMON_USAGE);


	vector<int> neighbour_distribution; neighbour_distribution.resize(100);
	double  metric_1 =0;
	double  metric_2=0;;
	double  metric_3=0;;
	double  metric_4=0;;
	double  step=0.1;;

	sa_me->	calibration_single_extenal_structure ( 
		probe_fragment,
		neighbour_distribution,
		metric_1,
		metric_2,
		metric_3,
		metric_4,
		step)	;		// будем видимо 0.1 использовать 

	string path_to_explanation_store = "D:/Didona/Test/calibration_single_extenal_structure_test";

	ofstream out ( path_to_explanation_store.c_str() );
	if ( ! out ) 	
	{
		out			<<	path_to_explanation_store  << " can't create " << endl;
		log_stream	<<	path_to_explanation_store  << " can't create " << endl;
		exit(-1);
	}
	out <<		metric_1 << endl;
	out <<		metric_2 << endl;
	out <<		metric_3 << endl;
	out << 		metric_4 <<  endl;
	out << 		"**************"<< endl;

	for (int kk=0;kk< 100;kk++)
		out << neighbour_distribution[kk] 		 << endl;

	delete sa_me ;
	
}
*/
void Sane_metrics_test::multiple_calibrating_test ()
{


	string path_to_explanation_store = "D:/Didona/Test/Tetrapeptide_multiple_calibrating_test";

	ofstream out ( path_to_explanation_store.c_str() );
	if ( ! out ) 	
	{
		out			<<	path_to_explanation_store  << " can't create " << endl;
		log_stream	<<	path_to_explanation_store  << " can't create " << endl;
		exit(-1);
	}

		int fragment_length_ = 5;

		double current_Phi,current_Psi,current_omega;

		current_Phi ;
		current_Psi ;
		current_omega;

		Model* model = new Model;

		//int len_seq = 5;
		model->join_aminoacid("Gly",1);	
		for (int ii=0;ii<fragment_length_-1;ii++)
			model->join_aminoacid("Gly",0);

		model->init_phi_psi_owega_backbone_set ();   // только раз использовать

		vector < Atom * > core_atoms = model->get_core_atoms();
		Atom * initial_atom = core_atoms.front();
		vector < double > matrix; matrix.resize(9);
		matrix[0]=1;
		matrix[4]=1;
		matrix[8]=1;
		initial_atom->set_rotation_matrix(matrix);
		initial_atom->set_x(0);
		initial_atom->set_y(0);
		initial_atom->set_z(0);

		string  name = string ("5a_20");   // вообще не нужно здезь Sane_metrics. Только Fragment_base
	
		Sane_metrics *sa_me = new  Sane_metrics   (  
			name,
			SANE_METRICS_COMMON_USAGE);



		double *probe_fragment  = new double [fragment_length_*9];

		double Phi_start = -180; 
		double Phi_end	 =  180; 

		double Psi_start = -180; 
		double Psi_end	 =  180; 

		current_omega = 180;

		double step = 1;
		double calibration_step = 0.1;

		current_Phi = Phi_start; 

		while ( current_Phi < Phi_end )
		{
			current_Psi = Psi_start; 
			while ( current_Psi < Psi_end )
			{
				for (int ii=0; ii < fragment_length_; ii++ )
				{

					model->set_phi(ii,current_Phi*Pythagorean_Number ()/180);
					model->set_psi(ii,current_Psi*Pythagorean_Number ()/180);
					model->set_ome(ii,current_omega*Pythagorean_Number()/180);
				}

				model->calc_cartesain_coordinates ( initial_atom );

				for (unsigned tt=0;tt<core_atoms.size(); tt++)
				{
					string pdb_atom_name = core_atoms[tt]->get_pdb_atom_name();
					if ( pdb_atom_name == "N" || pdb_atom_name == "CA" || pdb_atom_name == "C") {
						probe_fragment  [3*tt]=core_atoms[tt]->x();
						probe_fragment  [3*tt+1]=core_atoms[tt]->y();
						probe_fragment  [3*tt+2]=core_atoms[tt]->z();
					}
				}

				vector<int> neighbour_distribution; neighbour_distribution.resize(100);
				for (unsigned zz=0;zz<neighbour_distribution.size();zz++)
				{
					neighbour_distribution[zz]=0;
				}

				double  metric_1=0;
				double  metric_2=0;;
				double  metric_3=0;;
				double  metric_4=0;;


				sa_me->	calibration_single_extenal_structure ( 
					probe_fragment,
					neighbour_distribution,
					metric_1,
					metric_2,
					metric_3,
					metric_4,
					calibration_step)	;		// будем видимо 0.1 использовать 


				PutVaDouble (current_Phi,out,8,0,'l');
				PutVaDouble (current_Psi,out,8,0,'l');

				out << ":\t"; 

				PutVaDouble (metric_1,out,15,3,'l');
				PutVaDouble (metric_2,out,15,3,'l');
				PutVaDouble (metric_3,out,15,3,'l');
				PutVaDouble (metric_4,out,15,3,'l');

				out << ":\t"; 

				for (unsigned zz=0;zz<neighbour_distribution.size();zz++)
				{
					PutVa (neighbour_distribution[zz],out,8,0,'l');
				}
				
				out << endl; 


				current_Psi += step;
			}
			current_Phi += step;
			cout << current_Phi << endl;

		}

		delete [] probe_fragment ;
}
void Sane_metrics_test::
prepare_rama_plot ()
{

	string path_to_explanation_store = "D:/Didona/Test/ramachandran_plot_5_res.txt";
	ifstream in ( path_to_explanation_store.c_str() );
	if ( ! in ) 	
	{
		cout		<<	path_to_explanation_store  << " can't create " << endl;
		log_stream	<<	path_to_explanation_store  << " can't create " << endl;
		exit(-1);
	}

	string res1 = "D:/Didona/Test/rama_metrics_for_graph.txt";
	ofstream out1 ( res1 .c_str() );
	if ( ! out1 ) 	
	{
		cout		<<	res1 << " can't create " << endl;
		log_stream	<<	res1 << " can't create " << endl;
		exit(-1);
	}



	vector <int> neighbour_distribution; neighbour_distribution.resize(100);

	string current_line;
	vector <vector <double>> Main_metrics;
	Main_metrics.resize(360);
	for (int kk=0;kk<360;kk++)
		Main_metrics[kk].resize(360);


	while( getline( in  , current_line, '\n' ) )
	{
		if (current_line.size() == 0) 
			continue;
		if (   current_line[0] == '/'  || 
			   current_line[0] == '#'  || 
			   current_line[0] == ' '  || 
			   current_line[0] == '\n' || 
			   current_line[0] == '\0')
			continue;

		string dummy;
		double phi_d,psi_d;
		double metric_1,metric_2,metric_3,metric_4;

		istringstream ist (current_line);
		ist >> phi_d >> psi_d>> dummy >> metric_1 ; // >> metric_2 >>metric_3 >> metric_4;  
//		for (int ii=0; ii<50; ii++ )
//			ist >> neighbour_distribution[]

		int phi = (int) phi_d;
		int psi = (int) psi_d;

		Main_metrics[phi+180][psi+180] = metric_1;
		
	}						 
	PutVa("          ",out1,10,0,'l'); out1 <<'\t';
	for (int kk=0;kk<360;kk+=2)
	{
		PutVa(kk-180,out1,10,0,'l');
		out1 <<'\t';
	}
	out1 << endl;

	for (int ii=0;ii<360;ii+=2)
	{
		PutVa(ii-180,out1,10,0,'l');out1 <<'\t';
		for (int jj=0;jj<360;jj+=2)
		{
			PutVaDouble(Main_metrics[ii][jj] ,out1,10,2,'l');
			out1 <<'\t';
		}
		out1 << endl;

	}

}