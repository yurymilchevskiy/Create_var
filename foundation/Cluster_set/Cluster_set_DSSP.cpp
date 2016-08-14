#pragma warning( disable : 4786 )
#pragma warning( disable : 4018 )

#include "Cluster_set_DSSP.h"

#include "../Censorship.h"
#include "../Sheduler.h"

#include "../Fragment_base/Fragment_base_subtle.h" 
#include "../Zip_iteator.h"
#include "../Sane_metrics/Sane_metrics_DSSP.h"
#include "../Pair_int_double.h"
#include "../by_Qmol/kabsch_stolen.h"

#include "../CommonFunc.h"

#include "../Fragment_base/fill_up_fragment_torsion_angles.h"

#include "Single_cluster_record.h"

#include <iostream>
#include <cmath>

#include <cassert>

#include "../BioPolymerMechanics/foundation/Model.h"
#include "../BioPolymerMechanics/foundation/Core_iterator.h"
#include "../BioPolymerMechanics/foundation/Atom.h"

extern Censorship configuration;
extern ofstream log_stream;


using namespace std;

#ifndef GET_ANGLES_FROM_TUNE_FILE_H
#define GET_ANGLES_FROM_TUNE_FILE_H

void	get_angles_from_tune_file (
		vector <vector <double > >  & Phi_set,
		vector <vector <double > >  & Psi_set,
		vector <vector <double > >  & Omega_set,
		vector <string>				& conformation_name,
		const string				& data_filename);

#endif

Cluster_set_DSSP::
~Cluster_set_DSSP()
{

	 if ( fragment_base_  )
		 delete fragment_base_;
	 
	 if ( sa_me_  )
		delete sa_me_ ;

	 if ( sheduler_  )
		delete sheduler_ ;

	if (subsample_fragments_) 
	{
 		for ( unsigned ii=0; ii < subsample_size_; ii++  )
			delete [] subsample_fragments_[ii];
		delete [] subsample_fragments_;
	}

	if ( claster_motif_coordinates_ ) 		
	{
		for (int kk=0; kk<number_of_classes_; kk++)			delete [] claster_motif_coordinates_[kk];
		delete [] claster_motif_coordinates_;
		claster_motif_coordinates_=0;
	}

	if ( all_BASE_fragments_ != 0 )  
	{
		for ( int ii=0; ii < total_fragment_number_; ii++  )
			delete [] all_BASE_fragments_ [ii];
		delete [] all_BASE_fragments_;
	}
	all_BASE_fragments_= 0;
	total_fragment_number_=0;

}

Cluster_set_DSSP::Cluster_set_DSSP (
	const string & cluster_set_name,
	const string dssp_word,
	Сluster_set_DSSP_operating_mode run_mode): 
	
	cluster_set_mame_ 					(cluster_set_name),
	dssp_word_							(dssp_word),
//	sane_metrics_mame_					(0),
//	zip_factor_							(0),
//	shift_								(0),
	fragment_base_						(0),
	sa_me_								(0),
	subsample_size_					(0),
	subsample_fragments_			(0),
	sheduler_							(0),
	claster_motif_coordinates_ (0),
	all_BASE_fragments_ (0),
	total_fragment_number_ (0)
{
	string path_to_current_schedule = 
		configuration.option_meaning("Path_to_Cluster_set")  +  
		string(cluster_set_name) + string ("/") + 		string ("sheduler") ;

	 host_dir_ = 	configuration.option_meaning("Path_to_Cluster_set")  +  
			cluster_set_mame_ + string ("/");

	sheduler_					= new Sheduler  (path_to_current_schedule);	

	int pointer=0;
	int shift;
	int zip_factor ;
	int size_of_fragment_base;

	string sane_metrics_mame ;
	string fragment_base_subtle_name ;

	string small_bases_rejected_list = host_dir_ + string ("small_bases_rejected_list") ;
	ofstream r_stream ( small_bases_rejected_list .c_str(),ios::app  );
	if ( ! r_stream ) 	
	{
	cout				<<	small_bases_rejected_list << " can't create " << endl;
	log_stream			<<	small_bases_rejected_list << " can't create " << endl;
	exit(-1); 
	}

	int minimal_size_of_fragment_base;


	int kkkk=0;
	int number_of_choosen_record = 0;

	switch ( run_mode  )
	{
		case COMMON_USAGE_CLUSTER_DSSP_SET_MODE: 
			init_claster_motif() ;
			break;
		case FILL_UP_MODEL_CLUSTER_SET_DSSP_MODE:

		sane_metrics_mame_ = sheduler_->option_meaning("SANE_METRICS_NAME");
		sa_me_ = new  Sane_metrics_DSSP  (  
			sane_metrics_mame_,  // в новой версии это дирректория где собраны все базы Sane_metrics_DSSP, соответствующие разным DSSP words
			dssp_word_,
			SANE_METRICS_COMMON_DSSP_USAGE);

		zip_it_ = sa_me_->get_Zip_iteator();
		subsample_size_ = zip_it_->get_step_number();  // subsample size


		fragment_base_subtle_name = sa_me_->get_fragment_base_subtle_name();
		 zip_factor = zip_it_->get_zip_factor();
		 shift	    = zip_it_->get_shift();

	//	 maх_sample_size_
	

		fragment_base_  = new Fragment_base_subtle ( fragment_base_subtle_name , FRAGMENT_BASE_SUBTLE_COMMON_USAGE);



		//size_of_fragment_base =  fragment_base_->get_total_fragment_number () ;  
		 // теперь size_of_fragment_base  из sa_me_

		 base_indexes_for_dssp_words_ = sa_me_->get_base_indexes_for_dssp_words ();

		 size_of_fragment_base = base_indexes_for_dssp_words_.size();

		 minimal_size_of_fragment_base= atoi (sheduler_->option_meaning("MINIMAL_SIZE_OF_FRAGMENT_BASE").c_str() );
		

		 if ( size_of_fragment_base < minimal_size_of_fragment_base  )
		 {
			 r_stream << dssp_word << " " << size_of_fragment_base << endl;;
			 cout << dssp_word << " too small base " << size_of_fragment_base << endl;
			 return;
		 }

		length_  =  fragment_base_->get_fragment_length();

       pointer=0;


	   /*
			zip_it_ = new Zip_iteator (
				size_of_fragment_base,
				zip_factor,
				shift);	
*/

			subsample_fragments_ = new double* [subsample_size_] ;	

			for ( int ii=0; ii < subsample_size_; ii++  )
				subsample_fragments_[ii]  = new double [length_*9];

		//	заменим цикл на итератор Zip_iteator)
			zip_it_->rewind_back();

			number_of_choosen_record = 0;
			kkkk=0;
			
			while (1)
			{
				// только тут преогбразовать локальный индекс zip_it_ в глобальный индекс инде
				int curr_rec_number_in_fragment_base = zip_it_->current();
				fragment_base_->get_coord ( 
					base_indexes_for_dssp_words_[ curr_rec_number_in_fragment_base] , 
					subsample_fragments_[number_of_choosen_record] );

				number_of_choosen_record++;

				if ( zip_it_->has_next() )	zip_it_->next();
				else 					break;

				kkkk++;
			} 


			 
/*
			while (1)
			{
				int kk = zip_it_->current();
				fragment_base_->get_coord ( kk, subsample_fragments_[pointer++] );

				if ( zip_it_->has_next() )	zip_it_->next();
				else 					break;
			}
*/

			 host_dir_ = 	configuration.option_meaning("Path_to_Cluster_set")  +  
				cluster_set_mame_ + string ("/");

			optimize_clasterization();


			break;
		default :
			log_stream	<< "Inadmissible run mode for Cluster_Set class" << endl;
			cout		<< "Inadmissible run mode for Cluster_Set class" << endl;
			exit (-1);
	}

}


vector <int>  Cluster_set_DSSP::algorithm_2 (
	double upper_dist,				
	double lower_dist,				
	double step_dist,				
	int waited_claster_number,
	string & metrics_mode)
{

	double *w	= new double [length_*9];
	for ( int kk=0;kk<(length_*9);kk++)		w[kk] = 1;



// PREPARING METRICS. IT MAY BE COMPLICATED METRICS ( COMBINATION OF METRIS )  
	double *metrics_for_ordering	= new double [ subsample_size_*sizeof (double) ];
	warm_up_metrics( metrics_for_ordering,metrics_mode );


// CREATION AND SORTING VECTOR THAT ENABLE POINT TO INDEX BY VALUE
	vector < Pair_int_double > indexed_metrics; 
	int ii;
	for (  ii=0;ii<subsample_size_;ii++)
		indexed_metrics.push_back ( Pair_int_double (ii,metrics_for_ordering[ii]) );
	sort (indexed_metrics.begin(),indexed_metrics.end() );

	vector < int >  choosen_indexes;
	vector < Pair_int_double >	spatial_metrics; 
	vector < int >				forbidden_indexes;

	// test only
		double large_val = indexed_metrics.back().value();
		double small_val   = indexed_metrics.front().value(); 
	//***

		int test_first_chosen = indexed_metrics.back().index();
	choosen_indexes.push_back	( indexed_metrics.back().index());
	forbidden_indexes.push_back	( indexed_metrics.back().index());

	while ( choosen_indexes.size() <  waited_claster_number  ) 
	{
		spatial_metrics.clear();  


		for (  ii=0 ; ii< subsample_size_; ii++ ) 
		{
			bool run_flag = true;
			for (int jj=0;jj<forbidden_indexes.size();jj++) 
			{
				if ( ii == forbidden_indexes[jj] )
				{
					run_flag = false;
					break;
				}
			}

			if (! run_flag) 
				continue;

		
			for (  int kk=0; kk< choosen_indexes.size();kk++ ) 
			{

				bool error_flag=true;;
				double current_distance = kabsch_rmsd (
					subsample_fragments_[ii], 
					subsample_fragments_[choosen_indexes[kk]], 
					w, 
					3*length_, 
					error_flag);


				if ( current_distance < lower_dist)				// если слишком (lower_dist) близко к существ. кластерам
				{
					forbidden_indexes.push_back ( ii );
					run_flag = false;
					break;
				}

			}


			double Spatial_metrics_value = 0;
			//double Spatial_metrics_value = 1;
			int kk;
			for (  kk=0; kk< choosen_indexes.size();kk++ ) 
			{

				bool error_flag=true;;
				double current_distance = kabsch_rmsd (
					subsample_fragments_[ii],
					subsample_fragments_[choosen_indexes[kk]], 
					w,
					3*length_,
					error_flag);

				Spatial_metrics_value += sqrt(current_distance)/100; ;
			//	Spatial_metrics_value *= sqrt(current_distance); ;

			}
			Spatial_metrics_value *= metrics_for_ordering[ii];
			spatial_metrics.push_back(Pair_int_double (ii,Spatial_metrics_value) ) ;


		}

		sort (spatial_metrics.begin(),spatial_metrics.end() );
		int pretender_index = spatial_metrics.back().index(); 

		bool addition_flag = true;

		for (  int kk=0; kk< choosen_indexes.size();kk++ ) 
		{

			bool error_flag=true;;
			double current_distance = kabsch_rmsd (
					subsample_fragments_[pretender_index],
					subsample_fragments_[choosen_indexes[kk]], 
					w,
					3*length_,
					error_flag);



			if ( current_distance < lower_dist)				// если слишком (lower_dist) близко к существ. кластерам
			{
				forbidden_indexes.push_back ( pretender_index );
				addition_flag = false;
				break;
			}

		}

		if ( addition_flag )
		{
			choosen_indexes.push_back   ( pretender_index );
			forbidden_indexes.push_back	( pretender_index );
		}

		cout << choosen_indexes.size()  << endl;

		if( forbidden_indexes.size () > (subsample_size_ -3)  ) 
		{
			cout << " All chances exhausted !!! " << endl;
			break;
		}


	}
	
	delete [] metrics_for_ordering;
	
	delete [] w;

	return choosen_indexes;

}



void Cluster_set_DSSP::
warm_up_metrics( 
				double *metrics,
				const string & metrics_mode)
{
	// PREPARING METRICS. IT MAY BE COMPLICATED METRICS ( COMBINATION OF METRIS )  
//	string metrics_mode = clasterization_option_->option_meaning("METRICS_MODE_FOR_ORDERING");

	if		( metrics_mode == "SQUARE_ROOT_INVERSE_DISTANCE") 
		sa_me_->fill_metrics ( 0, metrics );	
	else if	( metrics_mode == "PLAIN_INVERSE_DISTANCE") 
		sa_me_->fill_metrics ( 1, metrics );	
	else if	( metrics_mode == "SQUARED_INVERSE_DISTANCE") 
		sa_me_->fill_metrics ( 2, metrics );	
	else if	( metrics_mode == "CUBED_INVERSE_DISTANCE") 
		sa_me_->fill_metrics ( 3, metrics );	
	else if	( metrics_mode == "FOURTH_POWER_INVERSE_DISTANCE") 
		sa_me_->fill_metrics ( 4, metrics );	
	else if ( metrics_mode == "MULTIPLIED_PLAIN_SQUARED_INVERSE_DISTANCE"	) 
	{
//		int subsample_size_ = fragment_base_->get_subsample_size_();

		double *metrics_1   = new double [ subsample_size_*sizeof (double) ];
		double *metrics_2   = new double [ subsample_size_*sizeof (double) ];

		sa_me_->fill_metrics ( 1, metrics_1 );	
		sa_me_->fill_metrics ( 2, metrics_2 );	

		for ( int ii=0; ii< subsample_size_; ii++ )
			metrics[ii] = metrics_1[ii] *metrics_2[ii] ;

		delete [] metrics_1;
		delete [] metrics_2;
	}
	else 
	{
		cout		<< "Strange record int options file : " << metrics_mode << endl;
		log_stream	<< "Strange record int options file : " << metrics_mode << endl;
		exit (1);
	}
}
void Cluster_set_DSSP::
regulate_choseen_index_subsample_version_ALL_BASE (
	const vector <int> & choosen_indexes,
	vector < Single_cluster_record > & claster_diversity,
	double & distance_scattering,
	double & quadrate_distance_scattering) 

{
	claster_diversity.resize(0);
	distance_scattering = 0; 
	quadrate_distance_scattering = 0;  

	vector < int >		neighbour_numbers; 					neighbour_numbers.					resize	( choosen_indexes.size() );
	vector < double >	distances_sum_to_origin; 			distances_sum_to_origin.			resize	( choosen_indexes.size() );
	vector < double >	distances_sum_to_origin_squared ; 	distances_sum_to_origin_squared .	resize	( choosen_indexes.size() );

	vector < vector < Pair_int_double > > pair_neighbour_index_distances;  pair_neighbour_index_distances.resize	( choosen_indexes.size() );

	double *cur_fragment = new double [length_*9];

	double *bs_fragment_lcopy = new double [length_*9];

	double *w	= new double [length_*9];
	int kk;
	for (  kk=0;kk<(length_*9);kk++)		w[kk] = 1;

	for (  int ii=0;ii< total_fragment_number_ ; ii+=10 ) 
	{
		 vector < Pair_int_double > indexed_dist_store;  
	 
	//	 memcpy(cur_fragment ,subsample_fragments_ [ii],length_*9*sizeof(double));
		 memcpy(cur_fragment ,all_BASE_fragments_[ii],length_*9*sizeof(double));
		 

		 for (int kk=0; kk< choosen_indexes.size();kk++ )
		 {
			memcpy(bs_fragment_lcopy,subsample_fragments_ [choosen_indexes[kk]],length_*9*sizeof(double));

 			bool error_flag=true;;
			int ttt = choosen_indexes[kk];
			double current_distance = kabsch_rmsd (
					cur_fragment,
					bs_fragment_lcopy,
					w,
					3*length_,
					error_flag);

			indexed_dist_store.push_back( Pair_int_double (kk,current_distance) );
		 }

 		sort (indexed_dist_store.begin(),indexed_dist_store.end() );
		int			nearest_index 	=	indexed_dist_store.front().index() ;
		double 		nearest_distance =	indexed_dist_store.front().value() ;
		neighbour_numbers				[ nearest_index]	++;
		distances_sum_to_origin			[ nearest_index]	+=  nearest_distance;
		distances_sum_to_origin_squared [ nearest_index]	+=  nearest_distance*nearest_distance;

		distance_scattering								+=  nearest_distance;
		quadrate_distance_scattering						+=  nearest_distance*nearest_distance;

		// index in database and distance to origin
		pair_neighbour_index_distances	[ nearest_index].push_back( Pair_int_double (ii,nearest_distance) ) ;

		indexed_dist_store.clear();
	 }

	for ( kk=0; kk < choosen_indexes.size();kk++ )
		sort(	pair_neighbour_index_distances[kk].begin(),	pair_neighbour_index_distances[kk].end() );


/// тут собака порылась ************************	
	for ( kk=0; kk < choosen_indexes.size();kk++ )
		claster_diversity.push_back( Single_cluster_record (
			choosen_indexes[kk],
			neighbour_numbers[kk],
			distances_sum_to_origin[kk],
			distances_sum_to_origin_squared[kk],
			pair_neighbour_index_distances [kk]) );
	
	sort (claster_diversity.begin(),
		  claster_diversity.end() );


//	for ( int  ii=0; ii < choosen_indexes.size(); ii++  )
//		delete [] bs_fragments [ii];
//	delete [] bs_fragments ;
	delete []bs_fragment_lcopy;

	delete [] w;
}



void Cluster_set_DSSP::
regulate_choseen_index_subsample_version( // fills properties for subset made be zip_it_ selection
	const vector <int> & choosen_indexes,
	vector < Single_cluster_record > & claster_diversity,
	double & distance_scattering,
	double & quadrate_distance_scattering) 

{
	claster_diversity.resize(0);
	distance_scattering = 0; 
	quadrate_distance_scattering = 0;  

	vector < int >		neighbour_numbers; 					neighbour_numbers.					resize	( choosen_indexes.size() );
	vector < double >	distances_sum_to_origin; 			distances_sum_to_origin.			resize	( choosen_indexes.size() );
	vector < double >	distances_sum_to_origin_squared ; 	distances_sum_to_origin_squared .	resize	( choosen_indexes.size() );

	vector < vector < Pair_int_double > > pair_neighbour_index_distances;  pair_neighbour_index_distances.resize	( choosen_indexes.size() );

	double *cur_fragment = new double [length_*9];

	double *bs_fragment_lcopy = new double [length_*9];

	double *w	= new double [length_*9];
	int kk;
	for (  kk=0;kk<(length_*9);kk++)		w[kk] = 1;

	for (  int ii=0;ii< subsample_size_ ; ii++ ) 
	{
		 vector < Pair_int_double > indexed_dist_store;  
	 
		 memcpy(cur_fragment ,subsample_fragments_ [ii],length_*9*sizeof(double));

		 for (int kk=0; kk< choosen_indexes.size();kk++ )
		 {
			memcpy(bs_fragment_lcopy,subsample_fragments_ [choosen_indexes[kk]],length_*9*sizeof(double));

 			bool error_flag=true;;
			int ttt = choosen_indexes[kk];
			double current_distance = kabsch_rmsd (
					cur_fragment,
					bs_fragment_lcopy,
					w,
					3*length_,
					error_flag);

			indexed_dist_store.push_back( Pair_int_double (kk,current_distance) );
		 }

 		sort (indexed_dist_store.begin(),indexed_dist_store.end() );
		int			nearest_index 	=	indexed_dist_store.front().index() ;
		double 		nearest_distance =	indexed_dist_store.front().value() ;
		neighbour_numbers				[ nearest_index]	++;
		distances_sum_to_origin			[ nearest_index]	+=  nearest_distance;
		distances_sum_to_origin_squared [ nearest_index]	+=  nearest_distance*nearest_distance;

		distance_scattering								+=  nearest_distance;
		quadrate_distance_scattering						+=  nearest_distance*nearest_distance;

		// index in database and distance to origin
		pair_neighbour_index_distances	[ nearest_index].push_back( Pair_int_double (ii,nearest_distance) ) ;

		indexed_dist_store.clear();
	 }

	for ( kk=0; kk < choosen_indexes.size();kk++ )
		sort(	pair_neighbour_index_distances[kk].begin(),	pair_neighbour_index_distances[kk].end() );


/// тут собака порылась ************************	
	for ( kk=0; kk < choosen_indexes.size();kk++ )
		claster_diversity.push_back( Single_cluster_record (
			choosen_indexes[kk],
			neighbour_numbers[kk],
			distances_sum_to_origin[kk],
			distances_sum_to_origin_squared[kk],
			pair_neighbour_index_distances [kk]) );
	
	sort (claster_diversity.begin(),
		  claster_diversity.end() );


//	for ( int  ii=0; ii < choosen_indexes.size(); ii++  )
//		delete [] bs_fragments [ii];
//	delete [] bs_fragments ;
	delete []bs_fragment_lcopy;

	delete [] w;
}


void Cluster_set_DSSP::
regulate_choseen_index_for_whole_base( // fills properties for subset made be zip_it_ selection
//	const vector <int> & choosen_indexes,
	vector < Single_cluster_record > & claster_diversity,
	double & distance_scattering,
	double & quadrate_distance_scattering) 

{

// choosen_indexes - relating to total ( not reduced base) !!!
	string path_to_current_protocol = configuration.option_meaning("Path_to_Cluster_set")  +  string(cluster_set_mame_) + string ("/") + 	dssp_word_ + 	string (".protocol") ;
	vector <int > choosen_indexes = pull_out_claster_origin_structure_list (path_to_current_protocol );


	claster_diversity.resize(0);
	distance_scattering = 0; 
	quadrate_distance_scattering = 0;  

	vector < int >		neighbour_numbers; 					neighbour_numbers.					resize	( choosen_indexes.size() );
	vector < double >	distances_sum_to_origin; 			distances_sum_to_origin.			resize	( choosen_indexes.size() );
	vector < double >	distances_sum_to_origin_squared ; 	distances_sum_to_origin_squared .	resize	( choosen_indexes.size() );

	vector < vector < Pair_int_double > > pair_neighbour_index_distances;  pair_neighbour_index_distances.resize	( choosen_indexes.size() );


	double *cur_fragment = new double [length_*9];
	double *bs_fragment_lcopy = new double [length_*9];


	double *w	= new double [length_*9];
	for ( int kk=0;kk<(length_*9);kk++)		w[kk] = 1;

	int total_fragment_number = fragment_base_->get_total_fragment_number();



	for (  int ii=0;ii< total_fragment_number ; ii++ )   
//	for (  int ii=0;ii< 100000 ; ii++ )
	{
		 vector < Pair_int_double > indexed_dist_store;  

		 fragment_base_->get_coord ( ii , cur_fragment  );
		 //memcpy(cur_fragment ,subsample_fragments_ [ii],length_*9*sizeof(double));

		 for (int kk=0; kk< choosen_indexes.size();kk++ )
		 {
			//memcpy(bs_fragment_lcopy,subsample_fragments_ [choosen_indexes[kk]],length_*9*sizeof(double));
			 memcpy(bs_fragment_lcopy,claster_motif_coordinates_[kk],length_*9*sizeof(double));

 			bool error_flag=true;;
			int ttt = choosen_indexes[kk];
			double current_distance = kabsch_rmsd (
					cur_fragment,
					//bs_fragments [kk], 
					bs_fragment_lcopy,
					w,
					3*length_,
					error_flag);

			indexed_dist_store.push_back( Pair_int_double (kk,current_distance) );
		 }

 		sort (indexed_dist_store.begin(),indexed_dist_store.end() );
		int			nearest_index 	=	indexed_dist_store.front().index() ;
		double 		nearest_distance =	indexed_dist_store.front().value() ;
		neighbour_numbers				[ nearest_index]	++;
		distances_sum_to_origin			[ nearest_index]	+=  nearest_distance;
		distances_sum_to_origin_squared [ nearest_index]	+=  nearest_distance*nearest_distance;

		distance_scattering								+=  nearest_distance;
		quadrate_distance_scattering						+=  nearest_distance*nearest_distance;

		// index in database and distance to origin
		pair_neighbour_index_distances	[ nearest_index].push_back( Pair_int_double (ii,nearest_distance) ) ;

		indexed_dist_store.clear();
	 }

		for (int  kk=0; kk < choosen_indexes.size();kk++ )
			sort(	pair_neighbour_index_distances[kk].begin(),	pair_neighbour_index_distances[kk].end() );


// тут собака порылась ************************	
	for ( int kk=0; kk < choosen_indexes.size();kk++ )
		claster_diversity.push_back( Single_cluster_record (
			choosen_indexes[kk],
			neighbour_numbers[kk],
			distances_sum_to_origin[kk],
			distances_sum_to_origin_squared[kk],
			pair_neighbour_index_distances [kk]) );
	
	sort (claster_diversity.begin(),
		  claster_diversity.end() );


	delete []	cur_fragment ;
	delete []	bs_fragment_lcopy ;
	delete []	w;
}




void Cluster_set_DSSP::
plain_claster_show_subsample_version ( 
	const string & output_file_name,
	vector < Single_cluster_record > & claster_diversity,
	double & distance_scattering,
	double & quadrate_distance_scattering )
{

//	int number_of_records_	= fragment_base_->get_number_of_records_();
	//int length_				= fragment_base_->get_lenght (); 


	ofstream output ( output_file_name.c_str());
	if ( ! output )	
	{	
		log_stream << "fill_up_database(): ERROR -  can't create rejected_files.protocol" << endl;
		cout       << "fill_up_database(): ERROR -  can't create rejected_files.protocol" << endl;
		exit (1);	
	}

	bool is_it_final_dssp_words = false;
	map <int, string >  index_in_base_to_dssp_word;
	
	string index_to_name_file = host_dir_ + dssp_word_ + string (".and_names ");
	ifstream  index_to_name_stream( index_to_name_file .c_str());
	if ( ! index_to_name_stream)	
	{	
		log_stream << "index_to_name_file  ERROR -  can't create "  << index_to_name_file<< endl;
		cout       << "index_to_name_file  ERROR -  can't create "  << index_to_name_file<< endl;
	}
	else 
	{
		is_it_final_dssp_words = true;

		string word;
		int index_in_base;

		while (index_to_name_stream >> index_in_base >> word)
		{
			index_in_base_to_dssp_word [index_in_base] = word;
		}
	}



	// choosen_indexes - relating to total ( not reduced base) !!!


//	number_of_records__ = get_number_of_records ();

//	int total_fragment_number = fragment_base_->get_total_fragment_number();

	PutVa ( claster_diversity.size(),			output,8,1 ,'l');	

	if ( regulate_choseen_index_mode_ == "SUBSAMPLE") 
	{
		PutVaDouble (distance_scattering /subsample_size_ ,			output,10,3 ,'l');	
		PutVaDouble (quadrate_distance_scattering/subsample_size_ ,	output,10,3 ,'l');	
	}
	else if ( regulate_choseen_index_mode_ == "ALL_BASE") 
	{
		PutVaDouble (distance_scattering /total_fragment_number_ ,			output,10,3 ,'l');	
		PutVaDouble (quadrate_distance_scattering/total_fragment_number_ ,	output,10,3 ,'l');	
	}



	output << endl;
	output << "___________________________________________________________________________" << endl;


	double *w	= new double [length_*9];
	for ( int kk=0;kk<(length_*9);kk++)		w[kk] = 1;

	double *cur_fragment = new double [length_*9];
	double *cur_fragment_test = new double [length_*9];

	int claster_diversity_size_minus_1 = claster_diversity.size()-1;
	for  ( int ii= claster_diversity_size_minus_1; ii >= 0; ii-- )
	{

//		fragment_base_->get_coord ( claster_diversity[ii].in_data_base_index(),	cur_fragment );

//		int test_index = claster_diversity[ii].in_data_base_index();

		int pre_global_index = zip_it_->get_global_index (claster_diversity[ii].in_data_base_index() );
		int global_index = base_indexes_for_dssp_words_ [pre_global_index];

		memcpy(cur_fragment ,subsample_fragments_[claster_diversity[ii].in_data_base_index()],length_*9*sizeof(double));

		fragment_base_->get_coord ( global_index , cur_fragment_test );

		vector < double > torsion_set;
		fill_up_fragment_torsion_angles ( 
			cur_fragment,		
			length_,
			torsion_set, 
			'd');

		vector < double > torsion_set_test;
		fill_up_fragment_torsion_angles ( 
			cur_fragment_test ,		
			length_,
			torsion_set_test, 
			'd');



//		string sequence				= cur_fragment->get_sequence();
//		string  pdb_chain_ID		= cur_fragment->get_pdb_chain_ID();

		string	 pdb_chain_ID;
		string	 fragment_sequence;
		int		 serial_number;
		string   pdb_resudue_number;
		char *record_ch	 = new char [ 1000* sizeof (char)] ;

	//	int pre_global_index = zip_it_->get_global_index (claster_diversity[ii].in_data_base_index() );

		//int global_index = base_indexes_for_dssp_words_ [pre_global_index];

		fragment_base_->fill_up_record_items ( 
			global_index  ,
			pdb_chain_ID,
			fragment_sequence,
			serial_number,
			pdb_resudue_number	);


		PutVa (  global_index ,			output,10,	5 ,'l');	

		if (is_it_final_dssp_words ) 
			PutVa ( index_in_base_to_dssp_word[global_index]	,			output,6,	5 ,'l');	

		PutVa ( pdb_chain_ID	,			output,6,	5 ,'l');	
		PutVa (fragment_sequence,								output,10,	1 ,'l');	
		PutVa (serial_number,								output,10,	1 ,'l');	
		PutVa (pdb_resudue_number,								output,10,	1 ,'l');	

		PutVa (claster_diversity[ii].neighbour_number(),						output,7,	1 ,'l');	
		PutVaDouble ( 100* ( ( double) claster_diversity[ii].neighbour_number() ) /subsample_size_,	output,10, 3,'l');	

		PutVaDouble ( claster_diversity[ii].distance_sum_to_origin()  /claster_diversity[ii].neighbour_number(),	output,10,3 ,'l');	
		PutVaDouble ( claster_diversity[ii].distance_sum_to_origin_squared()  /claster_diversity[ii].neighbour_number(),	output,10,3 ,'l');	

		for ( int kk =0; kk<torsion_set.size(); kk++ ) 
			PutVaDouble ( torsion_set_test[kk],output,10,1 ,'r');	


		output << endl;

	}
// *********************** Cluster_distance_matrix ********************

	vector < vector <double> > Cluster_distance_matrix;

	Cluster_distance_matrix.resize( claster_diversity.size());
	for (int tt=0;tt<claster_diversity.size();tt++)
		Cluster_distance_matrix[tt].resize( claster_diversity.size() );
	int tt;
	for (tt=0;tt<claster_diversity.size();tt++)
		Cluster_distance_matrix[tt][tt]=0;

	double *coord_1 = new double [ length_ * 9 * sizeof (double) ];
	double *coord_2 = new double [ length_ * 9 * sizeof (double) ];

	int counter = 0;
	int Class_SIZE= claster_diversity.size();
	for  ( int ii=  Class_SIZE-1 ; ii >= 0; ii-- )
	{
		int index1 =  claster_diversity[ii].in_data_base_index() ;

	//	int index1 =  claster_diversity[jj].in_data_base_index() ;
		int  pre_global_index_1= zip_it_->get_global_index ( index1);
		int global_index_1 = base_indexes_for_dssp_words_ [pre_global_index_1];
		fragment_base_->get_coord ( global_index_1 , coord_1 );


		PutVa (ii,		output,9,	1 ,'l');	
		PutVa (index1,		output,10,	2 ,'l');output << "|" ; 	

		for (int  kk =0 ; kk < counter; kk++ ) 
			output << "          ";
		counter++;


		for  ( int jj= ii ; jj >= 0; jj-- )
		{
			if (ii==jj) 
			{
				PutVa (0,		output,10,	3 ,'l');	
				continue;
			}


		int index2 =  claster_diversity[jj].in_data_base_index() ;
		int  pre_global_index_2= zip_it_->get_global_index ( index2);
		int global_index_2 = base_indexes_for_dssp_words_ [pre_global_index_2];
		fragment_base_->get_coord ( global_index_2 , coord_2 );





	//		memcpy(coord_1,subsample_fragments_[index1],length_*9*sizeof(double));
		//	memcpy(coord_2,subsample_fragments_[index2],length_*9*sizeof(double));


//			fragment_base_->get_coord(index1,coord_1);
//			fragment_base_->get_coord(index2,coord_2);




			bool error_flag=true;;
			double current_distance = kabsch_rmsd (
					coord_1,
					coord_2, 
					w,
					3*length_,
					error_flag);


			PutVa (current_distance,		output,10,	3 ,'l');
			
		//	int aa = ii - Class_SIZE + 1;
		//	int bb = jj - Class_SIZE + 1;
			int aa = ii;
			int bb = jj;

			Cluster_distance_matrix[aa][bb] = current_distance;
			Cluster_distance_matrix[bb][aa] = current_distance;

			cout <<  aa << "  " << bb << " " << Cluster_distance_matrix[aa][bb];
	
		}
		output << endl;
	}

//	ofstream cdm ( "Cluster_distance_matrix");
//	if ( ! cdm  )	
//	{	
//		log_stream << "can't create Cluster_distance_matrix" << endl;
//		cout       << "can't create Cluster_distance_matrix" << endl;
//		exit (1);	
//	}

	output << claster_diversity.size() << endl;
	for ( tt=0;tt<claster_diversity.size();tt++)
	{
		for (int pp=0;pp<claster_diversity.size();pp++)
			PutVaDouble (Cluster_distance_matrix[tt][pp],output ,12,8,'l');
		output << endl ;
	}

	
	delete [] cur_fragment ;
	delete [] cur_fragment_test ;


}

vector < vector <double> >  Cluster_set_DSSP::  // дублируется часть кода из plain_claster_show_subsample_version ну извитите. 
prepare_cluster_distance_matrix()  
{
		// *********************** Cluster_distance_matrix ********************

	vector < vector <double> > Cluster_distance_matrix;

	Cluster_distance_matrix.resize( number_of_classes_);
	for (int tt=0;tt<number_of_classes_;tt++)
		Cluster_distance_matrix[tt].resize( number_of_classes_ );
	int tt;
	for (tt=0;tt<number_of_classes_;tt++)
		Cluster_distance_matrix[tt][tt]=0;


//	string manually_setted_angles_flag = sheduler_->option_meaning("MANUALLY_SETTED_ANGLES_FILE");
//	if (manually_setted_angles_flag  != "UNKNOWN")  // все как было раньше


	double *coord_1 = new double [ fragment_length_ * 9 * sizeof (double) ];
	double *coord_2 = new double [ fragment_length_ * 9 * sizeof (double) ];

	int counter = 0;
	int Class_SIZE= number_of_classes_;


	double *w	= new double [fragment_length_*9];
	
	for ( int kk=0;kk<(fragment_length_*9);kk++)		w[kk] = 1;

	for  ( int ii=  Class_SIZE-1 ; ii >= 0; ii-- )
	{
		int index1 =  ii ;

		PutVa (ii,		log_stream,9,	1 ,'l');	
		PutVa (index1,		log_stream,10,	2 ,'l');log_stream << "|" ; 	

		for (int  kk =0 ; kk < counter; kk++ ) 
			log_stream << "          ";
		counter++;


		for  ( int jj= ii ; jj >= 0; jj-- )
		{
			if (ii==jj) 
			{
				PutVa (0,		log_stream,10,	3 ,'l');	
				continue;
			}


			int index2 =  jj ;


			memcpy(coord_1,claster_motif_coordinates_[index1],fragment_length_*9*sizeof(double));
			memcpy(coord_2,claster_motif_coordinates_[index2],fragment_length_*9*sizeof(double));


//			fragment_base_->get_coord(index1,coord_1);
//			fragment_base_->get_coord(index2,coord_2);




			bool error_flag=true;;
			double current_distance = kabsch_rmsd (
					coord_1,
					coord_2, 
					w,
					3*fragment_length_,
					error_flag);


			PutVa (current_distance,		log_stream,10,	3 ,'l');
			
		//	int aa = ii - Class_SIZE + 1;
		//	int bb = jj - Class_SIZE + 1;
			int aa = ii;
			int bb = jj;

			Cluster_distance_matrix[aa][bb] = current_distance;
			Cluster_distance_matrix[bb][aa] = current_distance;

			cout <<  aa << "  " << bb << " " << Cluster_distance_matrix[aa][bb];
	
		}
		log_stream << endl;
	}

//	ofstream cdm ( "Cluster_distance_matrix");
//	if ( ! cdm  )	
//	{	
//		log_stream << "can't create Cluster_distance_matrix" << endl;
//		cout       << "can't create Cluster_distance_matrix" << endl;
//		exit (1);	
//	}

	log_stream << number_of_classes_ << endl;
	for ( tt=0;tt<number_of_classes_;tt++)
	{
		for (int pp=0;pp<number_of_classes_;pp++)
			PutVaDouble (Cluster_distance_matrix[tt][pp],log_stream ,12,8,'l');
		log_stream << endl ;
	}

	

	delete[] w;


	return Cluster_distance_matrix;


}



void Cluster_set_DSSP::
optimize_clasterization () 
{

	// пока task_file_name одинаковый для всех dssp words
	string task_file_name				= host_dir_ + 	sheduler_->option_meaning("OPTIMIZATION_PARAMETERS_FILE");  
	string protocol_file_name			= host_dir_ + dssp_word_ + string(".")	+ sheduler_->option_meaning("PROTOCOL_FILE_NAME") ;
	string plain_claster_show_file_name = host_dir_ + dssp_word_ + string(".")	+ 	sheduler_->option_meaning("PLAIN_CLASTER_SHOW_FILE_NAME");
		
	Sheduler *optimize_clasterization_option	= new Sheduler  ( task_file_name  );

	int		waited_cluster_number			= atoi ( optimize_clasterization_option->option_meaning("WAITED_CLUSTER_NUMBER").		c_str() ) ;			
	double	upper_distance					= atof ( optimize_clasterization_option->option_meaning("UPPER_DISTANCE").		c_str() ) ;			
	double	lower_distance_max				= atof ( optimize_clasterization_option->option_meaning("LOWER_DISTANCE_MAX").	c_str() ) ;						
	double	lower_distance_min				= atof ( optimize_clasterization_option->option_meaning("LOWER_DISTANCE_MIN").	c_str() ) ;						
	double	distance_step					= atof ( optimize_clasterization_option->option_meaning("DISTANCE_STEP").		c_str() ) ;						
	string 	objective_parameter				= optimize_clasterization_option->option_meaning("OBJECTIVE_PARAMETER");
	string  metrics_mode_mode				= optimize_clasterization_option->option_meaning("METRICS_MODE_MODE");

	regulate_choseen_index_mode_	=  optimize_clasterization_option->option_meaning("REGULATE_CHOSEEN_INDEX_MODE");

	vector <string> metrics_mode_set;

	if ( metrics_mode_mode == "ALL")
	{
		metrics_mode_set.push_back("SQUARE_ROOT_INVERSE_DISTANCE");
		metrics_mode_set.push_back("PLAIN_INVERSE_DISTANCE");
		metrics_mode_set.push_back("SQUARED_INVERSE_DISTANCE");
		metrics_mode_set.push_back("MULTIPLIED_PLAIN_SQUARED_INVERSE_DISTANCE");
	}
	else
		metrics_mode_set.push_back(metrics_mode_mode);


	int counter_first = 0;
	double	distance_scattering_min				= 999999999999;
	double	quadrate_distance_scattering_min	= 999999999999;

	bool is_exist_analysis_yet = assign_best_parameters (
										&distance_scattering_min,
										&quadrate_distance_scattering_min,
										protocol_file_name );

	ofstream p_stream ( protocol_file_name .c_str(),ios::app  );
	if ( ! p_stream ) 	
	{
		p_stream			<<	protocol_file_name << " can't create " << endl;
		log_stream			<<	protocol_file_name << " can't create " << endl;
		exit(-1); 
	}
	 

	int step_number = ( lower_distance_max - lower_distance_min )/distance_step;

	bool improvement_flag = false;
	int counter = 0 ;
	for (int ii=0;ii < step_number; ii++)
	{
		double cu_lower_distance = lower_distance_min + ii * distance_step;

		for (int jj=0; jj<metrics_mode_set.size(); jj++)
		{

			improvement_flag = false;

			vector < int >  choosen_indexes = algorithm_2 (
				upper_distance,
				cu_lower_distance,
				distance_step,
				waited_cluster_number,
				metrics_mode_set[jj] );

			if ( waited_cluster_number != choosen_indexes.size() )
				continue;

			vector < Single_cluster_record >  claster_diversity;
			double  distance_scattering, quadrate_distance_scattering; 

			if ( regulate_choseen_index_mode_ == "SUBSAMPLE") 
				regulate_choseen_index_subsample_version ( 
					choosen_indexes,
					claster_diversity,	
					distance_scattering,
					quadrate_distance_scattering ) ;
			else if ( regulate_choseen_index_mode_ == "ALL_BASE") 
			{
				if (total_fragment_number_ == 0)
				{
					total_fragment_number_ = fragment_base_->get_total_fragment_number();
					all_BASE_fragments_ = new double*       [total_fragment_number_] ;	
					memset ( all_BASE_fragments_,0,sizeof(double*)	 *total_fragment_number_);
					for ( int ii=0; ii < total_fragment_number_; ii++  )
							 all_BASE_fragments_[ii]  = new double [length_*9];
					for (int ttt=0;ttt<total_fragment_number_;ttt++)
						fragment_base_->get_coord ( 
							ttt , all_BASE_fragments_[ttt] );
				}

				regulate_choseen_index_subsample_version_ALL_BASE ( 
					choosen_indexes,
					claster_diversity,	
					distance_scattering,
					quadrate_distance_scattering ) ;

			}
			else 
			{
				cout		<< "strange regulate_choseen_index_mode_ mode : " << regulate_choseen_index_mode_ << endl;
				log_stream  << "strange regulate_choseen_index_mode_ mode : " << regulate_choseen_index_mode_ << endl;
			}

			cout << quadrate_distance_scattering	/	subsample_size_ << endl; 

			if ( (quadrate_distance_scattering < quadrate_distance_scattering_min) || counter_first  == 0  )
			{
				counter_first  ++;
				improvement_flag = true;
				quadrate_distance_scattering_min = quadrate_distance_scattering ;
//				p_stream << distance_scattering << "  " << quadrate_distance_scattering_min << "  ";

				
				PutVaDouble (cu_lower_distance,		p_stream,10,3 ,'l');	
				PutVa (metrics_mode_set[jj],p_stream,35,30,'l');
				p_stream << "  ";

				PutVaDouble (distance_scattering				/	subsample_size_,		p_stream,10,3 ,'l');	
				PutVaDouble (quadrate_distance_scattering_min	/	subsample_size_,		p_stream,10,3 ,'l');	

				assert ( claster_diversity.size() ==  waited_cluster_number );

				PutVa ( waited_cluster_number , p_stream,5,	3 ,'l');
				
//				for (int kk=0; kk < waited_cluster_number; kk++	)
				for  ( int kk=claster_diversity.size()-1 ; kk >= 0; kk-- )
				{

					int pre_global_index = zip_it_->get_global_index (claster_diversity[kk].in_data_base_index() );
					int global_index = base_indexes_for_dssp_words_ [pre_global_index];

					PutVa ( global_index,p_stream,10,5,'l');	

//					PutVa ( zip_it_->get_global_index( claster_diversity[kk].in_data_base_index() ),			p_stream,10,	5 ,'l');	
				}
				p_stream << endl;

				plain_claster_show_subsample_version ( 
					plain_claster_show_file_name,
					claster_diversity,
					distance_scattering,
					quadrate_distance_scattering_min );

			}




		}
	}

	delete  optimize_clasterization_option;
}
bool Cluster_set_DSSP::
assign_best_parameters (
	double *distance_scattering_min,
	double *quadrate_distance_scattering_min,
	const string & protocol_file_name  )
{
	*distance_scattering_min			= 999999;
	*quadrate_distance_scattering_min	= 999999;

	ifstream p_stream ( protocol_file_name.c_str());

	if ( ! p_stream )	{	
		log_stream	 << "ERROR -  can't find binary file" << protocol_file_name<< endl;
		cout		 << "ERROR -  can't find binary file" << protocol_file_name<< endl;
		return false;	
	}



	string current_line;
	while( getline( p_stream  , current_line, '\n' ) )
	{
		if (   current_line[0] == '/'  || 
			   current_line[0] == '#'  || 
			   current_line[0] == ' '  || 
			   current_line[0] == '\n' || 
			   current_line[0] == '\0')
			continue;

		double current_distance_scattering, current_quadrate_distance_scattering;

		{
			string word; 

			istringstream ist (current_line);
			ist >>  word;
			ist >>  word;

			ist >>  current_distance_scattering >> current_quadrate_distance_scattering ;
			if ( current_quadrate_distance_scattering < (*quadrate_distance_scattering_min) )
			{
				*quadrate_distance_scattering_min	= current_quadrate_distance_scattering;
				*distance_scattering_min			= current_distance_scattering;
			}
		}
	}
	*quadrate_distance_scattering_min	= (*quadrate_distance_scattering_min)	* subsample_size_;
	*distance_scattering_min			= (*distance_scattering_min)  			* subsample_size_;

	return true;
}

void Cluster_set_DSSP:: 
init_claster_motif() 
{
//	string claster_origin_structure_file =  sheduler_->option_meaning("CLASTER_ORIGIN_STRUCTURE_LIST");


	string manually_setted_angles_flag = sheduler_->option_meaning("MANUALLY_SETTED_ANGLES_FILE");
	if (manually_setted_angles_flag  == "UNKNOWN")  // все как было раньше
	{
	
		string path_to_cluster_set_protocol_file_name = 
			configuration.option_meaning("Path_to_Cluster_set")  +  
			cluster_set_mame_ + string ("/") + string ("protocol");


		claster_motif_index_ = pull_out_claster_origin_structure_list (
									path_to_cluster_set_protocol_file_name);

		number_of_classes_ =  claster_motif_index_.size(); 

// ***
	string sane_metrics_mame = sheduler_->option_meaning("SANE_METRICS_NAME");
	

	sa_me_ = new  Sane_metrics_DSSP (  
		sane_metrics_mame,
		dssp_word_,
		SANE_METRICS_COMMON_DSSP_USAGE);

	string fragment_base_subtle_name = sa_me_->get_fragment_base_subtle_name();
	int zip_factor = sa_me_->get_zip_factor();
	int shift	   = sa_me_->get_shift();
	
	length_	   = sa_me_->get_length();
	
	fragment_base_  = new Fragment_base_subtle ( fragment_base_subtle_name , FRAGMENT_BASE_SUBTLE_COMMON_USAGE);
// ****

		fragment_length_ = fragment_base_->get_fragment_length();
	
		claster_motif_coordinates_ = new double* [number_of_classes_] ;	
		for ( int ii=0; ii < number_of_classes_; ii++  )
			claster_motif_coordinates_[ii]  = new double [fragment_length_*9];


		for (int kk=0;kk<number_of_classes_;kk++)
			fragment_base_->get_coord( claster_motif_index_ [kk], claster_motif_coordinates_[kk]);

	}
	else   // тут кластеры задаются вручную по двугранным углам заданным в соотв. файле
	{
		
		string test = sheduler_->option_meaning("FRAGNMENT_LENGTH");

		fragment_length_	= atoi(sheduler_->option_meaning("FRAGNMENT_LENGTH").c_str() ) ;
		number_of_classes_	= atoi(sheduler_->option_meaning("NUMBER_OF_CLASSES").c_str() ) ;

		// FIX не забыть учесть dssp_word
		string path_to_manually_setted_angles = configuration.option_meaning("Path_to_Cluster_set")  +  string(cluster_set_mame_) + string ("/") + 		string (manually_setted_angles_flag) ;

		Model* model = new Model;

		//int len_seq = 5;
		model->join_aminoacid("Gly",1);	
		for (int ii=0;ii<fragment_length_-1;ii++)
			model->join_aminoacid("Gly",0);

		vector <vector <double> >	Phi_set;
		vector <vector <double> >	Psi_set;
		vector <vector <double> >	Omega_set;
		vector <string> conformation_name;

		get_angles_from_tune_file (
			Phi_set,
			Psi_set,
			Omega_set,
			conformation_name,
			path_to_manually_setted_angles);

	//	vector < Atom * > core_atoms = model->get_core_atoms();
	//	Atom * initial_atom = core_atoms.front();


		vector <double> phi, psi, ome;
		phi.resize(fragment_length_);
		psi.resize(fragment_length_);
		ome.resize(fragment_length_);;

		model->init_phi_psi_owega_backbone_set ();   // только раз использовать


		int check_number_of_classes = Phi_set.size();

		assert (check_number_of_classes == number_of_classes_ );

		claster_motif_coordinates_ = new double* [number_of_classes_] ;	
		for ( int ii=0; ii < number_of_classes_; ii++  )
			claster_motif_coordinates_[ii]  = new double [fragment_length_*9];


	
		for (int kk=0;kk<Phi_set.size(); kk++)
		{
			for (int ii=0; ii < fragment_length_; ii++ )
			{

				double phi = Phi_set[kk][ii]; 
				double psi = Psi_set[kk][ii];
				double ome = Omega_set[kk][ii];

				model->set_phi(ii,Phi_set[kk][ii]*Pythagorean_Number ()/180);
				model->set_psi(ii,Psi_set[kk][ii]*Pythagorean_Number ()/180);
				model->set_ome(ii,Omega_set[kk][ii]*Pythagorean_Number()/180);
			}

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

			model->calc_cartesain_coordinates ( initial_atom );

			model->save_as("aaaa.ent");

			for (int tt=0;tt<core_atoms.size(); tt++)
				{
					string pdb_atom_name = core_atoms[tt]->get_pdb_atom_name();
					if ( pdb_atom_name == "N" || pdb_atom_name == "CA" || pdb_atom_name == "C") {
						claster_motif_coordinates_[kk][3*tt]=core_atoms[tt]->x();
						claster_motif_coordinates_[kk][3*tt+1]=core_atoms[tt]->y();
						claster_motif_coordinates_[kk][3*tt+2]=core_atoms[tt]->z();
						cout << claster_motif_coordinates_[kk][3*tt] << " ";
						cout << claster_motif_coordinates_[kk][3*tt+1] << " ";
						cout << claster_motif_coordinates_[kk][3*tt+2] << endl;
					}
				}
			}
	}

}



// Из предыдущей версии 

void Cluster_set_DSSP::subtle_claster_show_for_whole_base ( 
 	const string & output_file_name,
	vector < Single_cluster_record > & claster_diversity )
{

//	int number_of_records	= fragment_base_->get_number_of_records();
//	int length_				= fragment_base_->get_lenght (); 

	ofstream output ( output_file_name.c_str());
	if ( ! output )	
	{	
		log_stream << "fill_up_database(): ERROR -  can't create rejected_files.protocol" << endl;
		cout       << "fill_up_database(): ERROR -  can't create rejected_files.protocol" << endl;
		exit (1);	
	}

	int counter =0; 

//	for  ( int ii=claster_diversity_.size()-1 ; ii >= 0; ii-- )


////////////	double * cur_fragment = new double [number_of_records_] ;

	double *cur_fragment = new double [length_*9];

	for  ( int ii=claster_diversity.size()-1 ; ii >= 0; ii-- )
	{
		counter++ ;

		int test = claster_diversity[ii].in_data_base_index();

		fragment_base_->get_coord ( 
			claster_diversity[ii].in_data_base_index() , 
			cur_fragment );
//		Single_fragment * cur_fragment = get_fragment ( claster_diversity_[ii].in_data_base_index() )  ;
		vector < Pair_int_double >  pair_neighbour_index_distance = claster_diversity[ii].pair_neighbour_index_distance() ;
		
		output << endl;

		for ( int jj=0; jj<pair_neighbour_index_distance.size(); jj++ ) 
		{
			output << "CLASTER ";
			PutVa ( counter ,			output,6,	5 ,'l');	
			PutVaDouble ( pair_neighbour_index_distance[jj].value(),output,10,5 ,'r');	

			output << "  " ;

			int test1 = pair_neighbour_index_distance[jj].index();

			fragment_base_->get_coord ( 
				pair_neighbour_index_distance[jj].index() , 
				cur_fragment );

			string	 pdb_ID_chain_ID;
			string	 sequence;
			int		 serial_number;
			string   pdb_resudue_number;

			fragment_base_->fill_up_record_items ( 
				pair_neighbour_index_distance[jj].index() ,
				pdb_ID_chain_ID,
				sequence,
				serial_number,
				pdb_resudue_number);

			PutVa (  pair_neighbour_index_distance[jj].index() ,			output,10,	5 ,'l');	
			PutVa ( pdb_ID_chain_ID	,			output,6,	5 ,'l');	
			PutVa (sequence,					output,10,	1 ,'l');	

			PutVa (serial_number,					output,10,	1 ,'l');	
			PutVa (pdb_resudue_number,				output,10,	1 ,'l');			

			output << " | " ;
		
			vector < double > torsion;
			fill_up_fragment_torsion_angles ( 
				cur_fragment,		
				length_,
				torsion, 
				'd');
		
			for ( int kk =0; kk<torsion.size(); kk++ ) 
				PutVaDouble ( torsion[kk],output,10,3 ,'r');	

			output << endl;
		}
	}
	delete [] cur_fragment;
}
