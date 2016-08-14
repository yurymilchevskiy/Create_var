#pragma warning( disable : 4786 )

#include "Sane_metrics.h"

#include "../Fragment_base/Fragment_base_subtle.h"

#include "../Zip_iteator.h"

#include <direct.h>

//#include "Dunghill.h"

#include "../Censorship.h"
#include "../CommonFunc.h"
#include "../Sheduler.h"

#include "../by_Qmol/kabsch_stolen.h"

extern Censorship configuration;
extern ofstream log_stream;

//#include "nimble_distance_calc.h"

Sane_metrics::~Sane_metrics()
{
	if ( sheduler_ )				
		delete sheduler_;

	if ( metrics_sqrt_	) 
		delete [] metrics_sqrt_ ; 
    if ( metrics_1_		)
		delete [] metrics_1_; 
    if ( metrics_2_		)
		delete [] metrics_2_; 
    if ( metrics_3_		)
		delete [] metrics_3_; 
    if ( metrics_4_	    )
		delete [] metrics_4_; 

	sheduler_ = 0;

	metrics_sqrt_= 0;
	metrics_1_= 0;
	metrics_2_= 0;	
	metrics_3_= 0;
	metrics_4_= 0;


	if ( all_fragments_ != 0 )  
	{
		for ( int ii=0; ii < number_of_records_; ii++  )
			delete [] all_fragments_[ii];
		delete [] all_fragments_;
	}

	all_fragments_ = 0;

	if ( fragment_base_ != 0) 
		delete fragment_base_ ;

	if ( zip_it_ != 0) 
		delete zip_it_ ;
}

Sane_metrics::
Sane_metrics   (  
				// string & fragment_base_subtle_name,
				string & sane_metrics_name,
//				const int zip_factor,
//				const int shift,
				Sane_metrics_run_mode run_mode) :
	sane_metrics_name_ 			(sane_metrics_name),
//	zip_factor_		(zip_factor),
//	shift_			(shift),
	zip_factor_		(0),
	shift_			(0),
	fragment_base_  (0),
    sheduler_		(0),
	zip_it_			(0),
    metrics_sqrt_	(0),	
    metrics_1_		(0),		
    metrics_2_		(0),		
    metrics_3_		(0),		
    metrics_4_		(0),
	all_fragments_  (0)
	

{
	//fragment_base_subtle_name_ 		= 	fragment_base_subtle_name;

	sheduler_					= new Sheduler  (	
		configuration.option_meaning("Path_to_Sane_metrics")  +  
		sane_metrics_name_  + string ("/") + 
		string ("sheduler") ) ;


	fragment_base_subtle_name_ =		
		sheduler_->option_meaning ("FRAGMENT_BASE_SUBTLE_NAME").c_str() ;


	fragment_base_  = new Fragment_base_subtle ( fragment_base_subtle_name_ , FRAGMENT_BASE_SUBTLE_COMMON_USAGE);
	int size_of_fragment_base =  fragment_base_->get_total_fragment_number () ;  


//	zip_factor_ = atoi( sheduler_->option_meaning("ZIP_FACTOR").c_str() );
	zip_factor_ = 1 ; // FIX 
	shift_		= atoi( sheduler_->option_meaning("SHIFT").c_str() );
	a_val_		= atof( sheduler_->option_meaning("DEFAULT_DENOMINATOR").c_str() );

	zip_it_ = new Zip_iteator (
					size_of_fragment_base,
					zip_factor_,
					shift_);	

	number_of_records_ = zip_it_->get_step_number();
	
	length_  =  fragment_base_->get_fragment_length();

//	a_val_ = 0.10;   // FIX - надо бы покультурней 

	//prepare_host_dir();
//	bin_path_ = get_host_dir () +  string ("/") + 	
//		string ("sm.bin")  ;

	bin_path_ = configuration.option_meaning("Path_to_Sane_metrics")  +  
		sane_metrics_name_ + string ("/") +string ("sm.bin")  ;


	if		( run_mode == SANE_METRICS_FILL_UP )			fill_up ();
	else if ( run_mode == SANE_METRICS_COMMON_USAGE)		init ();
	else 
	{
		log_stream	<< "Fragment base " <<  bin_path_ << "error"  << endl;
		cout		 << "Fragment base " <<  bin_path_<< "error"  << endl;
		exit (1);
		
	}
}

void Sane_metrics::
allocate ()
{
		metrics_sqrt_  = new double [ number_of_records_*sizeof (double) ];
		metrics_1_     = new double [ number_of_records_*sizeof (double) ];
		metrics_2_		= new double [ number_of_records_*sizeof (double) ];	
		metrics_3_		= new double [ number_of_records_*sizeof (double) ];	
		metrics_4_		= new double [ number_of_records_*sizeof (double) ];	

 		memset ( metrics_sqrt_,				0,number_of_records_*sizeof (double) ) ;
		memset ( metrics_1_,				0,number_of_records_*sizeof (double) ) ;
		memset ( metrics_2_,				0,number_of_records_*sizeof (double) ) ;
		memset ( metrics_3_,				0,number_of_records_*sizeof (double) ) ;
		memset ( metrics_4_,				0,number_of_records_*sizeof (double) ) ;

		all_fragments_ = new double* [number_of_records_] ;	memset (all_fragments_,0,sizeof(double*)	 *number_of_records_);
		for ( int ii=0; ii < number_of_records_; ii++  )
			all_fragments_[ii]  = new double [length_*9];

}

void Sane_metrics::
prepare_host_dir()
{
	ostringstream ost ;
	ost << zip_factor_ << "_" << shift_;
	string local_name = ost.str();
	
	host_dir_ = 
		configuration.option_meaning("Path_to_Sane_metrics")  +  
		sane_metrics_name_ + string ("/") + local_name   ;

//	mkdir (host_dir_.c_str());
}


void Sane_metrics::
init ()
{
}

void Sane_metrics::
fill_up()
{


	_mkdir (host_dir_.c_str()); 

	ofstream 	out_stream(  bin_path_.c_str(),ios::binary );
    if ( ! out_stream )   	{	
		log_stream		<< " can't create " << bin_path_ << endl;		
		cout			<< " can't create " << bin_path_ << endl;		
		exit (1);		}

	allocate ();

	double *coord_array_ii = new double [length_*9];
	double *coord_array_jj = new double [length_*9];

//	for ( int kk=0; kk< number_of_records_; kk++ )				fragment_base_->get_coord ( kk, all_fragments_[kk] );


	int number_of_choosen_record = 0;
	int kkkk=0;
    while (1)
    {
		int curr_rec_number_in_fragment_base = zip_it_->current();
		fragment_base_->get_coord ( 
			curr_rec_number_in_fragment_base , 
			all_fragments_[number_of_choosen_record] );

		number_of_choosen_record++;

		if ( zip_it_->has_next() )	zip_it_->next();
		else 					break;

		kkkk++;
	} 

/*
	int number_of_choosen_record = 0;
	for ( int kk=shift_learning_position_; kk< number_of_records_; kk+=reduce_learning_base_size_  )
	{
		fragment_base_->get_coord ( kk, all_fragments_[number_of_choosen_record ] );
		number_of_choosen_record ++;
	}
*/

	double *w	= new double [(length_*9)];
	for ( int kk=0;kk<(length_*9);kk++)		w[kk] = 1;

	for ( int ii=0; ii< number_of_choosen_record ; ii++)
	{

//		fragment_base_->get_coord ( ii, coord_array_ii);


		for ( int  jj=0; jj< number_of_choosen_record; jj++)
		{

//			fragment_base_->get_coord ( jj, coord_array_jj);

			bool error_flag= true;

			double current_distance = kabsch_rmsd (all_fragments_[ii], all_fragments_[jj], w, (3*length_), error_flag);

			if (error_flag)
			{
				cout << ii << "\t"  << jj << "\t Error kabsch_rmsd " << endl;
				log_stream << ii << "\t   "  << jj << "\t Error kabsch_rmsd " << endl;
				continue;
			}

			metrics_sqrt_	[ii] += 1.0 /( a_val_ + sqrt (current_distance) ) ;
			metrics_1_		[ii] += 1.0 /( a_val_ + current_distance  ) ;
			metrics_2_		[ii] += 1.0 /( a_val_ + current_distance * current_distance ) ;
			metrics_3_		[ii] += 1.0 /( a_val_ + current_distance * current_distance *current_distance  ) ;
			metrics_4_		[ii] += 1.0 /( a_val_ + current_distance * current_distance * current_distance * current_distance  ) ;
//			cout << jj << endl;
//			if ( jj  == 37308 )				cin >> ddd;
		}
	//	if ( ii%100 == 0 )
			cout << ii << endl;
	}

	delete [] w;

	out_stream.write ( (char* ) metrics_sqrt_,number_of_records_*sizeof (double)  );
	out_stream.write ( (char* ) metrics_1_	 ,number_of_records_*sizeof (double)  );
	out_stream.write ( (char* ) metrics_2_	 ,number_of_records_*sizeof (double)  );
	out_stream.write ( (char* ) metrics_3_	 ,number_of_records_*sizeof (double)  );
	out_stream.write ( (char* ) metrics_4_	 ,number_of_records_*sizeof (double)  );

	delete [] coord_array_ii;
	delete [] coord_array_jj;

}

void Sane_metrics::
fill_metrics ( 
	const int power_index, 
	double *metrics )
{
	string bin_file_name = bin_path_;		

	ifstream 	in_stream(  bin_file_name.c_str(),ios::binary );
    if ( ! in_stream )   	{	
		log_stream		<< " can't create metrics file " << endl;		
		cout			<< " can't create metrix file " << endl;		
		exit (1);		}
	
	int position = power_index * number_of_records_*sizeof (double);
	in_stream.seekg(position,ios::beg);

	in_stream.read ( (char* ) metrics 	 ,number_of_records_*sizeof (double)  );

}

void Sane_metrics::
	calibration_single_extenal_structure ( 
		double *probe_fragment,
		vector<int> & neighbour_distribution,
		double & metric_1,
		double & metric_2,
		double & metric_3,
		double & metric_4,
		double step)			// будем видимо 0.1 использовать 
{
	int total_fragment_number = fragment_base_->get_total_fragment_number () ;

	double *current_fragment  = new double [length_*9];

	double *probe_fragment_copy  = new double [length_*9];

	double *w	= new double [(length_*9)];
	for ( int kk=0;kk<(length_*9);kk++)		w[kk] = 1;

	for ( int ii=0; ii< total_fragment_number; ii+=100)
	{

			memcpy(probe_fragment_copy,probe_fragment,length_*9*sizeof(double)    );

			fragment_base_->get_coord ( 
				ii, 
				current_fragment );

			bool error_flag= true;
			double current_distance = kabsch_rmsd (current_fragment, probe_fragment_copy, w, (3*length_), error_flag);

			if (error_flag)
			{
				cout << ii <<  "\t Error kabsch_rmsd " << endl;
				log_stream << ii <<  "\t Error kabsch_rmsd " << endl;
				continue;
			}

//			metrics_sqrt_	[ii] += 1.0 /( a_val_ + sqrt (current_distance) ) ;
			metric_1		 += 1.0 /( a_val_ + current_distance  ) ;
			metric_2		 += 1.0 /( a_val_ + current_distance * current_distance ) ;
			metric_3		 += 1.0 /( a_val_ + current_distance * current_distance *current_distance  ) ;
			metric_4		 += 1.0 /( a_val_ + current_distance * current_distance * current_distance * current_distance  ) ;

			double temp  = current_distance/step;
			int index = (int)  (current_distance/step);

			index = ( index <= 99 ) ? index : 99;
			neighbour_distribution[index] ++;
	}
	delete [] current_fragment;
	delete [] probe_fragment_copy;
	delete [] w;
}
