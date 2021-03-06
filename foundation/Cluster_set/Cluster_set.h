#ifndef CLUSTER_SET_H
#define CLUSTER_SET_H

//#include "Pair_int_double.h"
#include "�luster_set_operating_mode.h"

#include <vector>
#include <string>
#include <fstream>
#include <map>

#include "Single_cluster_record.h"

using namespace std;

class	Sheduler;

class	Zip_iteator;

class	Fragment_base_subtle;
class	Sane_metrics;

class  Cluster_set
{
public:

   ~Cluster_set();
   Cluster_set () {};


	Cluster_set::Cluster_set (
		const string & cluster_set_name_name,
		�luster_set_operating_mode run_mode);

  	
	vector < vector <double> >  prepare_cluster_distance_matrix() ;

	
	vector <int>  algorithm_2 (
		double upper_dist,				
		double lower_dist,				
		double step_dist,				
		int waited_claster_number,
		string & metrics_mode);


//	void	plain_claster_show ( const string & output_file_name	);



	void regulate_choseen_index_subsample_version( // fills properties for subset made by zip_it_ selection
		const vector <int> & choosen_indexes,
		vector < Single_cluster_record > & claster_diversity,
		double & distance_scattering,
		double & quadrate_distance_scattering) ;

	void plain_claster_show_subsample_version ( 
		const string & output_file_name,
		vector < Single_cluster_record > & claster_diversity,
		double & distance_scattering,
		double & quadrate_distance_scattering );


	void 	regulate_choseen_index( 
		const vector <int> & choosen_indexes,
		vector < Single_cluster_record > & claster_diversity,
		double & distance_scattering,
		double & quadrate_distance_scattering) ;


	void regulate_choseen_index_for_whole_base( //  choosen_indexes - relating to total ( not reduced base) !!!
	//	const vector <int> & choosen_indexes,
		vector < Single_cluster_record > & claster_diversity,
		double & distance_scattering,
		double & quadrate_distance_scattering) ;

	void plain_claster_show ( 
		const string & output_file_name,
		vector < Single_cluster_record > & claster_diversity,
		double & distance_scattering,
		double & quadrate_distance_scattering);



	void subtle_claster_show_for_whole_base ( 
 		const string & output_file_name,
		vector < Single_cluster_record > & claster_diversity );

	void	expand_clasterization_hystory_protocol ();

	void expand_clasterization_hystory_protocol (
		vector < Single_cluster_record > & claster_diversity,
		double & distance_scattering,
		double & quadrate_distance_scattering);

	bool raise_choosen_indexes 
		(const int			calculation_index,
		int					&	number_of_clasters,
		vector < int >		&	claster_index_in_base );

	void optimize_clasterization () ;

	void mutual_distance_for_BS_show ();

	void mutual_distance_for_BS_show ( 
		vector <int > & choosen_indexes,
		vector < vector <double> > & Cluster_distance_matrix);

// data for COMMON_USAGE_CLUSTER_SET_MODE runmode;	
	vector <int> get_claster_motif_index()  const {  return claster_motif_index_; }
	int number_of_classes() const {return number_of_classes_;} 
	int fragment_length  () const {return fragment_length_;};
	string get_sane_metrics_mame() const {return sane_metrics_mame_;} 
	string get_cluster_set_mame () const {return cluster_set_mame_;}


	double ** get_claster_motif_coordinates() const { return claster_motif_coordinates_;} ;


//	string get_cluster_set_mame () const { return cluster_set_mame_; }

// 
	Fragment_base_subtle	 * get_fragment_base() const {return fragment_base_;}

	string get_host_dir() const { return host_dir_;}

			void warm_up_metrics( 
				double *metrics,
				const string &metrics_mode);


	
private:

		bool 	assign_best_parameters (
			double *distance_scattering_min,
			double *quadrate_distance_scattering_min,
			const string & protocol_file_name  );


		string cluster_set_mame_;

		string sane_metrics_mame_;

		string host_dir_;					// ���������� ���������� Sane_metrics � ��� �� ��� ��� ������� ���������� ����

		Sheduler		*sheduler_;	
//		Sheduler		*clasterization_option_;						// �� ����� ���� ����� ��������� �������������
		
		Fragment_base_subtle	 * fragment_base_;
		Sane_metrics * sa_me_ ;

		
		Zip_iteator				*zip_it_;
		int						zip_factor_;
		int						shift_;

	//	int		number_of_records_;
		int		subsample_size_;		// sample generated by zip_it_ has this size
		int		length_;

		double **subsample_fragments_;
		



		map < int, int > zip_to_base_index_;


// data for COMMON_USAGE_CLUSTER_SET_MODE runmode;
		double **claster_motif_coordinates_;
		vector <int> claster_motif_index_;
		int number_of_classes_;
		int fragment_length_;


		void init_claster_motif();
};

bool raise_choosen_indexes (
	const string		& name_claster_store,
	const int			  calculation_index,
	int					& number_of_clasters,
	vector < int >		& claster_index_in_base );

vector <int > pull_out_claster_origin_structure_list (
	const string & exact_path_to);

void  make_claster_origin_list_and_names (
	const string & exact_path_to_store,
	const string & name,
	const string & extension,
	map < int, string > & cl_in_base_nu_to_dssp_word );



#endif

