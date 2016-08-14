#ifndef SANE_METRICS_DSSP_H
#define SANE_METRICS_DSSP_H

enum Sane_metrics_DSSP_run_mode
{
	SANE_METRICS_DSSP_FILL_UP,
	SANE_METRICS_DSSP_FILL_UP_SUBSET_BY_WHOLE_BASE,
	SANE_METRICS_COMMON_DSSP_USAGE
};

#include < string >
#include < vector >

using namespace std;

class Fragment_base_subtle;
class Sheduler;
class Zip_iteator;

class Sane_metrics_DSSP
{
public:

	Sane_metrics_DSSP   (  
				 string & fragment_base_subtle_name,
				 string & dssp_word,
		         Sane_metrics_DSSP_run_mode run_mode) ;

	~Sane_metrics_DSSP  ();

	void fill_metrics ( const int power_index, double *metrics );

	int get_number_of_records () const { return number_of_records_; }

	string get_host_dir () const { return host_dir_; }

	string get_fragment_base_subtle_name () const {
		return fragment_base_subtle_name_;} 

	string get_sane_metrics_name () const {		return sane_metrics_name_;} 

	int get_zip_factor() const { 
		return zip_factor_; 
	}

	vector <int> get_base_indexes_for_dssp_words () const {
		return base_indexes_for_dssp_words_;
	}
		

	int get_shift()		 const { return shift_; }
	int get_length()		 const { return length_; }

	Zip_iteator	* get_Zip_iteator()  {	 return zip_it_; }

	void 	calibration_single_extenal_structure (    // ПЕРЕЕЗЖДАЕТ В КЛАСС Ramachandrabpkop
		double *probe_fragment,
		vector<int> & neighbour_distribution,
		double & metric_1,
		double & metric_2,
		double & metric_3,
		double & metric_4,
		double step)	;		// будем видимо 0.1 использовать 

private:
	void fill_up();
	void fill_up_subset_by_whole_base ();
	void init ();

	string fragment_base_subtle_name_;  // имя базы Fragment_base_subtle. В этой базе уже есть данные о длинах фрагментов и проч.


	string sane_metrics_name_;

	string dssp_word_;

	void allocate ();
	
	void		prepare_host_dir();  
	string		bin_path_;
	string		host_dir_;

	Sheduler				*sheduler_;	

	Fragment_base_subtle	*fragment_base_;

	Zip_iteator				*zip_it_;
	int						zip_factor_;
	int						shift_;

	int			number_of_records_;   // число обабатываемых записей оставщихся после редукции Fragment_base_subtle в zip_factor_ раз

	double a_val_;

	int length_;

	double  *metrics_sqrt_;
	double  *metrics_1_;
	double  *metrics_2_;
	double  *metrics_3_;
	double  *metrics_4_;

	double **all_fragments_; 
	
	int total_fragment_number_ ; //= fragment_base_->get_total_fragment_number(); SIZE of fragment base
	double **all_BASE_fragments_; // all fragments from base


	string bridge_name_;

	vector <int > base_indexes_for_dssp_words_;   // index set in fragment base corresponding to chosen dssp_word

	int maх_sample_size_ ; // max allowed size for fragment base subset used to form metrics
};

#endif	