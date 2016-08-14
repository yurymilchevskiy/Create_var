#ifndef SANE_METRICS_H
#define SANE_METRICS_H

enum Sane_metrics_run_mode
{
	SANE_METRICS_FILL_UP,
	SANE_METRICS_COMMON_USAGE
};


#include < string >
#include < vector >

using namespace std;

//class Dunghill;
class Fragment_base_subtle;
class Sheduler;
class Zip_iteator;

class Sane_metrics
{
public:

//	Sane_metrics   (  Dunghill *fragment_base ); 
//	Sane_metrics   (  const string & name ); 

	Sane_metrics   (  
				 string & fragment_base_subtle_name,
//				const int zip_factor,
//				const int shift,
				Sane_metrics_run_mode run_mode) ;


	~Sane_metrics  ();

	void fill_metrics ( const int power_index, double *metrics );
	// power_index: 0 - sqrt,  the rest as is.

	int get_number_of_records () const { return number_of_records_; }

	string get_host_dir () const { return host_dir_; }

	string get_fragment_base_subtle_name () const {
		return fragment_base_subtle_name_;} 



		string get_sane_metrics_name () const {		return sane_metrics_name_;} 


	int get_zip_factor() const { 
		return zip_factor_; 
	}
	int get_shift()		 const { return shift_; }

	int get_length()		 const { return length_; }


/*	void Sane_metrics::	calibration_single_extenal_structure ( 
		double *probe_fragment,
		vector<int> neighbour_distribution,
		double & metric_1,
		double & metric_2,
		double & metric_3,
		double & metric_4,
		double step);			// будем видимо 0.1 использовать 
*/


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
	void init ();

	string fragment_base_subtle_name_;  // имя базы Fragment_base_subtle. В этой базе уже есть данные о длинах фрагментов и проч.


	string sane_metrics_name_;

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


//	int shift_learning_position_;
//	int reduce_learning_base_size_;

};

#endif	