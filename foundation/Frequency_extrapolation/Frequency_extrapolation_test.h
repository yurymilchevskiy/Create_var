#ifndef FREQUENCY_EXTRAPOLATION_TEST_H
#define FREQUENCY_EXTRAPOLATION_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  Frequency_extrapolation_test: public Simple_test
{
public:
	~Frequency_extrapolation_test();

    void run()
    {

		//	fill_up_setted_basis_structure();	// CREATES DATABASEAS
		//	fill_up();			// CREATES DATABASEAS
		//	first_test_extended();			// CREATES DATABASEAS



//**********************************************************
		//analyse_cluster_only_test();

	//	first_test();			// CREATES DATABASEAS


//       check_integrity ();
//		second_test();
/*		show_test();
		get_single_distance_record_test();
		prepare_jack_nife_chain_data_test ();
		*/
//	show_freq_data_test ();
		prepare_together_freq_data_test ();


	}
	void fill_up_setted_basis_structure();
	void fill_up();
	void first_test_extended();
	void check_integrity ();
	void first_test();
	void second_test();
	void show_test();
	void get_single_distance_record_test();
	void prepare_jack_nife_chain_data_test ();
	void show_freq_data_test ();
	void prepare_together_freq_data_test ();
	void analyse_cluster_only_test();

};

#endif
