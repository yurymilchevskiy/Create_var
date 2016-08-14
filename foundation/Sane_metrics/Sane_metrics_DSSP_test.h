#ifndef SANE_METRICS_DSSP_TEST_H
#define SANE_METRICS_DSSP_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif
 
#include <string>

class  Sane_metrics_DSSP_test: public Simple_test
{
public:
	~Sane_metrics_DSSP_test();

    void run()
    {
	//	show_values_for_dssp_words ();
		//first_constructor_test ();
	fill_up_subset_by_whole_base_test ();

		//create_all_existing_DSSP_words_metrics ();
	}
	void show_values_for_dssp_words ();
	void fill_up_subset_by_whole_base_test ();
	void first_constructor_test ();
	void create_all_existing_DSSP_words_metrics ();
};

#endif
