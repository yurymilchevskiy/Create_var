#ifndef SANE_METRICS_TEST_H
#define SANE_METRICS_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  Sane_metrics_test: public Simple_test
{
public:
	~Sane_metrics_test ();

    void run()
    {
	//	test_subbtle_error();

//		test_of_zip_integrity();   // еще и создает что положено

// очень полезно но лучше перенести к базисным структурам 
//		distance_to_canonical_structure_from_BS();


	//	calibration_single_extenal_structure();

		multiple_calibrating_test ();   // готовит карту рамачандрана
	//	prepare_rama_plot ();


	}
	void prepare_rama_plot ();
	void test_subbtle_error();
	void test_of_zip_integrity();
	void test1 ();
	void check_rmsd_calc(); 
	void distance_to_canonical_structure_from_BS();
	void calibration_single_extenal_structure();
	void multiple_calibrating_test ();
};

#endif
