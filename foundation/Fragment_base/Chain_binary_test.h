#ifndef CHAIN_BINARY_TEST1_H
#define CHAIN_BINARY_TEST1_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  Chain_binary_test: public Simple_test
{
public:
	~Chain_binary_test();

    void run()
    {
	//	remind_fragmnet_test ();

//		two_chain_distance_set_test();
//		test1 ();
//		get_fragment_test1 ();
//		EigenValues3D_test ();
//		fragment_to_principal_axes_test ();
//		kabsch_stolen_test();
//		kabsch_rmsd_test();

		center_of_mass_test();
	}
	void remind_fragmnet_test ();
	void two_chain_distance_set_test();
	void test1 ();
	void get_fragment_test1 ();
	void EigenValues3D_test ();
	void fragment_to_principal_axes_test ();
	void kabsch_stolen_test();
	void kabsch_rmsd_test();
	void center_of_mass_test();
};

#endif
