#ifndef CLUSTER_SET_TEST_H
#define CLUSTER_SET_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  Cluster_set_test: public Simple_test
{
public:
	~Cluster_set_test();

    void run()
    {
		pull_out_claster_origin_structure_list_test ( );  // Похоже правильный последний вариант
//		predict_dssp_eight_letter_by_sequence_test ();

	//	prepare_cluster_picture_and_description ();
	//	prepare_cluster_distance_matrix_test();

//		check_manually_settin_mode();
//		optimize_clasterization_test();   // ИМЕННО ТУТ ГЕНЕРЯТСЯ BS!!!!!

//		check_COMMON_USAGE_CLUSTER_SET_MODE_test ();

//		newlife_constructor_test();
//		constructor_test ();
//		optimize_clasterization_test ();

	//	mutual_distance_for_BS_show_test();
//		analyse_analyse_setted_regular_structure_presence();

	}

	void predict_dssp_eight_letter_by_sequence_test ();
	void prepare_cluster_picture_and_description ();
	void prepare_cluster_distance_matrix_test();
//	void check_get_length_etc();
	void check_manually_settin_mode();
	void check_COMMON_USAGE_CLUSTER_SET_MODE_test ();
	void optimize_clasterization_test();
	void newlife_constructor_test();
	void constructor_test ();
//	void optimize_clasterization_test ();
	void pull_out_claster_origin_structure_list_test ( );
	void mutual_distance_for_BS_show_test();
	void analyse_analyse_setted_regular_structure_presence();
};

#endif
