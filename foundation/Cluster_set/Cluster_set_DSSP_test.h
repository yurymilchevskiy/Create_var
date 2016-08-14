#ifndef CLUSTER_SET_DSSP_TEST_H
#define CLUSTER_SET_DSSP_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  Cluster_set_DSSP_test: public Simple_test
{
public:
	~Cluster_set_DSSP_test();

    void run()
    {


//  перетряхивает оч. большой файл subtle_clasters_show, который СОЗДАН КЛАССОМ Cluster_set.
// добавляет поля DSSP для каждого фрагмента. 
// Это чтоб смотрететь как убывает связь текщего кдастера и DSSP word по мере увеличения RMSD от базисной структуры
	peretrah_subtle_clasters_show ();



	//	make_claster_origin_list_and_names_test();
		// List_files_test (); не вышло ни хера
//		void predict_dssp_eight_letter_by_sequence_test ();

	//	optimize_clasterization_test ();

	//	handle_all_dssp_words ();

		//prepare_indexes_and_names_for_final_clasterization_test ();
	}
	void peretrah_subtle_clasters_show ();
	void make_claster_origin_list_and_names_test();
	void List_files_test ();
	void prepare_indexes_and_names_for_final_clasterization_test ();
	void handle_all_dssp_words ();
	void predict_dssp_eight_letter_by_sequence_test ();
	void optimize_clasterization_test();
	void newlife_constructor_test();
	void constructor_test ();
//	void optimize_clasterization_test ();
	void pull_out_claster_origin_structure_list_test ( );
	void mutual_distance_for_BS_show_test();
	void analyse_analyse_setted_regular_structure_presence();
};

#endif
