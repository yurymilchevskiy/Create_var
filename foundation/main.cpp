#pragma warning( disable : 4786 )

#include  "Censorship.h"

//#include "DisProt/DisProt_data/DisProt_data_test.h"
//#include "DisProt/DisProt_data/DisProt_data.h"

//#include "PdbSelectHandle/Get_pdbselect_chains_missing_in_current_set_test.h"
//#include "Chain_store/Chain_store_test.h"

//#include "Fragment_base/Fragment_base_subtle_test.h"
//#include "Sane_metrics/Sane_metrics_test.h"

//#include "Cluster_set/Cluster_set_test.h"	

//#include "Frequency_extrapolation/Frequency_extrapolation_test.h"	

//#include "CowardVariables/FrequencyVariables_test.h"

//#include "CowardVariables/CowardVariables_test.h"

//#include "CowardVariables/Distance_to_claster_variables_test.h"
//#include "Main_model/Abu_Maimonides_Rambam_test.h"

//#include "Statistical_utilits/Discriminant_stepwise_test.h"

//#include "Chain_store/DSSP_binary_test.h"

//#include "MutualDistances/MutualDistances_test.h"

//#include "Main_model/Mordehay_test.h"

//#include "Left_helix_analyser/Left_helix_analyser_test.h"
//#include "Left_helix_analyser/Left_helix_analyser_test.h"

//#include "Ramachandranplot/Ramachandranplot_test.h"

//#include "One-to-one-correspondence/One_to_one_correspondence_test.h"

//#include "Cartesian_by_cluster/Cartesian_by_cluster_test.h"

//#include "Alignment/Class_assignment_profile_test.h"
//#include "Alignment/Alignment_test.h"

//#include "Alignment/Cluster_dssp_interdependence.h"
//#include "Alignment/Cluster_dssp_interdependence_test.h"

#include <iostream>

#include "Pair_int_double.h" 

//#include "Chain_store/Fragment_base_DSSP_bridge_test.h"

//#include "Sane_metrics/Sane_metrics_DSSP_test.h"

//#include "Sane_metrics/Sane_metrics_DSSP_test.h"

//#include "Cluster_set/Cluster_set_DSSP_test.h"

//#include "Fragment_base/Chain_binary_test.h"
//#include "Fragment_base/Chain_binary.h"

#include "ELM/ELM_Test.h"


//#include "Alignment/Cluster_dssp_interdependence_PPII_test.h"
//#include "Alignment/Cluster_dssp_interdependence_PPII.h"


using namespace std;

Censorship configuration;
ofstream log_stream;

int main(int argc,char  **argv) 
{
	string  output_file ("log"); 
	log_stream.open (output_file.c_str()  ); // fix log or not log
	configuration.init("D:/Didona/config ") ;


	if (argc == 1)
	{

	ELM_Test elm_test_;
	elm_test_.run();

//	Cluster_set_DSSP_test	cluster_set_dssp_test_;
//	cluster_set_dssp_test_.run();


//	Sane_metrics_DSSP похоже на Sane_metrics (база данных на основании взаминых расстояний фрагментов строит метрики для кластеризации)
//  Но в выборке только структуры, которые соответствуют заданному DSSP слову (например: "GGGHH")

//		Sane_metrics_DSSP_test sane_metrics_DSSP_test_;
//		sane_metrics_DSSP_test_.run();

/// Делает итератор для прохода по подможествам fragment_base, соответствующих DSSP классам 
/// Например, бетта-слоям или изгибам. Собтвенно итератор нужен для последующей кластеризации 
/// 		внутри даннаго класса
//	Fragment_base_DSSP_bridge_test fragment_base_dssp_bridge_test_;
//	fragment_base_dssp_bridge_test_.run();

// Для предсказания вторичнеой структуры по предсказанным расстояниям до кластеров
//	Cluster_dssp_interdependence_test cluster_dssp_interdependence_test_;
//	cluster_dssp_interdependence_test_.run();

// Для предсказания вторичнеой структуры ВКЛЮЧАЯ PPII  по предсказанным расстояниям до кластеров
//	Cluster_dssp_interdependence_PPII_test cluster_dssp_interdependence_ppii_test_;
//	cluster_dssp_interdependence_ppii_test_.run();

// структурное выравнивание
//		Alignment_test alignment_test_;
//		alignment_test_.run();


// Создает базу данных предсказаний номеров классов. Из расстояний до кластеров вытаскивает номера 
// ближайших классов. Нужно для быстрого доступа при сканировании всей выборки.
// Файлы с предсказанными номерами классов в поддиректории <model_name>/Class_assignment_profile
//		Class_assignment_profile_test class_assignment_profile_test_;
//		class_assignment_profile_test_.run();

//		DisProt_data_test  disprot_data_test_;
//		disprot_data_test_.run();

//		Cartesian_by_cluster_test  cartesian_by_cluster_test_;
//		cartesian_by_cluster_test_.run();

//      Облом - продолжение в классе Cartesian_by_cluster_test
//		One_to_one_correspondence_test  one_to_one_correspondence_test_;
//		one_to_one_correspondence_test_.run();

// Построение карты рамачврдрана для фрагментов в отличие от одиночных остатков
//		Ramachandranplot_test ramachandranplot_test_;
//		ramachandranplot_test_.run();

// контроль наличия списка правильных выборки pdbselect 
	//	Get_pdbselect_chains_missing_in_current_set_test get_pdbselect_chains_missing_in_current_set_test_;
//		get_pdbselect_chains_missing_in_current_set_test_.run();

// Создается база данных белковых цепей.  Хранится в D:\Didona\Store\Chain_store\binary 
//		Chain_store_test chain_store_test_;
//		chain_store_test_.run();

// ИЗ базы набора цепей собираетвся база данных фрагментов для 
//		Fragment_base_subtle_test fragment_base_subtle_test_;
//		fragment_base_subtle_test_.run();

//	Sane_metrics база данных на основании взаминых расстояний фрагментов строит метрики для кластеризации
	//	Sane_metrics_test	sane_metrics_test_;
	//	sane_metrics_test_.run();
//
// Определяются номера структу в базе Fragment_base_subtle соответствующие базисным структурам
// !!!!!!!!!!!!!!!!!!!
	//	 Cluster_set_test cluster_set_test_;
	//	 cluster_set_test_.run();

// создаются базы данных по частотам встречаемости фрагментов
//	Frequency_extrapolation_test	frequency_extrapolation_test_;
//	frequency_extrapolation_test_.run();

// класс генерации переменных связанных со статистикой расстояний фрагментов от базисных структур
//	FrequencyVariables_test  frequencyvariables_test_;
//	frequencyvariables_test_.run();

// класс генерации переменных связанных физико-химическими свойствами аминокислот
//	CowardVariables_test	cowardvariables_test_;
//	cowardvariables_test_.run();

// класс для генерации набора зависимых переменныхC
//	Distance_to_claster_variables_test	distance_to_claster_variables_test_;
//	distance_to_claster_variables_test_.run();

// Главная модель lkz предсказания локальной структуры реализована в функциях этого класса
//	Abu_Maimonides_Rambam_test abu_maimonides_rambam_test_; 
//	abu_maimonides_rambam_test_.run();

// Tuning & testing stepwise discriminant analysis base properties
//	Discriminant_stepwise_test	discriminant_stepwise_test_;
//	discriminant_stepwise_test_.run(); 
//	discriminant_stepwise_test_.run(); 

	//DSSP_binary_test	dssp_binary_test_;
	//dssp_binary_test_.run();

//	MutualDistances_test	mutualdistances_test_;
//	mutualdistances_test_.run();

//	Mordehay_test   mordehay_test_;
//	mordehay_test_.run();

//	Left_helix_analyser_test left_helix_analyser_test_;
//	left_helix_analyser_test_.run();

//	Chain_binary_test chain_binary_test_;
//	chain_binary_test_.run();
	}
	else
	{
//		Command_line_option cmd_opt (argc, argv );
//		working_cad ( cmd_opt );
	}
	return 0;
}