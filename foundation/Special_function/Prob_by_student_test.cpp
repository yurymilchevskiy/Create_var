#pragma warning( disable : 4786 )

#include "Prob_by_student_test.h"

#include "../Censorship.h"


#include "special_function.h"

#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

extern Censorship configuration;
extern ofstream log_stream;

Prob_by_student_test::~Prob_by_student_test()
{
	cout << "Prob_by_student_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}


void Prob_by_student_test::first_test ()
{
	int njue = 20;

	for ( double t= -3; t<3; t+=0.1 )
	{
		double current_prob = prob_by_student (t, njue);
		cout << t << "\t "  << current_prob << endl;
	}

cout << "*************************" << endl;

//	int njue = 10;
//	double t_arr [8] = {0.260, 0.700,1.372,1.812,2.228, 2.764,3.169,4.587};

	double t_arr [8] = {0.257, 0.687,1.325,1.725,2.086, 2.528,2.845,3.850};


	for (int ii=0;ii<8;ii++)
	{
		double current_prob  = prob_by_student (t_arr [ii], njue);
		cout << t_arr [ii] << "\t "  << current_prob << endl;
	}



}

/*
0.26     0.199861			0.257    0.200196
0.7      0.500112			0.687    0.500028
1.372    0.799945			1.325    0.799889
1.812    0.899925			1.725    0.900052
2.228    0.949988			2.086    0.950004
2.764    0.980008			2.528    0.980001
3.169    0.989995			2.845    0.989992
4.587    0.999				3.85     0.999001
*/

