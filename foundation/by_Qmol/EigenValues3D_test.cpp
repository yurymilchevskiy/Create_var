#pragma warning( disable : 4786 )

#include "EigenValues3D_test.h"
#include "EigenValues3D.h"
#include <iostream>

#include "../Geometry_util/Geometry_util.h" 
#include "../Geometry_util/Macros.h"

#include "LinearAlgebra3D.h"

//void print_matrix ()

using namespace std;
using namespace Geometry_util;

extern ofstream log_stream;

EigenValues3D_test::
~EigenValues3D_test()
{
	cout << "EigenValues3D_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}

void EigenValues3D_test::
eigenvalues3d_test ()
{
	double a[3][3] =	{	{1,1,0},
							{1,1,2},
							{0,2,1} };
	double val[3];

	EigenValues3D( a, val );

}

void EigenValues3D_test::
EigenSystem3D_t ()
{

	double a[3][3] =	{	{1,1,0},
							{1,1,2},
							{0,2,1} };
	double val[3];

	double vect[3][3];

	bool true_flag = EigenSystem3D( a, vect, val );


	double r[9];
	double aa[9],vv[9];

	aa[0] = a[0][0];
	aa[1] = a[0][1];
	aa[2] = a[0][2];
	aa[3] = a[1][0];
	aa[4] = a[1][1];
	aa[5] = a[1][2];
	aa[6] = a[2][0];
	aa[7] = a[2][1];
	aa[8] = a[2][2];

	vv[0] = vect[0][0];
	vv[1] = vect[0][1];
	vv[2] = vect[0][2];
	vv[3] = vect[1][0];
	vv[4] = vect[1][1];
	vv[5] = vect[1][2];
	vv[6] = vect[2][0];
	vv[7] = vect[2][1];
	vv[8] = vect[2][2];



	multiplication_matrixes_3x3_by_3x3	(
		aa,vv, r ) ;

	double vvt[9],f[9];

   	vvt[0] =	vv[0];
	vvt[1] =	vv[3];
	vvt[2] =	vv[6];
	vvt[3] =	vv[1];
	vvt[4] =	vv[4];
	vvt[5] =	vv[7];
	vvt[6] =	vv[2];
	vvt[7] =	vv[5];
	vvt[8] =	vv[8];

	multiplication_matrixes_3x3_by_3x3	(
		vvt, r,f ) ;

}

