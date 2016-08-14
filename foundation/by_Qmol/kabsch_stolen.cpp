//#include "stdafx.h"
//#include "kabsch.h"
#include "LinearAlgebra3D.h"
//#include "Excepts.h"
#include <complex>
#include <iostream>

#include "EigenValues3D.h"


using namespace std;

//int EigenSystem3D(double a[3][3], double vect[3][3], double val[3]);
//void EigenValues3D(double a[3][3], double val[3]);

#define		ALIGNMENT_EPS	1.0e-10

// ref ... The reference vector
// vec ... The vector to be aligned to ref
// w   ... A vector of weights
//void kabsch(array<atom> &ref, array<float> &vec, array<float> &w)

void kabsch_stolen(double *ref, double *vec, double *w, const int num, bool & error_flag )
{
	// The number of atoms (i.e. the number of 3 space vectors)
//	int num = ref.Length();
	int i, j, k;

	int small = X;
	int med = X;
	int large = X;

	double R[3][3], RtR[3][3], U[3][3];
	double A[3][3], B[3][3];
	double eigenVal[3];
	double mag;

	for(i = X;i <= Z;i++)
		for(j = X;j <= Z;j++)
			R[i][j] = RtR[i][j] = 0.0;

	double coord[3];

	for(i = 0;i < num;i++){
		memcpy (coord,ref+3*i,3*sizeof (double));

		R[X][X] += w[i]*coord[X]*vec[3*i + X];
		R[X][Y] += w[i]*coord[X]*vec[3*i + Y];
		R[X][Z] += w[i]*coord[X]*vec[3*i + Z];

		R[Y][X] += w[i]*coord[Y]*vec[3*i + X];
		R[Y][Y] += w[i]*coord[Y]*vec[3*i + Y];
		R[Y][Z] += w[i]*coord[Y]*vec[3*i + Z];
 
		R[Z][X] += w[i]*coord[Z]*vec[3*i + X];
		R[Z][Y] += w[i]*coord[Z]*vec[3*i + Y];
		R[Z][Z] += w[i]*coord[Z]*vec[3*i + Z];
	}

	// RtR holds R^{transpose).R
	for(i = X;i <= Z;i++)
		for(j = X;j <= Z;j++)
			for(k = X;k <= Z;k++)
				RtR[i][j] += R[k][i]*R[k][j];

	// Compute the eigenvalues and eigen vectors of RtR.
	//RtR.Jacobi(A, eigenVal);
	if(EigenSystem3D(RtR, A, eigenVal) == false)
	{
		cout << "kabsch:: Unable to diagonalize matrix" << endl;
//		throw QMOL_ERROR("kabsch:: Unable to diagonalize matrix");
		error_flag = true;
		return;
	}

	// Test all of the eigen values to make sure that we have a valid
	// solution.Only skip structures where all of the atoms are colinear
	// (i.e. 2 zero eigenvalues).
	int nullsubspace = 0;

	for(i = X;i <= Z;i++){
		if(eigenVal[i] < ALIGNMENT_EPS){
			//throw QMOL_ERROR("Zero eigen value in kabsch: Colinear and "
			//	"coplanar molecules are not aligned.");
			nullsubspace ++;
/////			error_flag = true;
/////			return;
		}
	}

	if(nullsubspace > 1){
		//throw QMOL_ERROR("Zero eigen values in kabsch: Colinear molecules "	"are not currently allowed.");
		cout << "Zero eigen values in kabsch: Colinear molecules " <<	"are not currently allowed." << endl;
		error_flag = true;
		return;
	}

	// Sort the eigen values (in case we need to remove improper
	// an rotation).
	for(i = Y;i <= Z;i++){
		if(fabs(eigenVal[i]) < fabs(eigenVal[small]))
			small = i;

		if(fabs(eigenVal[i]) > fabs(eigenVal[large]))
			large = i;
	}

	if((large + 1)%3 == small)
		med = (large + 2)%3;
	else
		med = (large + 1)%3;


	double smallVec[3];
	double medVec[3];
	double largeVec[3];

	for(i = X;i <= Z;i++){
		medVec[i] = A[i][med];
		largeVec[i] = A[i][large];
	}

	// *************************************************************************
	// If you get this error message, you need to get the most recent
	// patch for the MS compiler!
	// fatal error C1001: INTERNAL COMPILER ERROR
    // (compiler file 'E:\8168\vc98\p2\src\P2\main.c', line 494)
    // Please choose the Technical Support command on the Visual C++
    // Help menu, or open the Technical Support help file for more information
	// *************************************************************************

	// The cross product defines the eigen vector cooresponding to
	// the smallest eigen value.
	CROSS(smallVec, largeVec, medVec);

	for(i = X;i <= Z;i++)
		A[i][small] = smallVec[i];
	
	/* Invert the eigenvalues as called for in kabsch */
	if(nullsubspace){
		eigenVal[large] = 1.0/sqrt(eigenVal[large]);
		eigenVal[med] = 1.0/sqrt(eigenVal[med]);

		if(eigenVal[large]*eigenVal[med] < 0.0)
			eigenVal[med] *= -1.0;
	}
	else{
		for(i = X;i <= Z;i++)
			eigenVal[i] = 1.0/sqrt(eigenVal[i]);

		if(eigenVal[X]*eigenVal[Y]*eigenVal[Z] < 0.0)
			eigenVal[small] *= -1.0;
	}

    for(i = X;i <= Z;i++)
		for(j = X;j <= Z;j++){
			B[i][j] = 0.0;
			for(k = X;k <= Z;k++)
    			B[i][j] += R[j][k]*A[k][i]*eigenVal[i];
		}
	
	for(i = X;i <= Z;i++){
		medVec[i] = B[med][i];
		largeVec[i] = B[large][i];
	}

	mag = 1.0/sqrt(DOT(medVec,medVec));

	medVec[X] *= mag;
	medVec[Y] *= mag;
	medVec[Z] *= mag;

	mag = 1.0/sqrt(DOT(largeVec,largeVec));

	largeVec[X] *= mag;
	largeVec[Y] *= mag;
	largeVec[Z] *= mag;

	// The cross product defines the eigen vector cooresponding to
	// the smallest eigen value.
	//smallVec = largeVec | medVec;
	CROSS(smallVec, largeVec, medVec);

	for(i = X;i <= Z;i++)
		B[small][i] = smallVec[i];

    for(i = X;i <= Z;i++)
		for(j = X;j <= Z;j++){
			U[i][j] = 0.0;
			for(k = X;k <= Z;k++)
				U[i][j] += A[j][k]*B[k][i];
		}

	for(i = 0;i < num;i++){
		double tmp[3];

		tmp[X] = U[X][X]*vec[3*i + X] + U[X][Y]*vec[3*i + Y] + 
			U[X][Z]*vec[3*i + Z];
		tmp[Y] = U[Y][X]*vec[3*i + X] + U[Y][Y]*vec[3*i + Y] + 
			U[Y][Z]*vec[3*i + Z];
		tmp[Z] = U[Z][X]*vec[3*i + X] + U[Z][Y]*vec[3*i + Y] + 
			U[Z][Z]*vec[3*i + Z];

		vec[3*i + X] = (float)tmp[X];
		vec[3*i + Y] = (float)tmp[Y];
		vec[3*i + Z] = (float)tmp[Z];
    }

	error_flag = false;

}

double rmsd(double *ref, double *vec, double *w, const int num)
{
	double _rmsd = 0;
	double total_mass = 0;

	for (int ii=0; ii<num; ii++ )
	{
		double tmp = ref[X] - vec[X];
		double tmpRmsd = tmp*tmp;

		tmp = ref[Y] - vec[Y];
		tmpRmsd += tmp*tmp;

		tmp = ref[Z] - vec[Z];
		tmpRmsd += tmp*tmp;

		_rmsd += w[ii]*tmpRmsd;

		total_mass += w[ii];
	}

	return (sqrt(_rmsd/total_mass));
}