#include "LinearAlgebra3D.h"
#include <cmath>
#include <complex>

using namespace std;

void EigenValues3D(double a[3][3], double val[3]);

// Compute the eigenvectors and eigenvales for a 3X3 symmetric matrix.
// Return true (i.e. 1) if the matrix has been successfully diagonalized
// or false (i.e. 0) if matrix is singular.
bool EigenSystem3D(double a[3][3], double vect[3][3], double val[3])
{
	// A small number
	const double eps = 1.0e-7;

	// The eigen vectors (returned as columns of the matrix)
	if(fabs(a[X][Z]) < eps){ // Don't divide by zero!
		
		if((fabs(a[X][Y]) < eps) || (fabs(a[Y][Z]) < eps)){
			if((fabs(a[X][Y]) < eps) && (fabs(a[Y][Z]) < eps)){
				// Matrix is diagonal!

				// Compute the eigenvalues
				val[X] = a[X][X];
				val[Y] = a[Y][Y];
				val[Z] = a[Z][Z];

				vect[X][X] = 1.0;
				vect[Y][X] = 0.0;
				vect[Z][X] = 0.0;

				vect[X][Y] = 0.0;
				vect[Y][Y] = 1.0;
				vect[Z][Y] = 0.0;

				vect[X][Z] = 0.0;
				vect[Y][Z] = 0.0;
				vect[Z][Z] = 1.0;

				// Success!
				return true;
			}
			else{
				if(fabs(a[X][Y]) < eps){
					// Compute the eigenvalues
					val[X] = a[X][X];
					val[Y] = 0.5*(a[Y][Y] + a[Z][Z] - 
						sqrt(a[Y][Y]*a[Y][Y] + 4.0*a[Y][Z]*a[Y][Z] - 
						2.0*a[Y][Y]*a[Z][Z] + a[Z][Z]*a[Z][Z]));
					val[Z] =  0.5*(a[Y][Y] + a[Z][Z] +
						sqrt(a[Y][Y]*a[Y][Y] + 4.0*a[Y][Z]*a[Y][Z] - 
						2.0*a[Y][Y]*a[Z][Z] + a[Z][Z]*a[Z][Z]));

					vect[X][X] = 1.0;
					vect[Y][X] = 0.0;
					vect[Z][X] = 0.0;

					vect[X][Y] = 0.0;
					vect[Y][Y] = -0.5*(-a[Y][Y] + a[Z][Z] +
						sqrt(a[Y][Y]*a[Y][Y] + 4.0*a[Y][Z]*a[Y][Z] - 
						2.0*a[Y][Y]*a[Z][Z] + a[Z][Z]*a[Z][Z]))/a[Y][Z];
					vect[Z][Y] = 1.0;

					vect[X][Z] = 0.0;
					vect[Y][Z] = -0.5*(-a[Y][Y] + a[Z][Z] -
						sqrt(a[Y][Y]*a[Y][Y] + 4.0*a[Y][Z]*a[Y][Z] - 
						2.0*a[Y][Y]*a[Z][Z] + a[Z][Z]*a[Z][Z]))/a[Y][Z];
					vect[Z][Z] = 1.0;
				}
				else{ //(fabs(a[Y][Z]) < eps)
				
					// Compute the eigenvalues
					val[X] = 0.5*(a[X][X] + a[Y][Y] - sqrt(a[X][X]*a[X][X] + 4.0*a[X][Y]*a[X][Y] - 
						2.0*a[X][X]*a[Y][Y] + a[Y][Y]*a[Y][Y]));
					val[Y] = 0.5*(a[X][X] + a[Y][Y] + sqrt(a[X][X]*a[X][X] + 4.0*a[X][Y]*a[X][Y] - 
						2.0*a[X][X]*a[Y][Y] + a[Y][Y]*a[Y][Y]));
					val[Z] = a[Z][Z];

					vect[X][X] = -0.5*(-a[X][X] + a[Y][Y] + sqrt(a[X][X]*a[X][X] + 4.0*a[X][Y]*a[X][Y] - 
						2.0*a[X][X]*a[Y][Y] + a[Y][Y]*a[Y][Y]))/a[X][Y];
					vect[Y][X] = 1.0;
					vect[Z][X] = 0.0;

					vect[X][Y] = -0.5*(-a[X][X] + a[Y][Y] - sqrt(a[X][X]*a[X][X] + 4.0*a[X][Y]*a[X][Y] - 
						2.0*a[X][X]*a[Y][Y] + a[Y][Y]*a[Y][Y]))/a[X][Y];
					vect[Y][Y] = 1.0;
					vect[Z][Y] = 0.0;

					vect[X][Z] = 0.0;
					vect[Y][Z] = 0.0;
					vect[Z][Z] = 1.0;
				}
			}
		}
		else{
			// Compute eigenvalues
			EigenValues3D(a, val);

			vect[X][X] = -(a[Y][Z])/(a[X][Y]) + ((a[Y][Y]-val[X])*(-val[X]+a[Z][Z]))/(a[X][Y]*a[Y][Z]);
			vect[Y][X] = -(-val[X]+a[Z][Z])/(a[Y][Z]);
			vect[Z][X] = 1.0;

			vect[X][Y] =  -(a[Y][Z])/(a[X][Y]) + ((a[Y][Y]-val[Y])*(-val[Y]+a[Z][Z]))/(a[X][Y]*a[Y][Z]);
			vect[Y][Y] = -(-val[Y]+a[Z][Z])/(a[Y][Z]);
			vect[Z][Y] = 1.0;

			vect[X][Z] =  -(a[Y][Z])/(a[X][Y]) + ((a[Y][Y]-val[Z])*(-val[Z]+a[Z][Z]))/(a[X][Y]*a[Y][Z]);
			vect[Y][Z] = -(-val[Z]+a[Z][Z])/(a[Y][Z]);
			vect[Z][Z] = 1.0;
		}
	}
	else{
		// Compute eigenvalues
		EigenValues3D(a, val);

		vect[X][X] = -(-val[X] + a[Z][Z])/(a[X][Z]) + 
					  (a[Y][Z]*(a[X][Z]*a[Y][Z]-a[X][Y]*(-val[X]+a[Z][Z])))/
					  (a[X][Z]*(-a[X][Y]*a[Y][Z]+a[X][Z]*(a[Y][Y]-val[X])));
		vect[Y][X] = -(a[X][Z]*a[Y][Z]-a[X][Y]*(-val[X]+a[Z][Z]))/
					  (-a[X][Y]*a[Y][Z]+a[X][Z]*(a[Y][Y]-val[X]));
		vect[Z][X] = 1.0;

		vect[X][Y] = -(-val[Y]+a[Z][Z])/(a[X][Z]) + 
					  (a[Y][Z]*(a[X][Z]*a[Y][Z]-a[X][Y]*(-val[Y]+a[Z][Z])))/
					  (a[X][Z]*(-a[X][Y]*a[Y][Z]+a[X][Z]*(a[Y][Y]-val[Y])));
		vect[Y][Y] = -(a[X][Z]*a[Y][Z]-a[X][Y]*(-val[Y]+a[Z][Z]))/
					  (-a[X][Y]*a[Y][Z]+a[X][Z]*(a[Y][Y]-val[Y]));
		vect[Z][Y] = 1.0;

		vect[X][Z] = -(-val[Z]+a[Z][Z])/(a[X][Z]) + 
					  (a[Y][Z]*(a[X][Z]*a[Y][Z]-a[X][Y]*(-val[Z]+a[Z][Z])))/
					  (a[X][Z]*(-a[X][Y]*a[Y][Z]+a[X][Z]*(a[Y][Y]-val[Z])));
		vect[Y][Z] = -(a[X][Z]*a[Y][Z]-a[X][Y]*(-val[Z]+a[Z][Z]))/
					  (-a[X][Y]*a[Y][Z]+a[X][Z]*(a[Y][Y]-val[Z]));
		vect[Z][Z] = 1.0;
	}

	// Normalize the eigenvectors
	double mag;

	mag = vect[X][X]*vect[X][X] + vect[Y][X]*vect[Y][X] + vect[Z][X]*vect[Z][X];
	
	mag = sqrt(mag);

	if(mag < eps)
		return false;

	mag = 1.0/mag;

	vect[X][X] *= mag;
	vect[Y][X] *= mag;
	vect[Z][X] *= mag;

	mag = vect[X][Y]*vect[X][Y] + vect[Y][Y]*vect[Y][Y] + vect[Z][Y]*vect[Z][Y];

	mag = sqrt(mag);

	if(mag < eps)
		return false;

	mag = 1.0/mag;
	
	vect[X][Y] *= mag;
	vect[Y][Y] *= mag;
	vect[Z][Y] *= mag;

	mag = vect[X][Z]*vect[X][Z] + vect[Y][Z]*vect[Y][Z] + vect[Z][Z]*vect[Z][Z];

	mag = sqrt(mag);

	if(mag < eps)
		return false;

	mag = 1.0/mag;
	
	vect[X][Z] *= mag;
	vect[Y][Z] *= mag;
	vect[Z][Z] *= mag;

	return true;
}
 

// Compute the eigenvalues for a 3x3 symmetric real matrix
void EigenValues3D(double a[3][3], double val[3])
{
	// Trace of the matrix a
	double trace_a = a[X][X] + a[Y][Y] + a[Z][Z];

	// Determinant of the matrix a
	double det_a = -a[X][Z]*a[X][Z]*a[Y][Y] + 2*a[X][Y]*a[X][Z]*a[Y][Z] - 
					a[X][X]*a[Y][Z]*a[Y][Z] - a[X][Y]*a[X][Y]*a[Z][Z] + 
					a[X][X]*a[Y][Y]*a[Z][Z];

	// Determinant of the cofactors of the matrix a
	double det_c = -a[X][Y]*a[X][Y] - a[X][Z]*a[X][Z] + a[X][X]*a[Y][Y] - 
					a[Y][Z]*a[Y][Z] + a[X][X]*a[Z][Z] + a[Y][Y]*a[Z][Z];

	// Intermediate constants
	complex<double> foo = 3.0*det_c - trace_a*trace_a;
	complex<double> bar = 27.0*det_a - 9.0*det_c*trace_a + 2.0*trace_a*trace_a*trace_a;

	complex<double> z = bar + sqrt(bar*bar + 4.0*foo*foo*foo);
	complex<double> z_1_3 = pow(z, 1.0/3.0);

	const complex<double> I(0, 1);
	const double two_1_3 = pow(2.0, 1.0/3.0);
	const double two_2_3 = pow(2.0, 2.0/3.0);
	const double three_1_2 = sqrt(3.0);

	complex<double> tmp;
	
	// The eigenvalues
	tmp = trace_a/3.0 - two_1_3*foo/(3.0*z_1_3) + z_1_3/(3.0*two_1_3);
	val[X] = tmp.real();

	tmp = trace_a/3.0 + (1.0 + I*three_1_2)*foo/(3.0*two_2_3*z_1_3) - 
		(1.0 - I*three_1_2)*z_1_3/(6.0*two_1_3);
	val[Y] = tmp.real();

	tmp = trace_a/3.0 + (1.0 - I*three_1_2)*foo/(3.0*two_2_3*z_1_3) - 
		(1.0 + I*three_1_2)*z_1_3/(6.0*two_1_3);
	val[Z] = tmp.real();
}