#ifndef FILL_UP_FRAGMENT_TORSION_ANGLES_H
#define FILL_UP_FRAGMENT_TORSION_ANGLES_H

#include < vector > 

using namespace std;

void fill_up_fragment_torsion_angles ( 
	 const double * coord_array,
	 const int length, 
	 vector < double > & torsion_set, 
	 const char radian_or_gradus_flag  );

#endif