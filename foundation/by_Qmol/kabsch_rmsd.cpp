
#include "kabsch_stolen.h"
#include "../Geometry_util/Fragmrent_tranformation.h"

double kabsch_rmsd (double *ref, double *vec, double *w, const int num, bool & error_flag)
{
	centroid ( ref, num/3); 
	centroid ( vec, num/3);

	kabsch_stolen	(ref, vec, w, num, error_flag );

	if (error_flag)
		return -1;
	
	double RMSD = rmsd	(ref, vec, w, num );
	return RMSD;
}