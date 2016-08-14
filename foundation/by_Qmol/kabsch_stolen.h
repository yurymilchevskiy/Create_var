#ifndef KABSCH_STOLEN_H
#define KABSCH_STOLEN_H

	void kabsch_stolen(double *ref, double *vec, double *w, const int num, bool & error_flag );
	double rmsd(double *ref, double *vec, double *w, const int num);
	double kabsch_rmsd (double *ref, double *vec, double *w, const int num, bool & error_flag);

#endif
