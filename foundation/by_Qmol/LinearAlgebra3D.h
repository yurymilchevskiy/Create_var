#ifndef __LINEAR_ALGEBRA_3D
#define __LINEAR_ALGEBRA_3D

#define		X	0
#define		Y	1
#define		Z	2

// Linear algebra macros:
// A.B (inner product)
#define		DOT(A, B)	(A[X]*B[X] + A[Y]*B[Y] + A[Z]*B[Z])

// AXB (outer product)
#define		CROSS(A,B,C) A[X] = B[Y]*C[Z] - B[Z]*C[Y]; \
						 A[Y] = B[Z]*C[X] - B[X]*C[Z]; \
						 A[Z] = B[X]*C[Y] - B[Y]*C[X]

#define		ROTATE(_RESULT, _MATRIX, _VECTOR)	\
	_RESULT[X] = (_MATRIX)[X][X]*(_VECTOR)[X] + (_MATRIX)[Y][X]*(_VECTOR)[Y] + (_MATRIX)[Z][X]*(_VECTOR)[Z]; \
	_RESULT[Y] = (_MATRIX)[X][Y]*(_VECTOR)[X] + (_MATRIX)[Y][Y]*(_VECTOR)[Y] + (_MATRIX)[Z][Y]*(_VECTOR)[Z]; \
	_RESULT[Z] = (_MATRIX)[X][Z]*(_VECTOR)[X] + (_MATRIX)[Y][Z]*(_VECTOR)[Y] + (_MATRIX)[Z][Z]*(_VECTOR)[Z]

#define		TRANSLATE(_RESULT,_R) \
	_RESULT[X] += (_R)[X]; \
	_RESULT[Y] += (_R)[Y]; \
	_RESULT[Z] += (_R)[Z]

#define		DIFF(_DELTA, _A, _B) \
	_DELTA[X] = (_A)[X] - (_B)[X]; \
	_DELTA[Y] = (_A)[Y] - (_B)[Y]; \
	_DELTA[Z] = (_A)[Z] - (_B)[Z]

#endif // __LINEAR_ALGEBRA_3D