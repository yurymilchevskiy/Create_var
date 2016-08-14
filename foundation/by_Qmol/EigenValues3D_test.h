#ifndef EIGENVALUES3D_TEST_H
#define EIGENVALUES3D_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  EigenValues3D_test: public Simple_test
{
public:
	~EigenValues3D_test();

    void run()
    {
		eigenvalues3d_test ();
		EigenSystem3D_t ();
	}
	void eigenvalues3d_test();
	void EigenSystem3D_t ();
};

#endif
