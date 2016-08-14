#ifndef CURLY_CALCULUS_SYSTEM_TEST_H
#define CURLY_CALCULUS_SYSTEM_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  Curly_calculus_system_test : public Simple_test
{
public:
	~Curly_calculus_system_test ();

    void run()
    {
		single_enough_test ();
		init_through_pointer_test ();
	}
	void single_enough_test ();
	void init_through_pointer_test ();
	};

#endif
