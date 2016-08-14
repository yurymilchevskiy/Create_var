#ifndef PROB_BY_STUDENT_TEST_H
#define PROB_BY_STUDENT_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class  Prob_by_student_test: public Simple_test
{
public:
	~Prob_by_student_test();

    void run()
    {
		first_test ();
	}
	void first_test ();
};

#endif
