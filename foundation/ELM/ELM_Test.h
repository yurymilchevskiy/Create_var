#ifndef ELM_TEST_H
#define ELM_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#endif

#include <string>

class ELM_Test : public Simple_test
{
public:
	~ELM_Test();

	void run()
	{
		dummy_preparing_ELM_data_for_DA();
	}
	void dummy_preparing_ELM_data_for_DA();
};

#endif
