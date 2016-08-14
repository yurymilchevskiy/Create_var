#ifndef FILL_UP_ACCEPTED_CHAIN_DATA_TEST_H
#define FILL_UP_ACCEPTED_CHAIN_DATA_TEST_H

#ifndef SIMPLE_TEST_H
#include "Simple_test.h"

#endif

#include <string>

class  Fill_up_accepted_chain_data_test: public Simple_test
{
public:
	~Fill_up_accepted_chain_data_test();

    void run()
    {
		first_simple_test();
	}
	void first_simple_test();

};

#endif
