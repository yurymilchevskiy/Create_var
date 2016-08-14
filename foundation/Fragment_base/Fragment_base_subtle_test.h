#ifndef FRAGMENT_BASE_SUBTLE_TEST_H
#define FRAGMENT_BASE_SUBTLE_TEST_H

#ifndef SIMPLE_TEST_H
#include "../Simple_test.h"

#include "../CommonFunc.h"

#include <iostream>

using namespace std;
#endif

#include <string>

class  Fragment_base_subtle_test: public Simple_test
{
public:
	~Fragment_base_subtle_test();

    void run()
    {
		fill_up_test	();
//		init_test		();
	}
	void fill_up_test	();
	void init_test		();
};

#endif
