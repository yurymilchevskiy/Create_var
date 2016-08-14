#include "Simple_test.h"
#include	<iostream>
#include	<string>

using namespace std;

Simple_test::Simple_test() : passed_(0), failed_(0)  { }
Simple_test::~Simple_test() {  passed_ = 0; failed_ = 0; }

void Simple_test::
run_test()
{
	run();
}
void   Simple_test::  
test_( const string & condition_str,	bool condition)
{
	if(condition)
		++passed_;
	else
	{
		++failed_;
		cout << "TEST NOT PASSED: " << condition_str << endl;
	}
}