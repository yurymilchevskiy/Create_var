#include "Curly_calculus_system_test.h"
#include "Curly_calculus_system.h"
#include <iostream>

Curly_calculus_system_test ::
~Curly_calculus_system_test ()
{
	cout << "Curly_calculus_system_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}

void Curly_calculus_system_test::
single_enough_test ()
{
	vector < int >  bases;
	bases.push_back(2);
	bases.push_back(3);
	bases.push_back(4);
	bases.push_back(5);

	Curly_calculus_system ob (bases);

	int number_of_elements = ob.get_number_of_elements ();

	for ( int ii=0; ii<number_of_elements; ii++) 
	{
		int cursor = ii;

		vector < int > array = ob.get_array_by_cursor  (cursor);

		for ( int jj=0; jj< array.size(); jj++) 
			cout << array [jj] << "\t";
		cout << endl;

		int inverse_cursor = ob.get_cursor_by_array    ( array ) ; 

		test_( "check get_array_by_cursor  & get_cursor_by_array  consistency",	cursor 	== inverse_cursor 	);

	}
}

void Curly_calculus_system_test::
init_through_pointer_test ()
{
	vector < int >  bases;
	bases.push_back(2);
	bases.push_back(3);
	bases.push_back(4);
	bases.push_back(5);

	Curly_calculus_system *ob = new  Curly_calculus_system (bases);


	int number_of_elements = ob->get_number_of_elements () ;

	for ( int ii=0; ii<number_of_elements; ii++) 
	{
		int cursor = ii;

		vector < int > array = ob->get_array_by_cursor  ( cursor);

		for ( int jj=0; jj< array.size(); jj++) 
			cout << array [jj] << "\t";
		cout << endl;

		int inverse_cursor = ob->get_cursor_by_array    ( array ) ; 

		test_( "check get_array_by_cursor  & get_cursor_by_array  consistency",	cursor 	== inverse_cursor 	);
	}

	delete  ob;


}
