#pragma warning( disable : 4786 )

#include "Fragment_base_subtle_test.h"
#include "Fragment_base_subtle.h"

#include "../CommonFunc.h"

#include <fstream>

using namespace std;

extern ofstream log_stream;

Fragment_base_subtle_test::
~Fragment_base_subtle_test()
{
	cout << "Fragment_base_subtle_test PASSED: " <<"  failed: " << get_failed() <<  "	passed: " << get_passed() << endl;
}
void Fragment_base_subtle_test::
fill_up_test ()
{
	//Fragment_base_subtle ob("5a",FRAGMENT_BASE_SUBTLE_FILL_UP);
	Fragment_base_subtle ob("10a",FRAGMENT_BASE_SUBTLE_FILL_UP);

	int total_fragment_number = ob.get_total_fragment_number();
}
void Fragment_base_subtle_test::
init_test  ()
{
//	Fragment_base_subtle ob("5a",FRAGMENT_BASE_SUBTLE_COMMON_USAGE);
//	Fragment_base_subtle ob("4a",FRAGMENT_BASE_SUBTLE_COMMON_USAGE);
	Fragment_base_subtle ob("10a",FRAGMENT_BASE_SUBTLE_COMMON_USAGE);
	cout << "Ok" << endl;
	int total_fragment_number = ob.get_total_fragment_number();
	//cout << "total_fragment_number for tetrapeptide:" << total_fragment_number << endl;
	cout << "total_fragment_number for set of 10-mer:" << total_fragment_number << endl;

}
