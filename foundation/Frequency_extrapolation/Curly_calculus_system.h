#ifndef CURLY_CALCULUS_SYSTEM_H
#define CURLY_CALCULUS_SYSTEM_H

#include <vector>

using namespace std;
using namespace std;

class Curly_calculus_system
{
public:
	Curly_calculus_system ( const vector < int > & base);
//	~Curly_calculus_system ( ) 

	vector < int >	get_array_by_cursor		( int cursor );
	int				get_cursor_by_array		( const vector <int> array );
	int				get_number_of_elements	() const { return number_of_elements_; }  ;

private:
	vector < int >  base_;
	int number_of_elements_;
};



#endif

