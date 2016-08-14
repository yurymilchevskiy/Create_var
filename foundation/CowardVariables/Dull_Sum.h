#ifndef DULL_SUM_H
#define DULL_SUM_H

#ifndef PROPERTY_H
#include "Property.h"
#endif

#include <string>
#include <vector>

using namespace std;

class  Dull_Sum: public Property
{
public:

	/*Dull_Sum(vector < Frequency_chain_constants > & pull_fcc,
		vector < vector < double > >   & Chain_Prime_Constants) :
		Property(pull_fcc, Chain_Prime_Constants) {};*/

	Dull_Sum(CowardVariables & cowa_store ) :
		Property( cowa_store) 		
	{}

	explicit Dull_Sum(
		CowardVariables & cowa_store,
		const string	& task_string);

	//explicit Dull_Sum		( const string & task_string  );

	Property* clone	( 
		const string & task_string) const;

	double    calc_value(const int   position_in_chain);

protected:
	
	string	variable_name_in_list_; // имя переменной для вызовов внутри списка

	int		left_border_ ;
	int		right_border_ ;

	double  power_ ;

	int		property_ID_;

	char	fabs_mode_;

	Dull_Sum(const Dull_Sum&);
	Dull_Sum& operator = (const Dull_Sum&);
};

#endif

