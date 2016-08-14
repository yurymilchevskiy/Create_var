#ifndef LADDER_CATCH_H 
#define LADDER_CATCH_H 

#ifndef BASE_COWARD_H
#include "Base_coward.h"
#endif

#include <string>
#include <vector>

using namespace std;


class  Ladder_catch: public Base_coward
{
public:

	Ladder_catch () {} ;

	explicit Ladder_catch ( 	const string			& task_string,
				map   < string, int	>	&	co_task_variable_name_to_index );

	//explicit Dull_Sum		( const string & task_string  );

    Base_coward*			clone	( const string & task_string, map   < string, int	>	&	co_task_variable_name_to_index  ) const;

	void    calc_value ( 
		const int   position_in_chain,  
			  int	var_set_cursor, 
		const		vector < vector < double > >   & Chain_Prime_Constants,
					vector < vector < double > >   & sophisticated_variables    )  ;
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

