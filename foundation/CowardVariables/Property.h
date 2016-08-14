#ifndef PROPERTY_H
#define PROPERTY_H

#include	<string>
#include	<vector>
#include	<map>



using namespace std;

#include "../Frequency_extrapolation/Frequency_chain_constants.h"

class Frequency_extrapolation ;
class CowardVariables;

class  Property
{
public:

	Property(CowardVariables & cowa_store):
		cowa_store_ (cowa_store)
	{}
	  
	virtual			~Property() {}; // тут вставил тело ф-ции  {} чтобф исчезли ошибки unresolved external 
	//virtual			Property* clone (
	//	const string & task_string, 
	//	map    < string, int >	& frequency_name_to_pull_index,
	//	map   < string, int >   * co_task_variable_name_to_index,
	//	map    < string, int >	* frequency_name_to_pull_index ) const					= 0;


	virtual Property* clone(const string & task_string	) const = 0;
	virtual		double   calc_value ( const int   position_in_chain  )   = 0;

	//double			get_value		() const { return value_; } 

	//bool			is_subsidiary () const { return is_subsidiary_; }

protected:	

	CowardVariables & cowa_store_;

	//double	value_;
	string	name_;
	//bool    is_subsidiary_;
};

#endif