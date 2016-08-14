#ifndef C_FOURIERSMOOTHED_H
#define C_FOURIERSMOOTHED_H

#ifndef PROPERTY_H
#include "Property.h"
#endif

#include <string>
#include <vector>

using namespace std;

class  c_FourierSmoothed: public Property
{
public:

	c_FourierSmoothed(CowardVariables & cowa_store) :
		Property(cowa_store)
	{}

	explicit c_FourierSmoothed(
		CowardVariables & cowa_store,
		const string	& task_string);


	Property* clone(
		const string & task_string) const;
		double    calc_value(const int   position_in_chain);

protected:


	string	name_of_property_ ;  // это лишняя переменная - потом убрать !!

	string	variable_name_in_list_; // имя переменной для вызовов внутри списка

	int		left_border_ ;
	int		right_border_ ;
	int		point_number_;

	double		period_lenfth_;
	int			period_number_;

	double		power_;

	int		property_ID_;

	c_FourierSmoothed(const c_FourierSmoothed&);
	c_FourierSmoothed& operator = (const c_FourierSmoothed&);
};

#endif