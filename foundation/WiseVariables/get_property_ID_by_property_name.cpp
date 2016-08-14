#pragma warning( disable : 4786 )

#include <string> 
#include <map> 
#include <cassert> 
#include <iostream> 
#include <fstream> 
#include <cstdlib> 

#include "property_constant.h"

using namespace std;

extern ofstream log_stream;

	static	map < string , int>  set_property_index_map ();

int get_property_ID_by_property_name  ( const string  & property_name )
{
	static	map < string , int> property_index_map  = set_property_index_map ();
	if ( property_index_map.find( property_name ) != property_index_map.end()  ) 
		return  property_index_map [ property_name ];
	else 
	{
	//	assert (0);
		cout << property_name << ": Illegal property!" ;
		log_stream << property_name << ": Illegal property!" ;
		exit (-1);
		return -1;
	}
}


map < string , int> set_property_index_map ()
{
	map < string , int> index_map  ;

	for (int ii=0;ii< number_of_all_base_variables  ;ii++) 
		index_map  [ string ( BaseVariablesNames [ii] ) ] = ii;

	return index_map ;

}