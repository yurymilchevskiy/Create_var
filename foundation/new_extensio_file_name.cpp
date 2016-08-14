#pragma warning( disable : 4786 )
 
#include "CommonFunc.h"

string new_extension_file_name ( const string & old_name, const  string & new_extension)
{
	string modyfied_name;

	for (int ii=0; ii< old_name.size(); ii++ ) 
	{
		if (old_name[ii] == '.' )
			break;
		modyfied_name += old_name[ii];
	}
	modyfied_name += '.';
	modyfied_name +=  new_extension;

	return modyfied_name ;
}