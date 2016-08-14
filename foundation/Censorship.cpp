#pragma warning( disable : 4786 )
	
#include "Censorship.h"

#include "taskset_get_option.h"

#include "CommonFunc.h"

extern ofstream log_stream;

void Censorship::
init ( const string & parm_source_file_name )
{
	ifstream  parm_stream( parm_source_file_name.c_str() );
	if ( ! parm_stream  )	
	{	
		cout		<< "Censorship: file " << parm_source_file_name   << " not found "  << endl;
		log_stream	<< "Censorship: file " << parm_source_file_name  << " not found "  << endl;
		exit (1);
	}
	 taskset_get_option ( parm_stream,CensorshipOption_  ) ;
}

void Censorship::
display_CensorshipOption_setted (
	   const string & path_to_option_show_file_name ) const
{

	ofstream  show_stream ( path_to_option_show_file_name .c_str() );
	if ( ! show_stream )	{	
		cout				<< "Can't create option show file "  << endl;
		log_stream			<< "Can't create option show file "  << endl;
		exit (1);
	}

	typedef map < string, string > MAP_STRING_STRING;
	MAP_STRING_STRING::const_iterator theIterator;

	int counter =0;
	for ( theIterator = CensorshipOption_.begin();theIterator != CensorshipOption_.end();theIterator ++ ) 
	{
		string key		=	(*theIterator).first ;
		string value	=   (*theIterator).second ;

		PutVa(key,		show_stream,25,24,'l');
		PutVa(value,	show_stream,35,34,'l');
		show_stream << endl;

	}
}

const string  & Censorship::
option_meaning ( const string & key ) const
{
    map < string, string  > ::const_iterator  theIterator;

	static const string unknown_string ("UNKNOWN");

	theIterator = CensorshipOption_.find(key) ;

	if ( theIterator !=  CensorshipOption_.end() ) 
		return theIterator->second ;
	else 
		return unknown_string;
}


