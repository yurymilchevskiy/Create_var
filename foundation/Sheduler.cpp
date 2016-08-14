#pragma warning( disable : 4786 )

#include "Sheduler.h"

#include "TaskSet.h"

#include "CommonFunc.h"

extern ofstream log_stream;


Sheduler ::
Sheduler ( const string & parm_source_file_name  ) 
{
	init ( parm_source_file_name );
}
void Sheduler::
init ( const string & parm_source_file_name )
{
	ifstream  parm_stream( parm_source_file_name.c_str() );
	if ( ! parm_stream  )	
	{	
		cout		<< "Sheduler: file " << parm_source_file_name   << " not found "  << endl;
		log_stream	<< "Sheduler: file " << parm_source_file_name  << " not found "  << endl;
		exit (1);
	}
	 taskset_get_option ( parm_stream,sheduler_option_  ) ;
}

void Sheduler::
display_sheduler_option_setted (
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
	for ( theIterator = sheduler_option_.begin();theIterator != sheduler_option_.end();theIterator ++ ) 
	{
		string key		=	(*theIterator).first ;
		string value	=   (*theIterator).second ;

		PutVa(key,		show_stream,25,24,'l');
		PutVa(value,	show_stream,35,34,'l');
		show_stream << endl;

	}
}

const string  & Sheduler::
option_meaning ( const string & key ) const
{
	//string  value = sheduler_option_[key ];

    
    map < string, string  > ::const_iterator  theIterator;

	static const string unknown_string ("UNKNOWN");

	theIterator = sheduler_option_.find(key) ;

	if ( theIterator !=  sheduler_option_.end() ) 
		return theIterator->second ;
	else 
		return unknown_string;
}


