#ifndef SHEDULER_H
#define SHEDULER_H

#include <string>
#include <fstream>
#include <map>
;
using namespace std; 

class  Sheduler
{
public:

	~Sheduler() {};
	Sheduler ( const string & parm_source_file_name  ) ;
	const string  & option_meaning ( const string & key )  const;
	
	void init ( const string & parm_source_file_name );
    void display_sheduler_option_setted (  const string & option_show_file_name ) const ;
private:
	map < string, string  > sheduler_option_; 

	bool is_initialysed_yet;

	Sheduler (const Sheduler&);
	void operator = (const Sheduler &);
};

#endif
	


