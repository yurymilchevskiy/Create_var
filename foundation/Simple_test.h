#ifndef SIMPLE_TEST_H
#define SIMPLE_TEST_H

#include	<string>

class Simple_test
{
public:
	Simple_test ();
	virtual ~Simple_test () ;
	void run_test();
protected:
    virtual void    run() = 0;
	void            test_( const std::string & condition_str,	bool condition);
	inline int		get_passed() {return passed_;}
	inline int		get_failed() {return failed_;}
private:
	int		passed_;
	int     failed_;
};


#endif