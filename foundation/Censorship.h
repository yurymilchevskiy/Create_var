// Class Censorship enables Restriction for acceptig DSSP data

#ifndef SCENSORSHIP_H
#define SCENSORSHIP_H
/*
enum SelectedChainStatus            
{
	ALL_CHAINS,
	STATED_SINGLE_CHAIN,
	STATED_SUBSET_OF_CHAIN
} ;

enum NonStandardResidueStatus
{
	PROHIBITED_NON_STANDARD_RESIDUE_AT_ALL,
	PROHIBITED_NON_STANDARD_RESIDUE_FOR_CURRENT_CHAIN,
	PERMITTED_NON_STANDARD_RESIDUE                   
} ;


enum ExperimentalMethodStatus
{
	PERMITTED_ANY_EXPERIMENTAL_METHOD,
	PERMITTED_X_RAY_ONLY,
	PERMITTED_NMR_ONLY
} ;
*/
/*enum Validity_Resume
{
	INVALID_UNDER_CURRENT_RESTRICTION,
	VALID_SUBSET_OF_CHAIN_ONLY,
	VALID_WHOLE_PROTEIN
};
*/


#include <string>
#include <fstream>
#include <map>

using namespace std; 

class  Censorship
{
public:
	Censorship():is_initialysed_yet(false) {};
	~Censorship() {}

    void display_CensorshipOption_setted (  const string & option_show_file_name ) const ;

	const string  & option_meaning		 ( const string & key )  const;
	
	void init ( const string & parm_source_file_name );

private:

	map < string, string  > CensorshipOption_; 
	bool is_initialysed_yet;

	Censorship (const Censorship&);
	void operator = (const Censorship &);

};

#endif
