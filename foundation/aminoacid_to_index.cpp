#include "aminoacid_to_index.h" 
#include <cassert>

int aminoacid_to_index ( const char aa )
{
	static	map <char, int> aa_index = set_index_by_aa();
	if ( aa_index.find(aa) != aa_index.end()  ) 
	return  aa_index  [ aa  ];
	else 
	{
//		assert(0);
		return -1;
	}
	
}

map <char, int>  set_index_by_aa() 
{
	map < char, int >  aa_index ;

    aa_index  [  'A'  ] = 0  ;
    aa_index  [  'C'  ] = 1  ;
    aa_index  [  'D'  ] = 2  ;
    aa_index  [  'E'  ] = 3  ;
    aa_index  [  'F'  ] = 4  ;
    aa_index  [  'G'  ] = 5  ;
    aa_index  [  'H'  ] = 6  ;
    aa_index  [  'I'  ] = 7  ;
    aa_index  [  'K'  ] = 8  ;
    aa_index  [  'L'  ] = 9  ;
    aa_index  [  'M'  ] = 10  ;
    aa_index  [  'N'  ] = 11  ;
    aa_index  [  'O'  ] = 12  ;
    aa_index  [  'P'  ] = 13  ;
    aa_index  [  'Q'  ] = 14  ;
    aa_index  [  'R'  ] = 15  ;
    aa_index  [  'S'  ] = 16  ;
    aa_index  [  'T'  ] = 17  ;
    aa_index  [  'V'  ] = 18  ;
    aa_index  [  'W'  ] = 19  ;
    aa_index  [  'X'  ] = 20  ;
    aa_index  [  'Y'  ] = 21  ;


	return aa_index ;
}

char index_to_aminoacid ( const int  index )
{
	static	map < int, char >  index_aa =  set_aa_by_index () ;
	if ( index_aa  .find(index) != index_aa .end()  ) 
	return  index_aa [ index  ];
	else 
	{
//		assert(0);
		return -1;
	}
	
}


map <int , char >  set_aa_by_index () 
{
	map < int, char >  index_aa ;

    index_aa [ 0   ] ='A' ;
    index_aa [ 1   ] ='C' ;
    index_aa [ 2   ] ='D' ;
    index_aa [ 3   ] ='E' ;
    index_aa [ 4   ] ='F' ;
    index_aa [ 5   ] ='G' ;
    index_aa [ 6   ] ='H' ;
    index_aa [ 7   ] ='I' ;
    index_aa [ 8   ] ='K' ;
    index_aa [ 9   ] ='L' ;
    index_aa [ 10  ] ='M'  ;
    index_aa [ 11  ] ='N'  ;
    index_aa [ 12  ] ='O'  ;
    index_aa [ 13  ] ='P'  ;
    index_aa [ 14  ] ='Q'  ;
    index_aa [ 15  ] ='R'  ;
    index_aa [ 16  ] ='S'  ;
    index_aa [ 17  ] ='T'  ;
    index_aa [ 18  ] ='V'  ;
    index_aa [ 19  ] ='W'  ;
    index_aa [ 20  ] ='X'  ;
    index_aa [ 21  ] ='Y'  ;


	return index_aa;
}


int get_size_aminoacid_set () 
{
	static	map <char, int> aa_index = set_index_by_aa();
	return  aa_index.size(); 
}


bool is_standard_aa ( const char aa )
{
	int index = aminoacid_to_index (aa) ;
	return ( index  >= 0 &&  index <= 21 ) ;
}

int get_virtual_residue_index ()
{
	return aminoacid_to_index ('O');
}

//get_aa_size ()