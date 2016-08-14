	//*************************************************************************
//*** Содержит ф-ии общего назначения 
// ***************************************************************************
#ifndef COMMON_FUNC_H
#define COMMON_FUNC_H

// Шаблон позволяет записывать переменные любых встроенных типов в любой выходной поток 
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <map>
#include <vector>

#define max(a,b) ((a>=b) ? a : b)
#define min(a,b) ((a <b) ? a : b)

using namespace std;

	inline double	Pythagorean_Number () { static double value = 4*atan (1.0); return value ; } 
	inline double	epsilon_float(){ return ( std::numeric_limits<float>::epsilon()  );}

	inline double	large_float  (){ return  1/epsilon_float() ;}

template <class T >
//void PutVa (T  Var,ofstream &output,int Length,> std::_Ios_Fmtflags  Precision,char AlignParm) {
void PutVa (
 T  Var,
 ostream &output,
 int Length, 
 int Precision,
 char AlignParm) {
//	output.flags(0) ;
	//long Align = (AlignParm=='l' || AlignParm=='L') ? ios::left : ios::right ;
//	std::_Ios_Fmtflags Align = (AlignParm=='l' || AlignParm=='L') ? ios::left : ios::right ;

	output << setw(Length) << setiosflags( ( AlignParm=='l' || AlignParm=='L') ? ios::left : ios::right ) << setprecision (Precision) << Var;
//	    setw(Length)  << 
//	    setiosflags(Align) << 
//	    setprecision(Precision) << 
//	    Var;
}

void PutVaDouble (double Var,ostream &output,int Length, int Precision,char AlignParm);

template <class T >
bool read_value_from_position_in_fixed_format_string (T & value, const string & line, const int start_position, const int width)
{
	string target = line.substr (start_position, width);
	istringstream ost (target);

	if ( ost >> value ) 
		return true ;
	else 
		return false;
}

// Returns positon of string origin that contains template
// if template was not found returns -1; 
int  find_word_position_in_file (ifstream & input_file,const string shablon );
// overloaded function contains last read line
int  find_word_position_in_file (ifstream & input_file,const string shablon, string  & current_line);


string new_extension_file_name	( const string & old_name, const  string & new_extension);
string new_extension_file_name_from_tail ( const string & old_name, const  string & new_extension);

string cut_extension_file_name	( const string & name );

string get_extension_file_name	( const string & name );
string get_extension_file_name_from_tail (const string & source_file_nake);

string get_base_file_name		( const string & source_file_nake);

bool get_line_from_stream_and_compare_to_templare 	( istream &input, const string template_string ) ;

// returns true if residue_name correspond to standard aminoacid residue
// UPPERCASE letter used for residue_name 
//bool is_standard_residue( string & residue_name );

map < string, string  >    Suck_up_options ( string & option_file_name);

// Calculates combination number of N by K
int Combination_N_by_K (const int N, const int K);

// Calculates correlation for vectors v1 & v2
double    calculate_correlation(vector <double> & v1,vector <double> & v2);
double    calculate_correlation(const double * v1,const double * v2,const int casenum);

vector < int > translate_to_another_calculus_system (const int base,const int length,const int cursor );
vector < long > translate_to_another_calculus_system_long (const int base,const int length,const long int cursor );

// перемнож. матр. m1 и m2 и помещает рез-т в mR;
//r1 - число стр. m1, c2 число столб. m2 ,c1_r2 - число столб. m1 = число стр. m2
void  MatrixMultiplication(
	const int r1,
	const int c1_r2,
	const int c2,\
	const vector <vector < double> > & m1,
	const vector <vector < double> > & m2, 
	vector <vector < double>  > & mR) ;

//string get_path_by_keyword (const string & OptionsFileName, const string & keyword);

// читает options из файла: KEY_WORD opt1 opt2 .... 
map < string, vector < string > >    Suck_up_complexive_options ( const string & option_file_name);

string word_toupper ( const string  word );

string word_tolower (const string  word );


// Sapienti sat 
vector < int > prepare_histogram ( 
   vector < double > metrics, 
   const int item_number );

vector < int > prepare_histogram ( 
   vector < double > metrics, 
   const int item_number,
   const vector < double > histogram_scale );

void copy_file ( const string & source , const string & dest ) ;


void print_symmetric_matrix (
	const double *matrix, 
	const int size,
	ofstream		& out );


void print_square_matrix (
	const double *matrix, 
	const int size,
	ofstream	& out );

void print_3x3_matris (
	double m[3][3],
	const int size,
	ofstream	& out );

#endif

