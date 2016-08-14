#include <cmath>
#include <iostream>

#include "special_function.h"

using namespace std;

	double betai(double a, double b, double x);

double prob_by_student (double t, int njue)
{
	double d_njue = (double) njue;
			
	double x = d_njue /(d_njue + t*t);
	double a = d_njue/2;
	double b = 0.5;
	
	double check = betai(a, b, x);

//	return 1.0 - check ;
	
	return ( t >= 0 ) ? 1 - check/2 : check/2 ;


}


double betai(double a, double b, double x)
//Returns the incomplete beta function Ix(a, b).
{
	double betacf(double a, double b, double x);
	double gammln(double xx);
	void nrerror(char error_text[]);
	double bt;

	if (x < 0.0 || x > 1.0) 
		cout << "Bad x in routine betai" << endl;
	if (x == 0.0 || x == 1.0) bt=0.0;
	else 
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) 
		return bt*betacf(a,b,x)/a;
	else 
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}


#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
double betacf(double a, double b, double x)
//Used by betai: Evaluates continued fraction for incomplete beta function by modi.ed Lentz’s
//method (§5.2).
{
	void nrerror(char error_text[]);
	double aa,c,d,del,h,qab,qam,qap;

	int m,m2;

	qab=a+b; //These q’s will be used in factors that occur	in the coe.cients (6.4.6). 
	qap=a+1.0;
	qam=a-1.0;
	c=1.0; //First step of Lentz’s method.
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d; //One step (the even one) of the recurrence.
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d; //Next step of the recurrence (the odd one).
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
	if (fabs(del-1.0) < EPS) break; // Are we done?
	}
	if (m > MAXIT) cout << "a or b too big, or MAXIT too small in betacf" << endl;
	return h;
}




double gammln(double xx)
//Returns the value ln[Ã(xx)] for xx > 0.
{
//Internal arithmetic will be done in double precision, a nicety that you can omit if .ve-.gure
//accuracy is good enough.
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,	24.01409824083091,-1.231739572450155,	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) 
		ser += cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}