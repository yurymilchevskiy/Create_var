#pragma warning( disable : 4786 ) 

#include <cstring>
#include <cstdio>
#include <iostream>
//#include <iomanip>
#include "CommonFunc.h"


#define MAX_WRD_LEN   100

void PutVaDouble (double Var,ostream &output,int Length, int Precision,char AlignParm) {

	char *Buff = new  char [MAX_WRD_LEN];
	memset (Buff,0,MAX_WRD_LEN*sizeof(char));
	char *Resu = new  char [MAX_WRD_LEN];
	memset (Resu ,0,MAX_WRD_LEN*sizeof(char));


	sprintf(Buff,"%16.12lf",Var);

	int len = strlen (Buff);

	int OnFlag = 0,OffFlag =0,PecisionFlag=0;
	int kk=0;
	for(int ii=0;ii<len ;ii++) {
		if (Buff[ii] == ' ')	continue;
		
		Resu[kk] = Buff[ii];
		
		if (Buff[ii]=='.') 
			OnFlag = 1;
		if (OnFlag)
			PecisionFlag++;
		
		if (PecisionFlag == Precision+2 ) 
			OffFlag = 1;
		
		if (OffFlag)
			Resu[kk] = '\0' ;
		kk++;
	}

	int LenRes = strlen (Resu);

	char *Final = new  char [MAX_WRD_LEN];
	memset (Final,0,MAX_WRD_LEN*sizeof(char));

	LenRes = (LenRes < Length) ? LenRes : Length;
	
	switch (AlignParm) {
		case 'L': case 'l':
			memcpy(Final,Resu,LenRes);
			memset(Final+LenRes,' ',Length-LenRes);
			
			break;
		case 'R': case 'r':
			memcpy(Final+Length-LenRes,Resu,LenRes);
			memset(Final,' ',Length-LenRes);
			break;
	}



	output << Final;


	delete [] Buff ;	delete [] Resu ;	delete [] Final;
}
