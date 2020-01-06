#include "headfile.h"
#include "global.h"
#include "myfun.h"
#include "mem_loc.h"
#define PAI 3.1415926535897

float mf8(FMATRIX loca, int a)
{
	float result;

	float x1, x2;
	x1=loca[0][a]+5;
	x2=loca[1][a]+5;

	result=pow(sin(2*PAI*x1),3)*sin(2*PAI*x2)/(pow(x1,3)*(x1+x2));;
	return (float) result;

}

float f6(FMATRIX loca, int a)
{
/*
	This is the f6 function as described in the Handbook of
	Genetic Algorithms, p.8
*/
	float num, denom, f6;
	float errorf6;

	num=(sin(sqrt((loca[0][a]*loca[0][a])+(loca[1][a]*loca[1][a]))))  *
		 (sin(sqrt((loca[0][a]*loca[0][a])+(loca[1][a]*loca[1][a])))) - 0.5;
	denom=(1.0 + 0.001 * ((loca[0][a] * loca[0][a]) + (loca[1][a]*loca[1][a]))) *
		(1.0 + 0.001 * ((loca[0][a] * loca[0][a]) + (loca[1][a]*loca[1][a])));
	f6= (float) 0.5 - (num/denom);
	errorf6=1 - f6;
	return (float)errorf6;
}

float sphere(FMATRIX loca, int a, int b)
{
	/* This is the familiar sphere model
		int a: index of particles   b:dimension */

	float result;
	int i;

	result=0.0;

	for (i=0;i<b;i++)
	{
		result += loca[i][a]*loca[i][a];
	}

	return (float)result;
}

float rosenbrock(FMATRIX loca, int a, int b)
{
	/* this is the Rosenbrock function
		a: index of the particles; b:dimension */

	int i;
	float result;

	result=0.0;

	for (i=1;i<b;i++)
	{
		result +=100.0*(loca[i][a]-loca[i-1][a]*loca[i-1][a])*(loca[i][a]-loca[i-1][a]*loca[i-1][a]) + (loca[i-1][a]-1)*(loca[i-1][a]-1);
	}

	return (float)fabs(result);
}

float rastrigrin(FMATRIX loca, int a, int b)
{
	/* This is the generalized Rastrigrin function
		a:index of the particles; b:dimension */

	int i;
	float result;

	result=0.0;

	for (i=0;i<b;i++)
	{
		result +=loca[i][a]*loca[i][a] - 10.0*cos(2.0*PAI*loca[i][a])+10.0;
	}

	return (float)result;
}

float griewank(FMATRIX loca, int a,int b)
{
	/* This is the generalized Griewank function
		a:index of the particles; b:dimension */

	int i;
	float result_s,result_p;

	result_s=0.0;
	result_p=1.0;

	for (i=0;i<b;i++)
	{
		result_s +=loca[i][a]*loca[i][a];
		result_p *=cos(loca[i][a]/sqrt(i+1));
	}
	result_s =result_s/4000.0 - result_p +1;

	return (float)result_s;
}

float ackley(FMATRIX loca, int a,int b)
{
	/* This is the ackley function
		a:index of the particles; b:dimension */

	int i;
	float result_s,result_p;

	result_s=0.0;
	result_p=0.0;

	for (i=0;i<b;i++)
	{
		result_s +=loca[i][a]*loca[i][a];
		result_p +=cos(2*PAI*loca[i][a]);
	}
	result_s = -20*exp(-0.2*sqrt(result_s/b))-exp(result_p/b)+20+exp(1);

	return (float)result_s;
}

float mf2(FMATRIX locabefore, int a, int b)
{
	int i;
	float result;

	float result_opt, value1,value2, value3;

	float result_cons1;
	float result_cons2;

	FVECTOR loca;
//	printf("test");

	FVectorAllocate(&loca,b);

    for (i=0; i<b; i++) {
	  
      loca[i] = locabefore[i][a]+5;
    }

	//solve optimum function
	value1=0;
	value2=1;
	value3=0;

    for (i=0; i<b; i++) {
      value1 += pow(cos(loca[i]),4.0);
    }
    for (i=0; i<b; i++) {
      value2 *= pow(cos(loca[i]),2.0);
    }
    for (i=0; i<b; i++) {
      value3 += i*pow(loca[i],2.0);
    }
    result_opt = -1*fabs((value1-2.0*value2)/sqrt(value3));

	//solve first cons function
    result_cons1 = 1;
    for (i=0; i<b; i++) {
      result_cons1 *= loca[i];
    }
    result_cons1 -= 0.75;

	//solve second cons function
    result_cons2 = 0;
    for (i=0; i<b; i++) {
      result_cons2 += loca[i];
    }
    result_cons2 -= 7.5*b;

	if (result_cons1<0 || result_cons2>0) {
		result = 0.75-result_cons1+150+result_cons2;
	} else {
		result = result_opt;
	}

//	printf("%f, %f, %f, %f\n", result_opt, result_cons1, result_cons2, result);

	return (float)result;
}

