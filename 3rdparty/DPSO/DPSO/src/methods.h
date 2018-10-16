#ifndef __METHODS_H__
#define __METHODS_H__

#include "define_t.h"

float meanTotalValues(FMATRIX , int , int );
FVECTOR meanValues(FMATRIX , int , int );
float squareSum(FVECTOR, int);
void outputMatrix(FILE*, FMATRIX, int, int);
void printMatrix(FMATRIX, int, int);
void outputVECTOR(FILE* , FVECTOR , int );
void printVECTOR(FVECTOR , int );
float totalDiversity(FMATRIX, int, int);
float totalMeans(FMATRIX , int , int );

#endif
