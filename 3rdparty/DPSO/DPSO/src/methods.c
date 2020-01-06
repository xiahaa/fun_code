#include "headfile.h"
#include "global.h"
#include "mem_loc.h"
#include "myfun.h"
#include "methods.h"

float meanTotalValues(FMATRIX values, int row, int col) {
	float meanValue = 0;
	int i, j;

	j=0;
    for (i=0; i<row; i++) {
	    for (j=0; j<col; j++) {
			meanValue += (fabs(values[i][j]));
		}
    }
	
    return meanValue/((float)(row*col));
}

FVECTOR meanValues(FMATRIX values, int row, int col) {
	FVECTOR meanValues;
	int i, j;
	FVectorAllocate(&meanValues, row);
    for (i=0; i<row; i++) {
	    for (j=0; j<col; j++) {
			meanValues[i] +=  log(fabs(values[i][j]));
		}
		meanValues[i] = meanValues[i]/col;
    }
    return meanValues;
}

float squareSum(FVECTOR v1, int len) {
	float totalSum = 0;
	int i;
	for (i=0; i<len; i++) {
	  totalSum += pow(v1[i],2.0);
	}
	return sqrt(totalSum);
}

float totalDiversity(FMATRIX individuals, int row, int col) {
	FVECTOR meanVal;
	FVECTOR devValues;
	int i, j;
    float diversities = 0;
	FVectorAllocate(&devValues, row);
 	meanVal = meanValues(individuals, row, col);
   for (i=0; i<col; i++) {
		for (j=0; j<row; j++) {
			devValues[j] = individuals[j][i]-meanVal[j];
		}
		diversities += squareSum(devValues, row);
    }
    diversities = diversities/col;
    return diversities;
}

float totalMeans(FMATRIX individuals, int row, int col) {
	int i;
    float diversities = 0;
    for (i=0; i<col; i++) {
		diversities += squareSum(individuals[i], row);
    }
    diversities = diversities/col;
    return diversities;
}

void printMatrix(FMATRIX values, int row, int col) {
	int i, j;
	for (i=0; i<col; i++) {
		for (j=0; j<row; j++) {
			printf("%f\t", values[j][i]);
		}
		printf("\n");
	}
	printf("\n");
}

void outputMatrix(FILE* file, FMATRIX values, int row, int col) {
	int i, j;
	for (i=0; i<col; i++) {
		for (j=0; j<row; j++) {
			fprintf(file, "%f\t", values[j][i]);
		}
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
}
void printVECTOR(FVECTOR values, int col) {
	int j;
	for (j=0; j<col; j++) {
		printf("%f\t", values[j]);
	}
	printf("\n");
}

void outputVECTOR(FILE* file, FVECTOR values, int col) {
	int j;
	for (j=0; j<col; j++) {
		fprintf(file, "%f\t", values[j]);
	}
	fprintf(file, "\n");
}
