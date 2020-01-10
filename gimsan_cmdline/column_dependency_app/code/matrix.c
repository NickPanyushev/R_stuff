#include "stdinc.h"
#include "matrix.h"

IntMat* constructIntMat(int numrows, int numcols) {
	IntMat *mat = (IntMat*) malloc(sizeof(IntMat));
	mat->numrows = numrows;
	mat->numcols = numcols;
	mat->vals = (int**) malloc(sizeof(int*) * numrows);
	for(int i = 0; i < numrows; i++) {
		mat->vals[i] = (int*) calloc(numcols, sizeof(int));
	}
	return mat;
}

DblMat* constructDblMat(int numrows, int numcols) {
	DblMat *mat = (DblMat*) malloc(sizeof(DblMat));
	mat->numrows = numrows;
	mat->numcols = numcols;
	mat->vals = (double**) malloc(sizeof(double*) * numrows);
	for(int i = 0; i < numrows; i++) {
		mat->vals[i] = (double*) calloc(numcols, sizeof(double));
	}
	return mat;
}

void destructIntMat(IntMat *intmat) {
	for(int i = 0; i < intmat->numrows; i++) {
		free(intmat->vals[i]);
	}
	free(intmat->vals);
	free(intmat);
}


void destructDblMat(DblMat *dblmat) {
	for(int i = 0; i < dblmat->numrows; i++) {
		free(dblmat->vals[i]);
	}
	free(dblmat->vals);
	free(dblmat);
}

