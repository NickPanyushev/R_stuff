#ifndef _MATRIX_H
#define _MATRIX_H

typedef struct {
	int **vals;
	int numrows;
	int numcols;
}IntMat;

typedef struct {
	double **vals;
	int numrows;
	int numcols;
}DblMat;

extern IntMat* constructIntMat(int numrows, int numcols);
extern DblMat* constructDblMat(int numrows, int numcols);
extern void destructIntMat(IntMat*);
extern void destructDblMat(DblMat*);

#endif

