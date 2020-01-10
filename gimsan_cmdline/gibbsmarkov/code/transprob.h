#ifndef _TRANSPROB_H
#define _TRANSPROB_H

#include "stdinc.h"

static const int MARKOV_ORDER_BOUND = 5;

typedef struct {
	int count0[NUMCHARS]; //0-th order
	int count1[NUMCHARS][NUMCHARS]; //p(a | b) is given by first[b][a]
	int count2[NUMCHARS][NUMCHARS][NUMCHARS];
	int count3[NUMCHARS][NUMCHARS][NUMCHARS][NUMCHARS];
	int count4[NUMCHARS][NUMCHARS][NUMCHARS][NUMCHARS][NUMCHARS]; //fourth order
	int count5[NUMCHARS][NUMCHARS][NUMCHARS][NUMCHARS][NUMCHARS][NUMCHARS]; //fourth order

	int maxorder; //max order counted/calculated
} TransCount; 

typedef struct {
	double zeroth[NUMALPHAS];
	double first[NUMALPHAS][NUMALPHAS]; //p(a | b) is given by first[b][a]
	double second[NUMALPHAS][NUMALPHAS][NUMALPHAS];
	double third[NUMALPHAS][NUMALPHAS][NUMALPHAS][NUMALPHAS];
	double fourth[NUMALPHAS][NUMALPHAS][NUMALPHAS][NUMALPHAS][NUMALPHAS];
	double fifth[NUMALPHAS][NUMALPHAS][NUMALPHAS][NUMALPHAS][NUMALPHAS][NUMALPHAS];

	int maxorder; //max order counted/calculated
} TransProb; 

extern double getTransProb(TransProb *transprob, int order, int alpha0, int alpha1, int alpha2, int alpha3, int alpha4, int alpha5);
extern void copyTransProb(TransProb *src, TransProb *dest);
extern void resetTransCount(TransCount *transcount, const int bgPseudocount);

extern void tallyToTransCount(TransCount *transcount, int *seq, int len);
extern void marginalizeTransCount(TransCount *transcount);
extern TransProb* normalizeTransCount(TransCount *transcount);

extern void printTransCount(FILE *fptr, TransCount *transcount, bool displayAll);
extern void printTransProb(FILE *fptr, TransProb *transprob, bool displayAll);

#endif
