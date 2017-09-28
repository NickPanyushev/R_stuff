#ifndef _RESULTS_H
#define _RESULTS_H

#include "col_depend.h"

typedef struct{
	double *log_fact; //precomputed log-factorial values
	int fact_capacity; 

	double hypgeoPvalThreshold; //default: 0.05 (threshold to display the p-values in the contingency tables)
} Results;

extern Results* constructAndSetResults(int numsites);
extern void displayKmerSet(FILE *fptr, IntMat *kmerset);
extern void displayColumnPairsResults(FILE *fptr, ColDepend *coldep);
extern void displayColumnPairsContingencyTable(FILE *fptr, ColDepend *coldep, Results *results);
extern void destructResults(Results *res);

#endif
