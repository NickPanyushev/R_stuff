#ifndef _ZOOPS_H
#define _ZOOPS_H

#include "stdinc.h"

typedef struct {
	bool isZoopsMode;
	double *zeroOccurrenceSamplingScore; //only defined for 1 <= numsites < numseqs 
	int numseqs; //total number of sequences (or maximum number of sites)
	double zoopsPseudoweight;
} Zoops; //lookup table for entropy sampling

//precompute regardless isZoops is true or false
extern Zoops* initZoops(bool isZoopsMode, int numseqs, double zoopsPseudoweight);

//natural log of the Beta function term
//only defined for 1 <= numsites <= numseqs
extern double getZoopsNormalizationTermLn(Zoops *zoops, int *sites, int *numValidSites);

extern int getNumSites(int *sites, int numseqs);

extern void nilZoops(Zoops *zoops);

#endif
