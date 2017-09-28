#ifndef _COL_DEPEND_H
#define _COL_DEPEND_H

#include "stdinc.h"
#include "random.h"
#include "matrix.h"

enum StopMode {ALPHA_MODE, STDEV_MODE};

typedef struct {
	int colIndex1;
	int colIndex2;
	int *colVect1;
	int *colVect2;
	int vectLen;

	//results
	double pval; //point estimate (MLE)
	double pvalLower; //lower bound of confidence interval
	double pvalUpper;

	double obsEnt; //observed entropy of the column-pairs
	int numObsEntGe; //number of random permutations greater than or equal to observed entropy 
	int numObsEntEq; //number of random permutations equal to observed entropy
	int numPerms; //total number of random permutations
} ColPair;

typedef struct {
	int initNumPerms;
	int capacityNumPerms; //maximum number of possible permutations allowed
	
	ColPair **pairs;
	int numPairs;
	int numAnalyzedCols;
	int motifspan;

	double alpha;
	double beta;
	double alphaBonferroni; //Bonferroni-corrected alpha-level;

	double lowerColInfoContent;
	double upperColInfoContent;

	enum StopMode stopMode;

	IntMat *kmerset; //a set of kmers
	IntMat *countmat;
	DblMat *freqmat;
} ColDepend; 

extern void iterColumnPairs(ColDepend *coldep);
extern void destructColPair(ColPair *pair);

#endif
