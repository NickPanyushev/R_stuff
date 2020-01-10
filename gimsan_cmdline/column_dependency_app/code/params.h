#ifndef _PARAMS_H
#define _PARAMS_H

#include "col_depend.h"
#include "dataset.h"
#include "siteloc.h"

typedef struct {
	char *fastaFilename;
	Dataset *fasta;

	char *sitelocFilename;
	Siteloc *siteloc;

	unsigned int randomSeed;

	enum StopMode stopMode;

	FILE *fptr;

	//ColDepend params
	int initNumPerms;
	int capacityNumPerms; //maximum number of possible permutations allowed
	int motifspan;
	double alpha;
	double beta;
	double lowerColInfoContent;
	double upperColInfoContent;
} Params; 

extern void printUsage();
extern void setColDependParams(Params *params, ColDepend *coldep);
extern Params* constructParams(int argc, char *argv[]);
extern void destructParams(Params *params);

#endif

