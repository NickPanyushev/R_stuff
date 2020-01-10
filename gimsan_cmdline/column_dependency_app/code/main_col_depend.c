/**
*
* This file contains the main() function for ColumnDependency. 
*
**/

#include "col_depend.h"
#include "params.h"
#include "results.h"

ColDepend* constructColDepend() {
	ColDepend *cd = (ColDepend*) malloc(sizeof(ColDepend));
	return cd;
}


void destructColDepend(ColDepend *coldep) {
	for(int i = 0; i < coldep->numPairs; i++) {
		destructColPair(coldep->pairs[i]);
	}
	free(coldep->pairs);
	destructIntMat(coldep->kmerset);
	destructIntMat(coldep->countmat);
	destructDblMat(coldep->freqmat);
	free(coldep);
}
//----------------------------------------------------------------------
// Main
//----------------------------------------------------------------------
int main(int argc, char *argv[]) {
	if(DEBUG0) {
		printf("WARNING: This is currently running under DEBUG mode.\n");
	}
	if(DEBUG1) {
		printf("Verbose mode.\n");
	}
	printf("Compiled on " __DATE__ " " __TIME__ "\n");
	printf("\n");

	printf("ChangeLog\n");
	printf("20080911. Added StDev mode \n");
	printf("20080608. Output with summary of pairs with significant dependency \n");
	printf("\n");

	if(argc <= 1) {
		printUsage();
	}

	ColDepend *coldep = constructColDepend();
	Params *params = constructParams(argc, argv);
	setColDependParams(params, coldep);

	iterColumnPairs(coldep);

	Results *results = constructAndSetResults(coldep->kmerset->numrows);
	displayColumnPairsContingencyTable(params->fptr, coldep, results);
	displayColumnPairsResults(params->fptr, coldep);
	displayKmerSet(params->fptr, coldep->kmerset);

	destructParams(params);
	destructResults(results);
	destructColDepend(coldep);

	return(0);
}

