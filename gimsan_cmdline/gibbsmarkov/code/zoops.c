#include "zoops.h"

static 
double betaln(double x, double y) {
	if(DEBUG0) {
		assert(x > 0.0 && y > 0.0);
	}
	return (lgamma(x) + lgamma(y) - lgamma(x+y));
}

static
double betalnTerm(Zoops *zoops, int numsites) {
	double pseudocount = zoops->zoopsPseudoweight * zoops->numseqs;
	return betaln( numsites + pseudocount, zoops->numseqs - numsites + pseudocount);
}

static
void setZeroOccurrenceSamplingScore(Zoops *zoops) {
	zoops->zeroOccurrenceSamplingScore[0] = NAN;

	for(int i = 1; i < zoops->numseqs; i++) {
		zoops->zeroOccurrenceSamplingScore[i] = exp(betalnTerm(zoops, i) - betalnTerm(zoops, i+1));
	}
	
}

//precompute regardless isZoops is true or false
Zoops* initZoops(bool isZoopsMode, int numseqs, double zoopsPseudoweight) {
	Zoops *zoops = (Zoops*) malloc(sizeof(Zoops));
	zoops->isZoopsMode = isZoopsMode;
	zoops->numseqs = numseqs;
	zoops->zoopsPseudoweight = zoopsPseudoweight;
	//zoops->normalizationTermLn = (double*) calloc(numseqs+1, sizeof(double));
	zoops->zeroOccurrenceSamplingScore = (double*) calloc(numseqs, sizeof(double));

	if(isZoopsMode) {
		setZeroOccurrenceSamplingScore(zoops);
	}
	else {
		for(int i = 0; i < zoops->numseqs; i++) {
			zoops->zeroOccurrenceSamplingScore[i] = NAN;
		}
	}

	return zoops;
}


int getNumSites(int *sites, int numseqs) {
	int count = 0;

	for(int i = 0; i < numseqs; i++) {
		if(sites[i]>= 0) {
			count++;
		}
	}

	if(DEBUG0) {
		if( count == 0 ) {
			fprintf(stderr, "Error: numsites is 0\n");
			exit(1);
		}
	}

	return count;
}

//natural log of the Beta function term
double getZoopsNormalizationTermLn(Zoops *zoops, int *sites, int *numValidSites) {
	int numsites = getNumSites(sites, zoops->numseqs);

	if(DEBUG0) {
		assert(1 <= numsites && numsites <= zoops->numseqs);
		assert(zoops->isZoopsMode);
	}

	double lnTerm = betalnTerm(zoops, numsites) - betalnTerm(zoops, 0);
	for(int i = 0; i < zoops->numseqs; i++) {
		if(sites[i] >= 0) {
			//numValidSites should never be 0 here
			lnTerm -= log((double)numValidSites[i]);

			if(DEBUG0) {
				assert(numValidSites[i] > 0) ;
			}
		}
	}
	return lnTerm;
}

void nilZoops(Zoops *zoops) {
	free(zoops->zeroOccurrenceSamplingScore);
	free(zoops);
}

