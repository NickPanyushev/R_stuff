#include "col_depend.h"
#include "binomial_distribution.h"
//#include <boost/math/distributions/binomial.hpp>



double computeDnaInfoContent(double *arr) {
	//arr must be of size NUMALPHAS
	double sum = 0.0;
	for(int i = 0; i < NUMALPHAS; i++) {
		sum += arr[i];
	}

	double ic = log((double)NUMALPHAS);
	for(int i = 0; i < NUMALPHAS; i++) {
		double p = arr[i] / sum;
		if(p > 1e-10) { //avoid log(0.0)
			ic += p * log((double)p);
		}
	}
	return (ic / log(2.0)); //covert base-2 log
}

int determineGoodPos(DblMat *freqmat, double lowerBound, double upperBound, bool isGoodPos[SPAN_CAPACITY]) {
	//reset isGoodColumn to false
	for(int i = 0; i < SPAN_CAPACITY; i++) {
		isGoodPos[i] = false;
	}

	int count = 0;
	for(int i = 0; i < freqmat->numrows; i++) { //numrows is the span of the count-matrix
		double ic = computeDnaInfoContent(freqmat->vals[i]); //vals[i] is a vector of length 4
		if(ic >= lowerBound && ic <= upperBound) {
			isGoodPos[i] = true;
			count++;
		}
	}
	return count;
}

void constructColumnPairs(ColDepend *coldep, bool isGoodPos[], int numGoodPos) {
	coldep->numAnalyzedCols = numGoodPos;
	coldep->numPairs = (numGoodPos * (numGoodPos-1)) / 2; //nchoosek(numGoodPos,2)
	coldep->pairs = (ColPair**) malloc(sizeof(ColPair*) * coldep->numPairs);

	coldep->alphaBonferroni = coldep->alpha / coldep->numPairs;//Bonferroni-corrected alpha cutoff

	int count = 0;
	for(int i = 0; i < coldep->motifspan; i++) {
		if(isGoodPos[i]) {
			for(int j = i+1; j < coldep->motifspan; j++) {
				if(isGoodPos[i] && isGoodPos[j]) {
					coldep->pairs[count] = (ColPair*) malloc(sizeof(ColPair));
					coldep->pairs[count]->colIndex1 = i;
					coldep->pairs[count]->colIndex2 = j;

					coldep->pairs[count]->vectLen = coldep->kmerset->numrows;
					coldep->pairs[count]->colVect1 = (int*) calloc(coldep->pairs[count]->vectLen, sizeof(int));
					coldep->pairs[count]->colVect2 = (int*) calloc(coldep->pairs[count]->vectLen, sizeof(int));
					for(int k = 0; k < coldep->pairs[count]->vectLen; k++) {
						coldep->pairs[count]->colVect1[k] = coldep->kmerset->vals[k][i];
						coldep->pairs[count]->colVect2[k] = coldep->kmerset->vals[k][j];
					}

					count++;
				}
			}
		}
	}

	if(DEBUG0) {
		assert(count == coldep->numPairs);
	}
	
}

//vect1 and vect2 with length that represents the number of sites
//Its values are in {0,1,2,3} that represents the alphabet of nucleotides
double computeRelativeEntropy(int *vect1, int *vect2, int len) {
	double freqvect1[NUMALPHAS];
	double freqvect2[NUMALPHAS];
	double freqmat[NUMALPHAS][NUMALPHAS];

	//set H0 values
	for(int a = 0; a < NUMALPHAS; a++) {
		freqvect1[a] = 0.0;
		freqvect2[a] = 0.0;
	}

	for(int i = 0; i < len; i++) {
		freqvect1[vect1[i]] += 1.0;
		freqvect2[vect2[i]] += 1.0;
	}
	for(int i = 0; i < NUMALPHAS; i++) {
		freqvect1[i] /= len;
		freqvect2[i] /= len;
	}

	//set H1 values
	for(int i = 0; i < NUMALPHAS; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			freqmat[i][j] = 0.0;
		}
	}
	for(int i = 0; i < len; i++) {
		freqmat[vect1[i]][vect2[i]] += 1.0;
	}
	for(int i = 0; i < NUMALPHAS; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			freqmat[i][j] /= len;
		}
	}

	//compute relative entropy
	double ent = 0.0;
	for(int i = 0; i < NUMALPHAS; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			if(freqmat[i][j] > 1e-10) { 
				//avoid log(0.0). If freqmat > 0, then vect1 > 0 and vect2 > 0.
				ent += freqmat[i][j] * ( log(freqmat[i][j])-log(freqvect1[i])-log(freqvect2[j]) );
			}
		}
	}

	if(DEBUG0) {
		double sum;
		sum = 0.0;
		for(int i = 0; i < NUMALPHAS; i++) {
			sum += freqvect1[i];
		}
		assert(fabs(sum - 1.0) < 0.0000001);
		sum = 0.0;
		for(int i = 0; i < NUMALPHAS; i++) {
			sum += freqvect2[i];
		}
		assert(fabs(sum - 1.0) < 0.0000001);
		sum = 0.0;
		for(int i = 0; i < NUMALPHAS; i++) {
			for(int j = 0; j < NUMALPHAS; j++) {
				sum += freqmat[i][j];
			}
		}
		assert(fabs(sum - 1.0) < 0.0000001);
	}
	return ent / log(2.0); // convert to base-2
}


//-------------------------------------------------------------------------------
// Vector functions
//-------------------------------------------------------------------------------
int* copyIntVect(int *intVect, int len ) {
	int *newvect = (int*) malloc(sizeof(int) * len);
	for(int i = 0; i < len; i++) {
		newvect[i] = intVect[i];
	}
	return newvect;
}

void permuteIntVect(int *intVect, int len) {
	//if(DEBUG1) {
	//	fprintf(stderr, "original vector:\n");
	//	for(int i = 0; i < len; i++) {
	//		fprintf(stderr, "%d ", intVect[i]);
	//	}
	//	fprintf(stderr, "\n");
	//}

	for(int i = 0; i < len; i++) {
		int rint = ((int)(Random()*(len - i))) + i;

		//swap element i and rint
		if(i != rint) {
			int temp = intVect[i];
			intVect[i] = intVect[rint];
			intVect[rint] = temp;
		}
	}
	//if(DEBUG1) {
	//	fprintf(stderr, "permuted vector:\n");
	//	for(int i = 0; i < len; i++) {
	//		fprintf(stderr, "%d ", intVect[i]);
	//	}
	//	fprintf(stderr, "\n");
	//}
}

//-------------------------------------------------------------------------------
// Random Permutation functions
//-------------------------------------------------------------------------------
//Count the number of random permutation greater than or equal to the observed entropy
void countRndPermGEObsEnt(double obsEnt, int *vect1, int *vect2, int vectLen, int numperms, int* ret_countGe, int* ret_countEq) {
	int countGe = 0;
	int countEq = 0;
	for(int i = 0; i < numperms; i++) {
		permuteIntVect(vect2, vectLen);
		double ent = computeRelativeEntropy(vect1, vect2, vectLen);
		
		if(ent >= obsEnt - DBL_EPSILON) {
			countGe++;
		}
		if(fabs(ent - obsEnt) < DBL_EPSILON) {
			countEq++;
		}
		
		if(DEBUG1) {
			//fprintf(stderr, "rel-ent = %0.6lf\n", ent);
		}
	}
	*ret_countGe = countGe;
	*ret_countEq = countEq;
}


void randomPermutationAndPerformMultiTest(ColPair *colpair, ColDepend *coldep) {
	colpair->obsEnt = computeRelativeEntropy(colpair->colVect1, colpair->colVect2, colpair->vectLen);

	int *vect1 = copyIntVect(colpair->colVect1, colpair->vectLen);
	int *vect2 = copyIntVect(colpair->colVect2, colpair->vectLen);
	int vectLen = colpair->vectLen;

	int numperms = coldep->initNumPerms;
	double bfnCutoff = coldep->alphaBonferroni;

	colpair->numPerms = 0;
	colpair->numObsEntGe = 0;
	colpair->numObsEntEq = 0;
	colpair->pval = NAN;
	colpair->pvalLower = NAN;
	colpair->pvalUpper = NAN;

	if(DEBUG0) {
		for(int i = 0; i < colpair->vectLen; i++) {
			if(vect1[i] < 0 || vect1[i] >= NUMALPHAS || vect2[i] < 0 || vect2[i] >= NUMALPHAS) {
				fprintf(stderr, "Error: value of vect1 and vect2 are not in range.");
				abort();
			}
		}
	}
	if(DEBUG1) {
		fprintf(stderr, "========================================================\n");
		fprintf(stderr, "Column pair (%d,%d)\n", colpair->colIndex1, colpair->colIndex2);
	}


	while(numperms <= coldep->capacityNumPerms) {
		int countGe = INT_MIN;
		int countEq = INT_MIN;
		countRndPermGEObsEnt(colpair->obsEnt, vect1, vect2, vectLen, numperms, &countGe, &countEq);

		colpair->numObsEntGe += countGe;
		colpair->numObsEntEq += countEq;
		colpair->numPerms += numperms;

		numperms = colpair->numPerms;

		//fprintf(stderr, "(%d) ", colpair->numPerms);
		//fprintf(stderr, "%d ", numperms);

		//computing early stopping condition
		colpair->pval = colpair->numObsEntGe / (double)numperms; //point estimate (MLE)
		//colpair->pvalLower = boost::math::binomial_distribution<>::find_lower_bound_on_p(numperms, colpair->numObsEntGe, (1.0-coldep->beta)/2);
		//colpair->pvalUpper = boost::math::binomial_distribution<>::find_upper_bound_on_p(numperms, colpair->numObsEntGe, (1.0-coldep->beta)/2);
		colpair->pvalLower = binom_distrib_find_lower_bound_on_p(numperms, colpair->numObsEntGe, (1.0-coldep->beta)/2);
		colpair->pvalUpper = binom_distrib_find_upper_bound_on_p(numperms, colpair->numObsEntGe, (1.0-coldep->beta)/2);
			

		fprintf(stderr, "trials=%d; successes=%d; alpha=%.5lg; lowp=%.5lg; highp=%.5lg; bfnCutoff=%.5lg\n",
			numperms,
			colpair->numObsEntGe, 
			(1.0-coldep->beta)/2,
			colpair->pvalLower,
			colpair->pvalUpper,
			bfnCutoff
		);

		if(coldep->stopMode == ALPHA_MODE) {
			if(colpair->pvalLower > bfnCutoff || bfnCutoff > colpair->pvalUpper) {
				break;
			}
		}
		else if(coldep->stopMode == STDEV_MODE) {
			double pMidP = (colpair->numObsEntGe - colpair->numObsEntEq/2.0) / (double)numperms;
			double estPeq = colpair->numObsEntEq / (double)numperms;
			double estPgt = (colpair->numObsEntGe - colpair->numObsEntEq) / (double)numperms;
			double estStd = sqrt( estPeq*(1-estPeq)/(2*numperms) + estPgt*(1-estPgt)/numperms - estPeq*estPgt/numperms );
			if (estStd*3 <= pMidP/2 && colpair->numObsEntGe > 0)  {
				break;
			}
		}
		else {
			fprintf(stderr, "Invalid stop mode\n");
			abort();
		}
	}

	free(vect1);
	free(vect2);
}

//-------------------------------------------------------------------------------
// Core loop
//-------------------------------------------------------------------------------
void iterColumnPairs(ColDepend *coldep) {
	bool isGoodPos[SPAN_CAPACITY];
	int numGoodPos = determineGoodPos(coldep->freqmat, coldep->lowerColInfoContent, coldep->upperColInfoContent, isGoodPos);

	constructColumnPairs(coldep, isGoodPos, numGoodPos);

	for(int i = 0; i < coldep->numPairs; i++) {
		randomPermutationAndPerformMultiTest(coldep->pairs[i], coldep);

		if(DEBUG1) {
			fprintf(stderr, "Index %d, numPerms %d\n", i, coldep->pairs[i]->numPerms);
		}
	}
}

void destructColPair(ColPair *pair) {
	free(pair->colVect1);
	free(pair->colVect2);
	free(pair);
}

