#include "results.h"
#include "symbols.h"

#define BUFFER_SIZE 100

Results* constructAndSetResults(int numsites){
	Results *res = (Results*) malloc(sizeof(Results));
	res->hypgeoPvalThreshold = 0.05; 
	res->fact_capacity = numsites + 10; //set a bit of extra space
	res->log_fact = (double*) malloc(sizeof(double) * (res->fact_capacity+1));

	res->log_fact[0] = 0.0;
	for(int i = 1; i < res->fact_capacity + 1; i++) {
		res->log_fact[i] = res->log_fact[i-1] + log((double)i);
	}
	return res;
}

void destructResults(Results *res) {
	free(res->log_fact);
	free(res);
}

void displayKmerSet(FILE *fptr, IntMat *kmerset) {
	fprintf(fptr, "Input set of %d k-mers (k=%d):\n", kmerset->numrows, kmerset->numcols);
	for(int i = 0; i < kmerset->numrows; i++) {
		for(int j = 0; j < kmerset->numcols; j++) {
			fprintf(fptr, "%c", numToChar(kmerset->vals[i][j]));
		}
		fprintf(fptr, "\n");
	}
	fprintf(fptr, "\n");
}

void displayColumnPairsResults(FILE *fptr, ColDepend *coldep) {
	//compute lowest point-estimate
	double lowestPval = DBL_MAX;
	for(int i = 0; i < coldep->numPairs; i++) {
		if(lowestPval > coldep->pairs[i]->pval) {
			lowestPval = coldep->pairs[i]->pval;
		}
	}
	fprintf(fptr, "==========================================================================================================\n");
	fprintf(fptr, "\n");
	fprintf(fptr, "Number of columns analyzed: %d (out of %d)\n", coldep->numAnalyzedCols, coldep->motifspan);
	fprintf(fptr, "Number of column-pairs: %d\n", coldep->numPairs);
	fprintf(fptr, "Bonferroni corrected alpha level: %.5lg\n", coldep->alphaBonferroni);
	fprintf(fptr, "Lowest point-estimated p-value: %.2le\n", lowestPval);
	fprintf(fptr, "\n");
	fprintf(fptr, "Column dependency p-values:\n");
	for(int i = 0; i < coldep->numPairs; i++) {
		ColPair *pair = coldep->pairs[i];
		//column-index starts at 1 instead of 0
		fprintf(fptr, "Column-pair (%2d,%2d) has estimated p-value %0.6lf [%0.6lf,%0.6lf] ", 
			pair->colIndex1 + 1, pair->colIndex2 + 1, pair->pval, pair->pvalLower, pair->pvalUpper);
		fprintf(fptr, "with ent=%5.3lf bits and nPerms=%d\n", pair->obsEnt, pair->numPerms);
	}
	fprintf(fptr, "\n");
}
//=========================================================================================
// Contingency table
//=========================================================================================

static
void computeObsCount(ColPair *pair, IntMat *obsCount) {
	for(int i = 0; i < NUMALPHAS; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			obsCount->vals[i][j] = 0;
		}
	}
	for(int i = 0; i < pair->vectLen; i++) {
		obsCount->vals[pair->colVect1[i]][pair->colVect2[i]]++;
	}
}

static
int round_pos_dbl(double d) {
	if(d < 0.0){
		fprintf(stderr, "Error: only positive numbers for rounding\n");
		abort();
	}
	else {
		return (int)(d + 0.5 + DBL_EPSILON);
	}
}

static
void computeExpCount(ColPair *pair, IntMat *expCount) {
	int numsites = pair->vectLen;
	for(int i = 0; i < NUMALPHAS; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			expCount->vals[i][j] = 0;
		}
	}
	int countvect1[NUMALPHAS];
	int countvect2[NUMALPHAS];
	for(int i = 0; i < NUMALPHAS; i++) {
		countvect1[i] = 0;
		countvect2[i] = 0;
	}
	for(int i = 0; i < pair->vectLen; i++) {
		countvect1[pair->colVect1[i]]++;
		countvect2[pair->colVect2[i]]++;
	}
	for(int i = 0; i < NUMALPHAS; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			expCount->vals[i][j] = round_pos_dbl(((double)countvect1[i]) * countvect2[j] / ((double)numsites));
		}
	}
}

static
double log_nchoosek(int n, int k, double *log_fact) {
	return log_fact[n] - log_fact[k] - log_fact[n-k];
}

static
double log_sum(double log_a, double log_b) {
	//compute log(a+b) from log(a) and log(b)
	return ((log_a > log_b) ?
		log_a + log(1+exp(log_b-log_a)) :
		log_b + log(1+exp(log_a-log_b)));
}

inline 
int max_int(int a, int b) {
	return (a > b ? a : b);

}
inline
int min_int(int a, int b) {
	return (a < b ? a : b);
}

static
void computeHypgeoPval(ColPair *pair, IntMat *obsCount, DblMat *hypgeoPval, Results *results) {
	//compute hypergeometric probability for the tail with less than or equal to 0.5

	//compute countvect
	int countvect1[NUMALPHAS];
	int countvect2[NUMALPHAS];
	for(int i = 0; i < NUMALPHAS; i++) {
		countvect1[i] = 0;
		countvect2[i] = 0;
	}
	for(int i = 0; i < pair->vectLen; i++) {
		countvect1[pair->colVect1[i]]++;
		countvect2[pair->colVect2[i]]++;
	}

	//compute p-values
	if(DEBUG1) {
		fprintf(stderr, "Hypergeometric p-values\n\n");
	}

	int numsites = pair->vectLen;
	for(int i = 0; i < NUMALPHAS; i++) {
		//precomputation
		double log_denominator = log_nchoosek(numsites, countvect1[i], results->log_fact);

		//iterate through each cell
		for(int j = 0; j < NUMALPHAS; j++) {

			//set limits and bounds
			int k_lower_bound = max_int(countvect1[i] + countvect2[j] - numsites, 0);
			int k_upper_bound = min_int(countvect1[i], countvect2[j]); 

			int leftind, rightind;
			//compare with the empirical expected value of countvect1[i]*countvect2[j]/(double)numsites
			double expectation = countvect1[i]*countvect2[j] / ((double)numsites);
			if(obsCount->vals[i][j] > expectation) {
				leftind = max_int(k_lower_bound, obsCount->vals[i][j]); 
				rightind = k_upper_bound;
			}
			else if(obsCount->vals[i][j] < expectation) {
				leftind = k_lower_bound;
				rightind = min_int(k_upper_bound, obsCount->vals[i][j]);
			}
			else {
				hypgeoPval->vals[i][j] = 0.5;
				continue;
			}

			//summing over log-probablity
			double log_prob = NINF; //this doesn't work in Windows for some reason
			for(int k = leftind; k <= rightind; k++) { //rightind is inclusive
				double log_numerator = log_nchoosek(countvect2[j], k, results->log_fact) 
					+ log_nchoosek(numsites - countvect2[j], countvect1[i]- k, results->log_fact);
				if(k == leftind) {
					log_prob = log_numerator - log_denominator;
				}
				else {
					log_prob = log_sum(log_prob, log_numerator - log_denominator);
				}
			}
			hypgeoPval->vals[i][j] = exp(log_prob);

			//finding the tail with less than 0.5
			//the expectation does not determine the tail because the distribution is not symmetric
			if(hypgeoPval->vals[i][j] > 0.5) {
				hypgeoPval->vals[i][j] = 1.0 - hypgeoPval->vals[i][j];
			}

			if(DEBUG0) {
				double tmp_log_prob = NINF;
				for(int k = k_lower_bound; k <= k_upper_bound; k++) {
					double log_numerator = log_nchoosek(countvect2[j], k, results->log_fact) 
						+ log_nchoosek(numsites - countvect2[j], countvect1[i]- k, results->log_fact);
					if(k == k_lower_bound) {
						tmp_log_prob = log_numerator - log_denominator;
					}
					else {
						tmp_log_prob = log_sum(tmp_log_prob, log_numerator - log_denominator);
					}
				}
				if(fabs(tmp_log_prob - log(1.0)) > 0.000000001) {
					fprintf(stderr, "Error: the sum of hypergeometric prob is not 1.0\n");
					fprintf(stderr, "%.4lg\n", exp(tmp_log_prob));
					abort();
				}
			}
			if(DEBUG1) {
				fprintf(stderr, "alphabets: %c%c\n", numToChar(i), numToChar(j));
				fprintf(stderr, "observed count: %d\n", obsCount->vals[i][j]);
				fprintf(stderr, "expectation: %lg\n", expectation);
				fprintf(stderr, "k_lower_bound: %d\n", k_lower_bound);
				fprintf(stderr, "k_upper_bound: %d\n", k_upper_bound);
				fprintf(stderr, "leftind: %d\n", leftind);
				fprintf(stderr, "rightind: %d\n", rightind);
				fprintf(stderr, "\n");
			}
		}
	}
}


void displayColumnPairsContingencyTable(FILE *fptr, ColDepend *coldep, Results *results) {
	int count = 0; 
	for(int p = 0; p < coldep->numPairs; p++) {
		ColPair *pair = coldep->pairs[p];
		if(pair->pvalUpper < coldep->alphaBonferroni) {
			count++;
		}
	}

	fprintf(fptr, "==========================================================================================================\n");
	fprintf(fptr, "\n");
	fprintf(fptr, "Hypergeometric p-values for statistically significant pairs (%d pairs)\n", count);
	fprintf(fptr, "\n");

	IntMat *obsCount = constructIntMat(NUMALPHAS, NUMALPHAS);
	IntMat *expCount = constructIntMat(NUMALPHAS, NUMALPHAS);
	DblMat *hypgeoPval = constructDblMat(NUMALPHAS, NUMALPHAS);
	
	char buffer[BUFFER_SIZE];

	for(int p = 0; p < coldep->numPairs; p++) {
		ColPair *pair = coldep->pairs[p];
		if(pair->pvalUpper < coldep->alphaBonferroni) {
			computeObsCount(pair, obsCount);
			computeExpCount(pair, expCount);
			computeHypgeoPval(pair, obsCount, hypgeoPval, results);//numrows is the number of sites

			//print overall p-value
			fprintf(fptr, "Column-pair (%2d,%2d) has estimated p-value %.2le [%.2le, %.2le] ", 
				pair->colIndex1 + 1, pair->colIndex2 + 1, pair->pval, pair->pvalLower, pair->pvalUpper);
			fprintf(fptr, "with ent=%5.3lf bits and nPerms=%d\n", pair->obsEnt, pair->numPerms);

			//print header for details
			fprintf(fptr, "%-3s %7s %7s %7s %7s %4s %10s %10s %10s %10s %-s\n",
				"O/E", "A   ", "C   ", "G   ", "T   ", " || ", "A   ", "C   ", "G   ", "T   ", "hypergeometric p-value");

			for(int i = 0; i < NUMALPHAS; i++) {
				//display observe/expected table
				fprintf(fptr, "%-3c ", numToChar(i));
				for(int j = 0; j < NUMALPHAS; j++) {
					int ret = _snprintf(buffer, BUFFER_SIZE, "%d/%-3d", obsCount->vals[i][j], expCount->vals[i][j]);
					if(ret == -1) {
						fprintf(stderr, "Error writing to buffer in displayColumnPairsContingencyTable()\n");
						abort();
					}
					fprintf(fptr, "%7s ", buffer);
				}
				fprintf(fptr, "%4s ", " || ");

				//display hypergeometric p-values table
				for(int j = 0; j < NUMALPHAS; j++) {
					double hgPval = hypgeoPval->vals[i][j] * 2; //two-sided
					char sign; 
					if(obsCount->vals[i][j] > expCount->vals[i][j]) {
						sign = '+';
					}
					else if(obsCount->vals[i][j] < expCount->vals[i][j]) {
						sign = '-';
					}
					else {
						sign = ' ';
					}

					if(hgPval < results->hypgeoPvalThreshold / (NUMALPHAS*NUMALPHAS)) {
						int ret = _snprintf(buffer, BUFFER_SIZE, "%.1le%c", hgPval, sign); 
						if(ret == -1) {
							fprintf(stderr, "Error writing to buffer in displayColumnPairsContingencyTable()\n");
							abort();
						}
						fprintf(fptr, "%10s ", buffer);
					}
					else {
						fprintf(fptr, "%10s ", "");
					}
					if(DEBUG1) {
						fprintf(stderr, "Column-pair (%2d,%2d), ", pair->colIndex1 + 1, pair->colIndex2 + 1);
						fprintf(stderr, "%c%c pval: %.5lg%c\n", numToChar(i), numToChar(j), hgPval, sign); 
					}
				}
				fprintf(fptr, "%c\n", numToChar(i));
			}
			fprintf(fptr, "\n");
		}
	}
	
	destructIntMat(obsCount);
	destructIntMat(expCount);
	destructDblMat(hypgeoPval);
}
