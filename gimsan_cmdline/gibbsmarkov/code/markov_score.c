#include "markov_score.h"


//--------------------------------------------------------------------------------------------
// Supplementary functions 
//
// Note: Do not add any function that requires headers other than dataset.h to avoid bloating!
//--------------------------------------------------------------------------------------------

double computeEntropyPerSite(Markov *markov, Dataset *data, Zoops *zoops, int *sites, int seqind, int span) {
	int count[SPAN_CAPACITY][NUMALPHAS];
	for(int i = 0; i < span; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			count[i][j] = 0;
		}
	}

	for(int n = 0; n < data->numseqs; n++) {
		if(sites[n] >= 0) {
			for(int i = 0; i < span; i++) {
				count[i][data->seqs[n][i + sites[n]]]++;
			}
		}
	}

	//int numsites = getNumSites(sites, data->numseqs);

	if(DEBUG0) {
		for(int i = 0; i < span; i++) {
			int sum = 0;
			for(int j = 0; j < NUMALPHAS; j++) {
				sum += count[i][j];
			}
			assert(sum == getNumSites(sites, data->numseqs));
		}
	}

	//make larger pwm for GAP_CHAR
	double lpwm[SPAN_CAPACITY][NUMCHARS];
	for(int m = 0; m < span; m++) {
		double sum = 0.0;
		for(int a = 0; a < NUMALPHAS; a++) {
			lpwm[m][a] = count[m][a];
			sum += lpwm[m][a];
		}
		for(int a = 0; a < NUMALPHAS; a++) {
			lpwm[m][a] /= sum;
		}

		//set gap chars to be 0.0. (But there is already a sanity check to throw exception for these)
		lpwm[m][GAP_CHAR] = 0.0; 
	}

	if(DEBUG0) {
		if(sites[seqind] < 0) {
			fprintf(stderr, "Error: invalid site in computeEntropyPerSite");
			exit(1);
		}
	}
	double cursor = 1.0;
	for(int m = 0; m < span; m++) {
		cursor *= lpwm[m][ data->seqs[seqind][sites[seqind]+m] ];
	}
	double entropy = log((double)cursor) - log((double)markov->bgscore[span][seqind][sites[seqind]]);

	if(zoops->isZoopsMode) {
		//numValidSites should never be 0 here, because there was already a check for sites[seqind]>=0
		entropy -= log((double)data->numValidSites[span][seqind]); //normalization for ZOOPS

		if(DEBUG0) {
			assert(data->numValidSites[span][seqind] > 0);
		}
	}

	return entropy;
}

static
double computeEntropy(Markov *markov, Dataset *data, Zoops *zoops, int *sites, int span, double pseudoweight) {
	if(DEBUG0) {
		for(int n = 0; n < data->numseqs; n++) {
			if(sites[n] >= 0) {
				assert(sites[n] + span - 1 < data->seqlen[n]);
			}
		}
		//ensure that all positions considered does not have GAP_CHAR (N, W, etc)
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				for( int m =0; m < span; m++) {
					if(data->seqs[i][(m+sites[i])] >= NUMALPHAS) {
						fprintf(stderr, "Error: GAP_CHAR at computeEntropy!\n");
						fprintf(stderr, "seq %d, site %d, alpha %d", i, sites[i], data->seqs[i][(m+sites[i])]);
						exit(1);
					}
				}
				if(!data->isValidSite[span][i][sites[i]]) {
					fprintf(stderr, "Error: invalid site in computeEntropy: span %d, seqind %d, sites %d\n", span, i, sites[i]);
					exit(1);
				}
			}
		}
	}

	int count[SPAN_CAPACITY][NUMALPHAS];
	for(int i = 0; i < span; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			count[i][j] = 0;
		}
	}

	for(int n = 0; n < data->numseqs; n++) {
		if(sites[n] >= 0) {
			for(int i = 0; i < span; i++) {
				count[i][data->seqs[n][i + sites[n]]]++;
			}
		}
	}

	int numsites = getNumSites(sites, data->numseqs);

	if(DEBUG0) {
		for(int i = 0; i < span; i++) {
			int sum = 0;
			for(int j = 0; j < NUMALPHAS; j++) {
				sum += count[i][j];
			}
			assert(sum == getNumSites(sites, data->numseqs));
		}
	}

	//make larger pwm for GAP_CHAR
	double lpwm[SPAN_CAPACITY][NUMCHARS];
	for(int m = 0; m < span; m++) {
		double sum = 0.0;
		for(int a = 0; a < NUMALPHAS; a++) {
			lpwm[m][a] = count[m][a] + (data->ntFreq[a] * pseudoweight * numsites);
			sum += lpwm[m][a];
		}
		for(int a = 0; a < NUMALPHAS; a++) {
			lpwm[m][a] /= sum;
		}

		//set gap chars to be 0.0. (But there is already a sanity check to throw exception for these)
		lpwm[m][GAP_CHAR] = 0.0; 
	}

	double entropy = 0.0;
	for(int i = 0; i < data->numseqs; i++) {
		if(sites[i] >= 0) {
			double cursor = 1.0;
			for(int m = 0; m < span; m++) {
				cursor *= lpwm[m][ data->seqs[i][sites[i]+m] ];
			}

			entropy += log((double)cursor) - log((double)markov->bgscore[span][i][sites[i]]);
		}
	}

	//ZOOPS term
	if(zoops->isZoopsMode) {
		entropy += getZoopsNormalizationTermLn(zoops, sites, data->numValidSites[span]);
	}

	return entropy;
}

// Generalization of "entropy score" under k-th order Markov background.
double computeEntropy(Markov *markov, Dataset *data, Zoops *zoops, int *sites, int span) {
	//double pseudocount[NUMALPHAS] = {0.0, 0.0, 0.0, 0.0};
	double pseudoweight = 0.0;
	return computeEntropy(markov, data, zoops, sites, span, pseudoweight);
}



//static 
//int getIlrNormFactor(Dataset *data, int seqind, int span) {
//	if(data->useRevcompl) {
//		return data->seqlen[seqind] - 2 * span + 2;
//	}
//	else {
//		return data->seqlen[seqind] - span + 1;
//	}
//}

//Generalization of ILR under k-th order markov background
double computeIlr(Markov *markov, Dataset *data, Zoops *zoops, int *sites, int span) {

	if(DEBUG0) {
		for(int n = 0; n < data->numseqs; n++) {
			if(sites[n] >= 0) {
				assert(sites[n] + span - 1 < data->seqlen[n]);
			}
		}
		//ensure that all positions considered does not have GAP_CHAR (N, W, etc)
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				for( int m =0; m < span; m++) {
					if(data->seqs[i][(m+sites[i])] >= NUMALPHAS) {
						fprintf(stderr, "Error: GAP_CHAR at computeIlr!\n");
						fprintf(stderr, "seq %d, site %d, alpha %d", i, sites[i], data->seqs[i][(m+sites[i])]);
						exit(1);
					}
				}
			}
		}
	}

	//create count matrix
	int count[SPAN_CAPACITY][NUMCHARS];
	for(int i = 0; i < span; i++) {
		for(int j = 0; j < NUMCHARS; j++) {
			count[i][j] = 0;
		}
	}
	for(int n = 0; n < data->numseqs; n++) {
		if(sites[n] >= 0) {
			for(int i = 0; i < span; i++) {
				count[i][data->seqs[n][i + sites[n]]]++;
			}
		}
	}

	int numsites = getNumSites(sites, data->numseqs);

	double ilrval = 0.0;
	//for(int i = 0; i < data->numseqs; i++) {
	//	if(sites[i] >= 0) {
	//		double sum = 0.0;
	//		for(int j = 0; j < data->seqlen[i] - span + 1; j++) {
	//			double cursor = data->isValidSite[span][i][j];
	//			for(int m = 0; m < span; m++) {
	//				cursor *= ((double)count[m][ data->seqs[i][j+m] ]) / numsites;
	//			}
	//			cursor /= markov->bgscore[span][i][j];
	//			sum += cursor;
	//		}
	//		ilrval += log((double)sum) - log((double)data->numValidSites[span][i]);

	//		if(DEBUG0) {
	//			assert(data->numValidSites[span][i] > 0); //since sites[i] >= 0
	//		}
	//	}
	//}
	//ZOOPS term
	//if(zoops->isZoopsMode) {
	//	ilrval += getZoopsNormalizationTermLn(zoops, sites, data->numValidSites[span]);
	//}

	for(int i = 0; i < data->numseqs; i++) {
		if(data->numValidSites[span][i] > 0 && sites[i]>= 0) {
			double sum = 0.0;
			for(int j = 0; j < data->seqlen[i] - span + 1; j++) {
				double cursor = data->isValidSite[span][i][j];
				for(int m = 0; m < span; m++) {
					cursor *= ((double)count[m][ data->seqs[i][j+m] ]) / numsites;
				}
				cursor /= markov->bgscore[span][i][j];
				sum += cursor;
			}
			ilrval += log((double)sum) - log((double)data->numValidSites[span][i]);
		}
	}

	return ilrval; 
}

//Generalization of ILR under k-th order markov background
double computeIlr(Markov *markov, Dataset *data, double** pwm, int span) {

	//verify that sum of count is the number of sequences
	if(DEBUG0) {
		for(int i = 0; i < span; i++) {
			double sum = 0.0;
			for(int j = 0; j < NUMALPHAS; j++) {
				sum += pwm[i][j];
			}
			assert(fabs(sum - 1.0) < 0.0000000001);
		}
	}

	//make larger pwm for GAP_CHAR
	double lpwm[SPAN_CAPACITY][NUMCHARS];
	for(int m = 0; m < span; m++) {
		for(int a = 0; a < NUMALPHAS; a++) {
			lpwm[m][a] = pwm[m][a];
		}
		lpwm[m][GAP_CHAR] = 0.0;
	}

	double ilrval = 0.0;
	for(int i = 0; i < data->numseqs; i++) {
		if(data->numValidSites[span][i] > 0) {
			double sum = 0.0;
			for(int j = 0; j < data->seqlen[i] - span + 1; j++) {
				double cursor = (double) data->isValidSite[span][i][j];
				for(int m = 0; m < span; m++) {
					cursor *= lpwm[m][ data->seqs[i][j+m] ];
				}
				cursor /= markov->bgscore[span][i][j];
				sum += cursor;
			}
			ilrval += log((double)sum) - log((double)data->numValidSites[span][i]);
		}
	}

	return ilrval; 
}
