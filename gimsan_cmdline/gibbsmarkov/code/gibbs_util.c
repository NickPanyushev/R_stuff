#include "gibbs_util.h"


//----------------------------------------------------------------------
// RunNode/RunSet
//----------------------------------------------------------------------
RunNode* createRunNode(int runId, Dataset *data, int initspan, int maxspan) {
	RunNode *rnode = (RunNode*) malloc(sizeof(RunNode));
	rnode->countmat = initProfile(initspan, maxspan);
	rnode->pswm = initProfile(initspan, maxspan);
	rnode->next = NULL;
	rnode->runId = runId;
	rnode->sites = (int*) malloc( data->numseqs * sizeof(int));
	return rnode;
}

void nilRunNode(RunNode *rnode) {
	free(rnode->sites);
	nilProfile(rnode->countmat);
	nilProfile(rnode->pswm);
	free(rnode);
}

RunSet* createRunSet() {
	RunSet *rset = (RunSet*) malloc(sizeof(RunSet));
	rset->head = NULL;
	rset->len = 0;
	return rset;
}

void nilRunSet(RunSet *rset) {
	RunNode* garbage;
	while(rset->head != NULL) {
		garbage = rset->head;
		rset->head = rset->head->next;
		nilRunNode(garbage);
	}
	free(rset);
}


//----------------------------------------------------------------------
// Scoring functions
//----------------------------------------------------------------------

char* scoreMetricToStr(enum ScoreMetric metric) {
	if(metric == ILR) {
		return "ILR";
	}
	else if(metric == CLR){
		return "CLR";
	}
	else {
		return "ERROR";
	}
}

//----------------------------------------------------------------------
// Matrix/Profile updates
//----------------------------------------------------------------------

//Warning: Not designed for profiles with gaps
void updateCountmatFromSites(Profile *countmat, int *sites, Dataset *data) {
	if(DEBUG0) {
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				if(sites[i] + countmat->span - 1 >= data->seqlen[i]) {
					fprintf(stderr, "Error: sites out of range at updateCountMat.\n");
					fprintf(stderr, "sites[%d] %d\n", i, sites[i]);
					exit(1);
				}
			}
		}

		//ensure that all positions considered does not have GAP_CHAR (N, W, etc)
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				for(int m =0; m < countmat->span; m++) {
					assert(data->seqs[i][(m+sites[i])] < NUMALPHAS);
				}
			}
		}
	}

	for(int m =0; m < countmat->span; m++) {
		for(int a = 0; a < NUMALPHAS; a++) {
			countmat->mat[m][a] = 0.0;
		}
	}
	for(int i = 0; i < data->numseqs; i++) {
		if(sites[i] >= 0) { //non-empty sites in ZOOPS model
			for(int m =0; m < countmat->span; m++) {
				countmat->mat[m][ ((int)data->seqs[i][(m+sites[i])]) ] += 1.0;
			}
		}
	}
	if(DEBUG0) {
		for(int m = 0; m < countmat->span; m++) {
			double sum = 0.0;
			for(int a = 0; a < NUMALPHAS ; a++) {
				sum += countmat->mat[m][a];
			}
			if(fabs(sum - getNumSites(sites, data->numseqs)) > 0.000001) {
				fprintf(stderr, "Error: countmat inconsistent at updateCountmat().\n");
				exit(1);
			}
		}

		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				if(sites[i] + countmat->span - 1 >= data->seqlen[i]) {
					fprintf(stderr, "Error: sites out of range at updateCountMat.\n");
					exit(1);
				}
			}
		}
	}

}


//Warning: Not designed for profiles with gaps
boolean validCountmatWithSites(Profile *countmat, int *sites, Dataset *data) {
	int cmat[SPAN_CAPACITY][NUMALPHAS];

	if(DEBUG0) {
		//ensure that all positions considered does not have GAP_CHAR (N, W, etc)
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				for(int m =0; m < countmat->span; m++) {
					assert(data->seqs[i][(m+sites[i])] < NUMALPHAS);
				}
			}
		}
	}

	for(int  m =0; m < countmat->span; m++) {
		for(int a = 0; a < NUMALPHAS; a++) {
			cmat[m][a] = 0;
		}
	}
	for(int i = 0; i < data->numseqs; i++) {
		if(sites[i] >= 0) {
			for(int  m =0; m < countmat->span; m++) {
				cmat[m][ ((int)data->seqs[i][(m+sites[i])]) ] ++;
			}
		}
	}

	for(int m =0; m < countmat->span; m++) {
		for(int a = 0; a < NUMALPHAS; a++) {
			if(fabs((double)cmat[m][a] - countmat->mat[m][a]) > 0.0001) {
				return FALSE;
			}
		}
	}
	return TRUE;
}

//----------------------------------------------------------------------
// add/remove/set sites
//----------------------------------------------------------------------

void setRandomSites(int *sites, int initspan, Dataset *data, Zoops *zoops) {
	int numsites = 0;
	for (int i = 0; i < data->numseqs; i++) {
		bool hasSite = false;

		if(!zoops->isZoopsMode) {
			//OOPS
			hasSite = true;
		}
		else {
			//ZOOPS
			if( data->numValidSites[initspan][i] == 0) {
				hasSite = false;
			}
			else if( numsites == 0) {
				//force this motif to have at least one site
				hasSite = true; 
			}
			else if(Random() < 0.50) {
				hasSite = true;
			}
			else {
				hasSite = false;
			}
		}

		if(hasSite) {
			numsites++;
			int range = (data->seqlen[i] - initspan + 1);
			do {
				sites[i] = (int)(Random() * range);
			} while (!data->isValidSite[initspan][i][sites[i]]);
			
			if(DEBUG0) {
				if(sites[i] < 0 || sites[i] >= range) {
					fprintf(stderr, "Error: random generator error at setRandomSites()");
					exit(1);
				}
			}
		}
		else {
			sites[i] = EMPTY_SITE;
		}
	}
}


int addSite(int newsite, int seqind, int *sites, int numsites, Profile *countmat, Dataset *data) {
	if(DEBUG0) {
		assert(sites[seqind] == EMPTY_SITE);
	}
	if(DEBUG0) {
		//ensure that all positions considered does not have GAP_CHAR (N, W, etc)
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				for( int m =0; m < countmat->span; m++) {
					if(data->seqs[i][(m+sites[i])] >= NUMALPHAS) {
						fprintf(stderr, "seq %d, site %d, alpha %d", i, sites[i], data->seqs[i][(m+sites[i])]);
						exit(1);
					}
				}
			}
		}
		if(newsite >= 0) {
			for( int m =0; m < countmat->span; m++) {
				if(data->seqs[seqind][(m+newsite)] >= NUMALPHAS) {
					fprintf(stderr, "seq %d, site %d, alpha %d", seqind, newsite, data->seqs[seqind][(m+newsite)]);
					exit(1);
				}
			}
		}
	}

	if(newsite >= 0) { 
		//non-empty site - modify countmat, modify sites, and increment numsites
		int m;
		for( m =0; m < countmat->span; m++) {
			countmat->mat[m][ ((int)data->seqs[seqind][(m+newsite)]) ] += 1.0;
		}
		sites[seqind] = newsite;

		return (numsites + 1);
	}
	else {
		//do nothing
		sites[seqind] = EMPTY_SITE;
		return numsites;
	}
}

int removeSite(int seqind, int *sites, int numsites, Profile *countmat, Dataset *data) {
	if(DEBUG0) {
		//ensure that all positions considered does not have GAP_CHAR (N, W, etc)
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				for( int m =0; m < countmat->span; m++) {
					assert(data->seqs[i][(m+sites[i])] < NUMALPHAS);
				}
			}
		}
	}

	if(sites[seqind] >= 0) { 
		//non-empty site - modify countmat, modify sites, decrement numsites
		for(int m =0; m < countmat->span; m++) {
			countmat->mat[m][ ((int)data->seqs[seqind][(m+sites[seqind])]) ] -= 1.0;
		}
		sites[seqind] = EMPTY_SITE;
		return (numsites - 1);
	}
	else { 
		//do nothing
		sites[seqind] = EMPTY_SITE;
		return numsites;
	}
}


//2008-11-05 Algorithm from Uri
//This algorithm chooses the best sequence before choosing the best location within a sequence.
//It chooses the minimum according to the posterior P(Y_i = 0) over all each sequence i.
//This is eqivalent to choosing the highest average likelihood ratio of 
//each sequence (without using the ZOOPS term).
void findBestSitesByAvgLR(Markov *markov, Profile *freqmat, int *sites, Dataset *data) {
	//construct scoremat
	double scoremat[SPAN_CAPACITY][NUMCHARS]; //NUMCHARS include GAP_CHAR
	for(int x = 0; x < freqmat->span; x++) {
		for(int y = 0; y < NUMALPHAS; y++) {
			scoremat[x][y] = freqmat->mat[x][y];
		}
	}
	for(int x = 0; x < freqmat->span; x++) {
		scoremat[x][GAP_CHAR] = 0.0;
	}

	int old_numsites = getNumSites(sites, data->numseqs);
	double *scores = (double*) malloc(sizeof(double) * data->numseqs);

	//Step 1: compute avg-LR for each sequence and compute best location for each sequence
	for(int i=0; i < data->numseqs; i++) {
		if(data->numValidSites[freqmat->span][i] == 0) {
			scores[i] = 0.0;
			sites[i] = EMPTY_SITE;
		}
		else {
			double sum = 0.0;
			int maxind = -1;
			double maxscore = -DBL_MAX;
			for(int j=0; j < data->seqlen[i]- freqmat->span +1; j++) {
				double scoreCur = (double) data->isValidSite[freqmat->span][i][j] ;
				for(int m=0; m < freqmat->span; m++) {
					scoreCur *= scoremat[m][ data->seqs[i][j+m] ];
				}
				scoreCur /= markov->bgscore[freqmat->span][i][j];
				sum += scoreCur;
				if(maxscore < scoreCur) {
					maxscore = scoreCur;
					maxind = j;
				}
			}
			sites[i] = maxind;
			scores[i] = log((double)sum) - log((double)data->numValidSites[freqmat->span][i]); 

			if(DEBUG0) {
				if(!isValidConcatPos(data, i, maxind, freqmat->span)) {
					fprintf(stderr, "Error: Invalid site used as best site.\n");
					exit(1);
				}
			}
		}
	}

	//Step 2: choose best scoring sequences (highest avg LR)
	//This routine is probably faster using a heap, but we are only doing this once
	//througout the whole execution.
	//for(int k = 0; k < data->numseqs - numsites; k++) {
	while(getNumSites(sites, data->numseqs) > old_numsites) {
		//Eliminate the worst sequences
		double minval = DBL_MAX;
		int minind = -1;
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				if(scores[i] < minval) {
					minval = scores[i];
					minind = i;
				}
			}
		}
		
		//remove the site
		sites[minind] = EMPTY_SITE;

		if(DEBUG0) {
			if(minind == -1) {
				fprintf(stderr, "Error: cannot choose best scoring sequences.");
				exit(1);
			}
		}
	}

	if(getNumSites(sites, data->numseqs) != old_numsites) {
		fprintf(stderr, "Error: number of new sites is not the same as number of old sites!");
		exit(1);
	}

	free(scores);
}



//Find best sites with fixed freqmat and "number of sites"
//May have problems with number of valid sites == 0 (11/28/2008)
//void findBestSites(Markov *markov, Profile *freqmat, int *sites, Dataset *data) {
//	double scoremat[SPAN_CAPACITY][NUMCHARS]; //NUMCHARS include GAP_CHAR
//	for(int x = 0; x < freqmat->span; x++) {
//		for(int y = 0; y < NUMALPHAS; y++) {
//			scoremat[x][y] = freqmat->mat[x][y];
//		}
//	}
//	for(int x = 0; x < freqmat->span; x++) {
//		scoremat[x][GAP_CHAR] = 0.0;
//	}
//
//	int numsites = getNumSites(sites, data->numseqs);
//	double *scores = (double*) malloc(sizeof(double) * data->numseqs);
//
//	//Step 1: pick best site among all sequences
//	for(int i=0; i < data->numseqs; i++) {
//		int maxind = -1; 
//		double maxval = -DBL_MAX;
//		for(int j=0; j < data->seqlen[i]- freqmat->span +1; j++) {
//			double scoreCur = (double) data->isValidSite[freqmat->span][i][j] ;
//			for(int m=0; m < freqmat->span; m++) {
//				scoreCur *= scoremat[m][ data->seqs[i][j+m] ];
//			}
//			scoreCur /= markov->bgscore[freqmat->span][i][j];
//			if(scoreCur > maxval) {
//				maxval = scoreCur;
//				maxind = j;
//			}
//		}
//		sites[i] = maxind;
//		scores[i] = log((double)maxval) - log((double)data->numValidSites[freqmat->span][i]); //normalization term for ZOOPS score
//
//		if(DEBUG0) {
//			if(!isValidConcatPos(data, i, maxind, freqmat->span)) {
//				fprintf(stderr, "Error: Invalid site used as best site.\n");
//				exit(1);
//			}
//		}
//	}
//
//	//Step 2: choose best scoring sequences
//	//This routine is probably faster using a heap, but we are only doing this once
//	//througout the whole execution.
//	//Eliminate the worst sequences
//
//	for(int k = 0; k < data->numseqs - numsites; k++) {
//		double minval = DBL_MAX;
//		int minind = -1;
//		for(int i = 0; i < data->numseqs; i++) {
//			if(sites[i] >= 0) {
//				if(scores[i] < minval) {
//					minval = scores[i];
//					minind = i;
//				}
//			}
//		}
//		
//		//remove the site
//		sites[minind] = EMPTY_SITE;
//
//		if(DEBUG0) {
//			if(minind == -1) {
//				fprintf(stderr, "Error: cannot choose best scoring sequences.");
//				exit(1);
//			}
//		}
//	}
//
//	if(getNumSites(sites, data->numseqs) != numsites) {
//		fprintf(stderr, "Error: number of new sites is not the same as number of old sites!");
//		exit(1);
//	}
//
//	free(scores);
//}


//----------------------------------------------------------------------
// Phase Shifting
//----------------------------------------------------------------------
enum PhaseShift {NO_SHIFT, LEFT, RIGHT}; 

/**
* Uses Metropolis algorithm on entropy of the shifted column
* 
* param: 
*   countmat - count matrix to apply phase shift
*   sites - sites to apply phase shift
**/

void attemptPhaseShift(Markov *markov, Dataset *data, Profile *countmat, int *sites) {
	//equal chance for left or right
	enum PhaseShift shift = ((Random() < 0.5) ? LEFT : RIGHT); 

	double newcol[NUMALPHAS]; //newly attempted column to be used when shifted
	double *oldcol = NULL; //old column to be remove out when shifted
	for(int a = 0; a < NUMALPHAS; a++) {
		newcol[a] = 0.0;
	}

	double newscore = 0.0;
	double oldscore = 0.0;

	if(shift == LEFT) {
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				int newcol_pos = sites[i] - 1;
				if(newcol_pos >= 0 && data->isValidSite[countmat->span][i][sites[i]-1]) {
					newcol[data->seqs[i][newcol_pos]] += 1.0;
				}
				else {
					//out of bound -- cannot shift left
					return;
				}

				//subtracting background
				newscore -= log((double)markov->bgscore[countmat->span][i][sites[i]-1]);
				oldscore -= log((double)markov->bgscore[countmat->span][i][sites[i]]);
			}
		}
		//last column is removed during a left-shift
		oldcol = countmat->mat[countmat->span-1]; 
	}
	else if(shift == RIGHT) {
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) {
				int newcol_pos = sites[i] + countmat->span;
				if( newcol_pos < data->seqlen[i] && data->isValidSite[countmat->span][i][sites[i]+1]) {
					newcol[data->seqs[i][newcol_pos]] += 1.0;
				}
				else {
					//out of bound -- cannot shift right
					return;
				}

				//subtracting background
				newscore -= log((double)markov->bgscore[countmat->span][i][sites[i]+1]);
				oldscore -= log((double)markov->bgscore[countmat->span][i][sites[i]]);
			}
		}
		//first column is removed during a right-shift
		oldcol = countmat->mat[0];
	}

	for(int a = 0; a < NUMALPHAS; a++) {
		if(newcol[a] > 0.5) { //non-zero
			newscore += newcol[a] * log((double)newcol[a]);
		}
		if(oldcol[a] > 0.5) {
			oldscore += oldcol[a] * log((double)oldcol[a]);
		}
	}

	//if newent > oldent, then shift
	//if newent < oldent, then shift with prob exp(newent-oldent) 
	//[could introduce a temperature here]
	shift = (Random() < exp(newscore - oldscore) ? shift : NO_SHIFT);

	if(shift != NO_SHIFT){
		//diffShift is +1 for right shift; -1 for left shift
		int diffShift = (shift == LEFT ? -1 : +1);

		//shift count matrix
		if(shift == LEFT) {
			shiftProfileLeft(countmat, newcol);
		}
		else if(shift == RIGHT) {
			shiftProfileRight(countmat, newcol);
		}

		//shift sites
		for(int i = 0; i < data->numseqs; i++) {
			if(sites[i] >= 0) { //non-empty sites
				sites[i] = sites[i] + diffShift; 
			}
		}
	}

	if(DEBUG0) {
		int numsites = getNumSites(sites, data->numseqs);
		double sumold = 0.0;
		double sumnew = 0.0;
		for(int a = 0; a < NUMALPHAS; a++) {
			sumold += oldcol[a];
			sumnew += newcol[a];
		}
		assert(fabs(sumold - numsites) < 0.0000000001);
		assert(fabs(sumnew - numsites) < 0.0000000001);
	}

	if(DEBUG1) {
		fprintf(stderr, "shift with probability: %lf\n", exp(newscore - oldscore));
		fprintf(stderr, "newscore: %.6lg, oldscore: %.6lg\n", newscore, oldscore);
	}
}


//----------------------------------------------------------------------
// Combining isValid and bgscore
//----------------------------------------------------------------------
double*** buildValidBgscore(int **isValidSite[SPAN_CAPACITY], double ***bgscore, int minspan, int maxspan, Dataset *data) {
	double*** validBgscore;
	validBgscore = (double***) malloc(sizeof(double**) * SPAN_CAPACITY);
	for(int m = minspan; m <= maxspan; m++) {
		validBgscore[m] = (double**) malloc(sizeof(double*) * data->numseqs);
		for(int i=0; i < data->numseqs; i++) {
			validBgscore[m][i] = (double*) malloc(sizeof(double) * data->seqlen[i]);
		}
	}
	for(int m = minspan; m <= maxspan; m++) {
		for(int i=0; i < data->numseqs; i++) {
			for(int j = 0; j < data->seqlen[i]; j++) {
				if(isValidSite[m][i][j]) {
					validBgscore[m][i][j] = 1.0 / bgscore[m][i][j];
				}
				else {
					validBgscore[m][i][j] = 0.0;
				}
			}
		}
	}
	return validBgscore;
}
void nilValidBgscore(double ***validBgscore, int minspan, int maxspan, int numseqs) {
	for(int m = minspan; m <= maxspan; m++) {
		for(int i=0; i < numseqs; i++) {
			free(validBgscore[m][i]);
		}
		free(validBgscore[m]);
	}
	free(validBgscore);
}

