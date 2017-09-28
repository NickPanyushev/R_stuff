#include "entsamp.h"
#include "psprior.hpp"

enum IterMode { TRIAL, NORMAL }; //see updateOneIter



//----------------------------------------------------------------------
// Entropy Sampler functions
//----------------------------------------------------------------------

static
EntropySampTable* buildEntropySampTable(Dataset *data) {
	EntropySampTable* table = (EntropySampTable*) malloc(sizeof(EntropySampTable));
	table->maxcount = data->numseqs - 1; //maxcount in lookup table
	table->numalphas = NUMALPHAS;

	table->scoremat = (double**) malloc(sizeof(double*) * table->numalphas);
	for(int i = 0; i < table->numalphas; i++) {
		table->scoremat[i] = (double*) malloc(sizeof(double) * (table->maxcount+1));
	}

	for(int a = 0; a <  table->numalphas; a++) {
		table->scoremat[a][0] = 0.0; //modified for markov

		for(int c = 1; c <= table->maxcount; c++) {
			//modified for markov
			//(c+1) * log(c+1) - c*log(c) 
			table->scoremat[a][c] = (c+1.0)*log((double)(c+1.0)) - c*log((double)c);
		}
	}

	double maxval = -DBL_MAX;
	for(int a = 0; a <  table->numalphas; a++) {
		for(int c = 0; c <= table->maxcount; c++) {
			if(maxval < table->scoremat[a][c]) {
				maxval = table->scoremat[a][c];
			}
		}
	}

	//normalize
	for(int a = 0; a <  table->numalphas; a++) {
		for(int c = 0; c <= table->maxcount; c++) {
			table->scoremat[a][c] = exp( table->scoremat[a][c] - maxval);
		}
	}

	return table;
}

static
void nilEntropySampTable(EntropySampTable *table) {
	int i;
	for( i = 0; i < table->numalphas; i++) {
		free(table->scoremat[i]);
	}
	free(table->scoremat);
	free(table);
}


//----------------------------------------------------------------------
// Core functions
//----------------------------------------------------------------------
/**
* One iteration of Gibbs-sampling. This samples a new site for each sequence.
* 
* param: 
*   countmat - count matrix at the beginning of this iteration
*   sites - sites at the beginning of this iteration
*   mode - TRIAL is for quick sample and the entropy score is not computed; 
*          NORMAL otherwise
* return: 
*   returns the entropy score of the best pswm/sites.
**/
static
double updateOneIter(Profile *countmat, int *sites, Gibbs *gibbs, enum IterMode mode)
{
	//double *posScore = gibbs->posScore; //positional score
	double **table = gibbs->entsampTable->scoremat; //lookup table with dim [alpha, count]
	double scoremat[SPAN_CAPACITY][NUMCHARS]; //NUMCHARS include GAP

	//These entries for GAP_CHAR can be set to anything. "isValidSite" already 
	//makes positional scores be 0.0.
	for(int m = 0; m < countmat->span; m++) {
		scoremat[m][GAP_CHAR] = 0.0; 
	}

	//initial numsites for this iteration
	int numsites = getNumSites(sites, gibbs->data->numseqs);

	for(int i = 0; i < gibbs->data->numseqs; i++) {
		//double *bgscore = gibbs->markov->bgscore[countmat->span][i]; //new for markov
		//int *isValidSite = gibbs->data->isValidSite[countmat->span][i];
		double *validBgscore = gibbs->validBgscore[countmat->span][i];
		int *seq = gibbs->data->seqs[i]; 

		if(gibbs->data->numValidSites[countmat->span][i] == 0) {
			//Cannot pick a site by default (non-OOPS mode)
			sites[i] = EMPTY_SITE;
		}
		else {
			//remove site and remove site count from count-matrix
			numsites = removeSite(i, sites, numsites, countmat, gibbs->data); //return new number of sites

			////Choosing sampling technique
			//new for markov
			//does not divide by the usual 0th-order background
			if(gibbs->samptyp == CLR_SAMP) {
				for(int m = 0; m < countmat->span; m++) {
					for(int a = 0; a < gibbs->numalphas; a++) {
						scoremat[m][a] = table[a][(int)countmat->mat[m][a]];
					}
				}
			}
			else {
				for(int m = 0; m < countmat->span; m++) {
					double sum = 0.0;
					for(int a = 0; a < NUMALPHAS; a++) {
						scoremat[m][a] = countmat->mat[m][a] + (gibbs->data->ntFreq[a] * numsites * gibbs->pseudoweight);
						sum += scoremat[m][a];
					}
					for(int a = 0; a < NUMALPHAS; a++) {
						scoremat[m][a] /= sum;
					}
				}
				if(DEBUG0) {
					for(int m = 0; m < countmat->span; m++) {
						double sum = 0.0;
						for(int a = 0; a < NUMALPHAS; a++) {
							sum += scoremat[m][a];
						}
						if(fabs(sum - 1.0) > 0.000001) {
							fprintf(stderr, "Error: each column of scoremat does not sum up to 1.0\n");
							exit(1);
						}
					}
				}
			}

			//Given that there is a site
			double sumPosScore = 0.0;
			for(int j = 0; j < gibbs->data->seqlen[i] - countmat->span + 1; j++) {
				//Using "register" for fgscore and m is signficantly faster in my test
				//Most of the computational power goes into the following few lines
				//register double fgScore = (double) isValidSite[j]; //0 if invalid site; 1 if valid site
				register double fgScore = 1.0; 
				for(register int m = 0; m < countmat->span; m++) {
					fgScore *= scoremat[m][seq[j+m]];
				}
				//posScore[j] =  fgScore / bgscore[j];
				//sumPosScore += posScore[j];
				//sumPosScore += isValidSite[j] * fgScore / bgscore[j];
				sumPosScore += validBgscore[j] * fgScore;
				gibbs->cumsumArr[j] = sumPosScore;
			}

			//Score for no motif occurrence in sequence i
			double zeroOccurScore = 0.0;
			//the 'numsites' here is the number of sites excluding the current sequence i
			if(gibbs->zoops->isZoopsMode && numsites >= 3) { 
				//score is 0.0 if not in ZOOPS mode
				//why 3?
				//prevents from having less than 3 sites
				zeroOccurScore = gibbs->zoops->zeroOccurrenceSamplingScore[numsites] * gibbs->data->numValidSites[countmat->span][i];
			}

			if(DEBUG0) {
				//numsites is strictly less than numseqs because a site has been removed
				assert(numsites >= 0 && numsites < gibbs->data->numseqs);
			}

			//Sampling step
			int newsite;
			double drand = Random() * (sumPosScore + zeroOccurScore);
			if(gibbs->zoops->isZoopsMode && drand > sumPosScore) {
				newsite = EMPTY_SITE;
			}
			else {
				//double cumsum = 0.0;
				////if drand is a epsilon bigger than sum, then uses the last site
				//newsite = gibbs->data->seqlen[i] - countmat->span; 
				//for(int j = 0; j < gibbs->data->seqlen[i] - countmat->span + 1; j++) {
				//	cumsum += posScore[j];
				//	if(drand <= cumsum) {
				//		newsite = j;
				//		break;
				//	}
				//}

				//use binary search for cumsumArr
				//if cumsumArr[j-1] < drand <= cumsumArr[j], then the newsite should be j
				int cumsumArrLen = gibbs->data->seqlen[i] - countmat->span + 1;

				//if drand is a epsilon bigger than sum, then uses the last site
				newsite = gibbs->data->seqlen[i] - countmat->span; //default

				if(drand <= gibbs->cumsumArr[0]) {
					newsite = 0; //since no j-1 in this special case
				}
				else {
					int left = 0; 
					int right = cumsumArrLen - 1; //inclusive
					while(left + 1 < right) {
						int mid = left + (right - left) / 2;
						if(gibbs->cumsumArr[mid] < drand) {
							left = mid;
						}
						else {
							right = mid;
						}
					}
					newsite = right;
				}
			}
			//add new site and add site count to count-matrix
			numsites = addSite(newsite, i, sites, numsites, countmat, gibbs->data); //return new number of sites

			if(DEBUG1) {
				//fprintf(stderr, "seq: %3d, numsites: %3d, newsite: %4d, sumPos: %7lg, zeroOccur: %7lg, drand: %7lg\n", i, numsites, newsite, sumPosScore, zeroOccurScore, drand);
			}
		}
	}

	double entropy = -DBL_MAX; //entropy of this iteration
	if(mode == NORMAL) {
		entropy = computeEntropy(gibbs->markov, gibbs->data, gibbs->zoops, sites, countmat->span); //new for markov
	}


	return entropy;
}



/**
* One run of Gibbs-sampling. A run is given a set of starting positions (sites).
* 
* param: 
*   rnode - run-node that contains the new starting positions (sites)
*   initspan - initial span of the run
*   mode - TRIAL is for quick sample and the entropy score is not computed; 
*          NORMAL otherwise
**/

void runIters(RunNode *rnode, int initspan, Gibbs *gibbs) {
	int iterplat = 0; 
	int iters = 0;
	double score;
	double maxval = -DBL_MAX;

	//initialize temporary storage
	Profile *countmat = initProfile(initspan, gibbs->span);

	int *sites = (int*) malloc(gibbs->data->numseqs * sizeof(int));

	//set sites as random starting positions
	for(int i = 0; i < gibbs->data->numseqs; i++) {
		sites[i] = rnode->sites[i];
	}

	//the best are stored in rnode
	Profile *bestcountmat = rnode->countmat;
	int *bestsites = rnode->sites;

	updateCountmatFromSites(countmat, sites, gibbs->data);

	//Trial samplings
	for(int i = 0; i < gibbs->trialIters; i++) {
		updateOneIter(countmat, sites, gibbs, TRIAL);
		iters++;
	}


	//Normal samplings (the "real" samplings)
	while(1) {
		//rapid-convergence condition
		if( gibbs->useRapidConv && iterplat >= gibbs->iterPlateauLen) {
			break;
		}
		//fixed-iteration condition
		if( !gibbs->useRapidConv && iters >= gibbs->numFixIters) {
			break;
		}

		if(DEBUG1) {
			fprintf(stderr, "\nRun %d, iter %d, span %d\n", rnode->runId, 
				iters, countmat->span);
		}

		//sampling section
		score = updateOneIter(countmat, sites, gibbs, NORMAL);

		if(score > maxval) {
			maxval = score;
			rnode->score = score;
			copyProfile(bestcountmat, countmat);
			for(int i = 0; i < gibbs->data->numseqs; i++) {
				bestsites[i] = sites[i];
			}
			iterplat = 0;
		}
		else {
			iterplat++;
		}

		if(DEBUG1) {
			fprintf(stderr, "score %.4lf\n", score);
			for(int i = 0; i < gibbs->data->numseqs; i++) {
				fprintf(stderr, "sites[%d] %d\n", i, sites[i]);
			}
			printProfile(stderr, countmat);
		}

		//phase shift
		if(Random() < gibbs->phaseShiftFreq) {
			attemptPhaseShift(gibbs->markov, gibbs->data, countmat, sites);
		}

		if(DEBUG0) {
			assert(validCountmatWithSites(countmat, sites, gibbs->data));
			if(!gibbs->zoops->isZoopsMode) {
				assert(getNumSites(sites, gibbs->data->numseqs) == gibbs->data->numseqs);
			}
		}

		iters++; 
	}

	gibbs->totalIters += iters;

	free(sites);
	nilProfile(countmat);
}


//----------------------------------------------------------------------
// Constructor/Destructor
//----------------------------------------------------------------------

void initGibbsStructs(Gibbs *gibbs) {
	gibbs->randseed = labs(gibbs->randseed);
	if(gibbs->randseed == 0) {
		gibbs->randseed = 1; //cannot start with seed 0
	}
	sRandom(gibbs->randseed);

	if(DEBUG1) {
		fprintf(stderr, "Random seed: %u\n", gibbs->randseed);
	}

	gibbs->data = openDataset(gibbs->fastafile, gibbs->useRevcompl, gibbs->span, gibbs->span);
	if(gibbs->samptyp == CLR_SAMP) {
		gibbs->pseudoweight = 0.0;
	}

	//background
	Dataset *bgfile;
	if(gibbs->useBgfileOption) {
		//If reverse-complement is used for input, then reverse-complement is used
		//to estimate background probabilities. Background prob would be symmetric.
		bgfile = openBackgroundData(gibbs->bgfilename, gibbs->useRevcompl);
	}
	else {
		gibbs->bgfilename = gibbs->fastafile;
		bgfile = gibbs->data;
	}
	//Wndbgm *wndbgm = NULL;
	//if(gibbs->bgdimtyp == BGDIM_WND) {
	//	wndbgm = initWndbgm(gibbs->markovOrder, gibbs->markovModelsFilename, gibbs->bgWndSize);
	//	gibbs->numModelsForWndbgm = wndbgm->numModels; //save it for display only
	//}
	gibbs->markov = initMarkov(
			bgfile,
			gibbs->data,
			gibbs->markovOrder,
			gibbs->span,
			gibbs->span,
			gibbs->bgmodel,
			gibbs->bgdimtyp
			);
	//if(wndbgm != NULL) {
	//	nilWndbgm(wndbgm);
	//}
	if(gibbs->useBgfileOption) {
		nilDataset(bgfile); //only used to estimate probabilities
	}

	//PSP
	if(gibbs->pspfilename != NULL) {
		Psprior *psp = openPsprior(gibbs->pspfilename, gibbs->useRevcompl, gibbs->span);
		incorporatePspriorToBgscore(gibbs->markov, gibbs->data, psp);
		nilPsprior(psp);
	}

	//other options
	//gibbs->posScore = (double*) malloc(gibbs->data->maxseqlen * sizeof(double));
	gibbs->cumsumArr = (double*) malloc(gibbs->data->maxseqlen * sizeof(double));

	//combining isValid and bgscore
	gibbs->validBgscore = buildValidBgscore(gibbs->data->isValidSite, gibbs->markov->bgscore, gibbs->span, gibbs->span, gibbs->data);

	gibbs->runset = createRunSet();
	gibbs->entsampTable = buildEntropySampTable(gibbs->data);
	gibbs->zoops = initZoops(gibbs->isZoopsMode, gibbs->data->numseqs, gibbs->zoopsPseudoweight);

	if(!gibbs->isZoopsMode) {
		for(int i = 0; i < gibbs->data->numseqs; i++) {
			if(gibbs->data->numValidSites[gibbs->span][i] <= 0) {
				printf("Error: in OOPS mode, sequence %d does not have a site\n", i);
				exit(1);
			}
		}
	}

}

void nilGibbs(Gibbs *gibbs ) {
	//free(gibbs->posScore);
	free(gibbs->cumsumArr);
	nilEntropySampTable(gibbs->entsampTable);
	nilRunSet(gibbs->runset) ;
	nilValidBgscore(gibbs->validBgscore, gibbs->span, gibbs->span, gibbs->data->numseqs);
	nilDataset(gibbs->data);
	nilMarkov(gibbs->markov);
	nilZoops(gibbs->zoops);
	free(gibbs);
}
