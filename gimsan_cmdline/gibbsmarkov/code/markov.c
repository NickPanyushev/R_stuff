#include "markov.h"

//--------------------------------------------------------------------------------------------
// Calculate seqtransprob
//
//--------------------------------------------------------------------------------------------


//Count to compute transition probabilities
//== Parameters
// bgfile is the background sequences used to estimate input FASTA
static
TransProb* calculateSeqTransProbAggregate(int maxorder, Dataset *bgfile, const int bgPseudocount) {
	TransCount *transcount = (TransCount*) malloc(sizeof(TransCount));
	transcount->maxorder = maxorder;

	//tally to transcount
	resetTransCount(transcount, bgPseudocount);
	for(int i = 0; i < bgfile->numseqs; i++) {
		//counts are accumulated for all sequences
		tallyToTransCount(transcount, bgfile->seqs[i], bgfile->seqlen[i]);
	}

	//compute transprob by normalization
	TransProb *transprob = normalizeTransCount(transcount);
	free(transcount);
	return transprob;
}

//Only calculate the transition probabilites for the given sequence-index
static
TransProb* calculateSeqTransProbForCurrentSeq(int seqind, int maxorder, Dataset *bgfile, const int bgPseudocount) {
	TransCount *transcount = (TransCount*) malloc(sizeof(TransCount));
	transcount->maxorder = maxorder;

	//tally to transcount
	resetTransCount(transcount, bgPseudocount);
	tallyToTransCount(transcount, bgfile->seqs[seqind], bgfile->seqlen[seqind]);

	//compute transprob by normalization
	TransProb *transprob = normalizeTransCount(transcount);
	free(transcount);
	return transprob;
}


//Count to compute transition probabilities
//This is done in a naive algorithm. It's much faster to count using the highest order,
//and then sum up for the smaller order.
//== Parameters
// numseqs is the number of sequences of input FASTA.
// bgfile is the background sequences used to estimate input FASTA

static
TransProb** calculateSeqTransProb(int maxorder, Dataset *bgfile, enum BackgroundDimType bgdimtyp, int numseqs) {
	const int BG_PSEUDOCOUNT = 1;

	if(DEBUG0) {
		assert(maxorder <= MARKOV_ORDER_BOUND);
		assert(maxorder >= 0);
	}

	TransProb **seqtransprob = (TransProb**) malloc(sizeof(TransProb*) * numseqs);

	if(bgdimtyp == BGDIM_AGGREGATE) {
		TransProb *transprob = calculateSeqTransProbAggregate(maxorder, bgfile, BG_PSEUDOCOUNT);

		//memory copy to the first element
		seqtransprob[0] = (TransProb*) malloc(sizeof(TransProb));
		seqtransprob[0]->maxorder = maxorder;
		copyTransProb(transprob, seqtransprob[0]);

		//copy-by-reference to the other objects
		for(int i = 1; i < numseqs; i++) {
			//fprintf(stderr, "numseqs=%d, i=%d, size=%d\n", numseqs, i, (int)sizeof(TransProb));
			seqtransprob[i] = seqtransprob[0];
		}

		for(int i = 0; i < numseqs; i++) {
			copyTransProb(transprob, seqtransprob[i]);
		}
		free(transprob);
	}
	else if(bgdimtyp == BGDIM_NUMSEQS) {
		if(numseqs != bgfile->numseqs) {
			fprintf(stderr, "Error: Under -bgdim_numseqs mode, the number of sequences of ");
			fprintf(stderr, "input FASTA %d and bgfile %d must be the same.\n", numseqs, bgfile->numseqs);
		}
		for(int i = 0; i < numseqs; i++) {
			seqtransprob[i] = calculateSeqTransProbForCurrentSeq(i, maxorder, bgfile, BG_PSEUDOCOUNT);
		}
	}
	else {
		fprintf(stderr, "Invalid BackgroundDimType at calculateTransProb()\n");
		exit(1);
	}

	return seqtransprob;
}

static
void printSeqTransProb(FILE *fptr, TransProb **seqtransprob, int numseqs) {
	for(int i = 0; i < numseqs; i++) {
		fprintf(fptr, "Sequence %d:\n", i);
		printTransProb(fptr, seqtransprob[i], true);
		//fprintf(fptr, "\n");
	}
}

//--------------------------------------------------------------------------------------------
// Accessor functions
//
//--------------------------------------------------------------------------------------------



//Get transition probability based on sequence and site indices.
//Semantically, p(s(2); s(1), s(0), s(-1)) =  p(s(2); s(1), s(0))
//Two additional features within minBeginDependInd:
// (1) Takes care of N and W characters in the dataset. beginDependInd automatically set as MAX of the
//     closest previous "N" character and beginDependInd
// (2) Double strand. beginDependInd automatically set to MAX{beginDependInd, start-of-revcompl}
// (3) beginDependInd can be negative, but it means the same thing as 0.
//Param:
//	beginDependInd - this is the smallest site-index that the background depends on
//		e.g. fnp0(m+w+2, beginDepend = m+w) should compute for p(m+w+2; m+w+1, m+w)
//			for a 4th order markov chain, instead of p(m+w+2; m+w+1, m+w, m+w-1, m+w-2)
//		e.g. fnp0(m+w, beginDepend = m+w) should compute p(m+w)
static
double getTransProb(Markov *markov , Dataset *data, int seqind, int siteind,
			 int beginDependInd)
{
	if(beginDependInd < markov->minBeginDependInd[seqind][siteind]) { //take max()
		beginDependInd = markov->minBeginDependInd[seqind][siteind];
	}

	if(DEBUG0) {
		//check that siteind is within bound
		assert(siteind < data->seqlen[seqind]);
		assert(siteind >=0);
		//check that beginDependInd <= siteind
		if(beginDependInd > siteind) {
			fprintf(stderr, "Error: beginDependInd > siteind\n");
			fprintf(stderr, "beginDependInd %d, siteind %d\n", beginDependInd, siteind);
			exit(1);
		}
		assert(beginDependInd >= 0);
	}


	//Set markov order.
	//Since beginDependInd >= 0,  condition p(s(2); s(1), s(0), s(-1)) =  p(s(2); s(1), s(0)) would
	//be satisfied using the technique.
	int order; //min{maxorder, siteind - beginDepend}
	if(markov->maxorder > siteind - beginDependInd) {
		order = siteind - beginDependInd;
	}
	else {
		order = markov->maxorder;
	}

	//fprintf(stderr, "seqind %d, siteind %d, beginDependInd %d, order %d\n", seqind, siteind, beginDependInd, order);

	return markov->posTransprob[order][seqind][siteind];

}

//--------------------------------------------------------------------------------------------
// Precompute Positional-Transition probability functions
//
//--------------------------------------------------------------------------------------------


static
int getActualMarkovOrder(int **minBeginDependInd, int seqind, int pos, int ord) {
	int i = seqind;
	int j = pos;

	//note that ord == 3 may have results of markov-order of 0 if j-1 is a bad character
	//j-ord may be negative, but getTransProb will take care of it to mean 0.
	//getTransProb() will output the correct probability as long as j is not a bad character.

	//Set markov order.
	//Condition p(s(2); s(1), s(0), s(-1)) =  p(s(2); s(1), s(0)) would be satisfied using the technique.
	int actualOrder; //adjusted markov-order in terms of bad-character, rev-compl, leftmost edge.
	if(ord > j - minBeginDependInd[i][j]) { //take min()
		actualOrder = j - minBeginDependInd[i][j];
	}
	else {
		actualOrder = ord;
	}
	if(DEBUG0) {
		if(actualOrder < 0) {
			fprintf(stderr, "minBeginDependInd[%d][%d] = %d\n", i, j, minBeginDependInd[i][j]);
			fprintf(stderr, "FATAL ERROR: seq %d site %d seems to be a bad character at getActualMarkovOrder(). actualOrder: %d\n", i, j, actualOrder);
		}
	}
	return actualOrder;
}

static
void precomputePosTransprob(Markov *markov, Dataset *data, TransProb **seqtransprob) {
	for(int ord = 0; ord <= markov->maxorder; ord++) {
		for(int i = 0; i < markov->numseqs; i++) {
			TransProb *transprob = seqtransprob[i];
			for(int j = 0; j < markov->seqlen[i]; j++) {
				if(data->isBadPos[i][j]) {
					markov->posTransprob[ord][i][j] = 0.0;
				}
				else {
					int actualOrder = getActualMarkovOrder(markov->minBeginDependInd, i, j, ord);

					//determine the previous alphabets of the word
					int alphas[MARKOV_ORDER_BOUND + 1];
					for(int a = 0; a <= actualOrder; a++) {
						alphas[a] = data->seqs[i][j - a];
					}

					//compute transition probability
					if(actualOrder == 0) {
						markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], -1, -1, -1, -1, -1);
					}
					else if(actualOrder == 1) {
						markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], -1, -1, -1, -1);
					}
					else if(actualOrder == 2) {
						markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], -1, -1, -1);
					}
					else if(actualOrder == 3) {
						markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], alphas[3], -1, -1);
					}
					else if(actualOrder == 4) {
						markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], alphas[3], alphas[4], -1);
					}
					else if(actualOrder == 5) {
						markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], alphas[3], alphas[4], alphas[5]);
					}
					else {
						printf("Error: invalid markov order at precomputePosTransprob(): %d\n", actualOrder);
						exit(1);
					}
				}
			}
		}
	}
}
/*
static
void precomputePosTransprobInWindows(Wndbgm *wndbgm, Markov *markov, Dataset *data) {
	if(wndbgm->bgWndSize < 1) {
		fprintf(stderr, "Error: bgWndSize must be at least 1.\n");
		exit(1);
	}
	if(wndbgm->maxorder != markov->maxorder) {
		fprintf(stderr, "Error: inconsistent maxorder between wndbgm and markov\n");
		exit(1);
	}

	for(int ord = 0; ord < markov->maxorder + 1; ord++) {
		for(int i = 0; i < markov->numseqs; i++) {
			for(int j = 0; j < markov->seqlen[i]; j++) {
				if(data->isBadPos[i][j]) {
					markov->posTransprob[ord][i][j] = 0.0;
				}
				else {
					int startpos = j - wndbgm->bgWndSize/2;
					startpos = (startpos >= 0 ? startpos : 0);
					int endpos = j + wndbgm->bgWndSize/2;
					endpos = (endpos <= markov->seqlen[i] ? endpos : markov->seqlen[i]); //endpos is non-inclusive within loop

					int count[NUMCHARS] = {0, 0, 0, 0, 0};
					for(int k = startpos; k < endpos; k++) {
						count[data->seqs[i][k]]++;
					}


					if(ord == 0) {
						//I am normalizing over ACGT only; one could normalize over all chars including GAP_CHAR
						int total = count[0] + count[1] + count[2] + count[3];
						markov->posTransprob[ord][i][j] = ((double)count[data->seqs[i][j]]) / total;
					}
					else {
						int actualOrder = getActualMarkovOrder(markov->minBeginDependInd, i, j, ord);

						//determine the previous alphabets of the word
						int alphas[MARKOV_ORDER_BOUND + 1];
						for(int a = 0; a <= actualOrder; a++) {
							alphas[a] = data->seqs[i][j - a];
						}

						TransProb *transprob = chooseModelForWndbgm(wndbgm, count[0], count[1], count[2], count[3]);

						//compute transition probability
						if(actualOrder == 0) {
							markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], -1, -1, -1, -1, -1);
						}
						else if(actualOrder == 1) {
							markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], -1, -1, -1, -1);
						}
						else if(actualOrder == 2) {
							markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], -1, -1, -1);
						}
						else if(actualOrder == 3) {
							markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], alphas[3], -1, -1);
						}
						else if(actualOrder == 4) {
							markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], alphas[3], alphas[4], -1);
						}
						else if(actualOrder == 5) {
							markov->posTransprob[ord][i][j] = getTransProb(transprob, actualOrder, alphas[0], alphas[1], alphas[2], alphas[3], alphas[4], alphas[5]);
						}
						else {
							printf("Error: invalid markov order at precomputePosTransprob(): %d\n", actualOrder);
							exit(1);
						}
					}
				}
			}
		}
	}

}
*/


//--------------------------------------------------------------------------------------------
// Precompution functions
//
//--------------------------------------------------------------------------------------------


//Its entries should be monotonically increasing across the sequence.
static
void computeMinBeginDependInd(Dataset *data, int **minBeginDependInd) {

	//reverse-complement
	//If the site is beyond the border after concatenation, the beginDependInd should
	//be at least the starting position of the backward strand.
	for(int i = 0; i < data->numseqs; i++) {
		if(data->useRevcompl) {
			int back_start = data->seqlen[i]/2; //starting position for backward strand
			//forward strand
			for(int j = 0; j < back_start; j++) {
				minBeginDependInd[i][j] = 0;
			}
			for(int j = back_start; j < data->seqlen[i]; j++) {
				minBeginDependInd[i][j] = back_start;
			}
		}
		else {
			for(int j = 0; j < data->seqlen[i]; j++) {
				minBeginDependInd[i][j] = 0;
			}
		}
	}

	for(int i = 0; i < data->numseqs; i++) {
		int index = 0;
		for(int j = 0; j < data->seqlen[i]; j++) {
			//For N and W characters. It should depend on the next character at this point.
			//It's a bit strange that j depends on j+1 here, but this actually is a good
			//idea as the code at getTransProb will catch such an inconsistent error and abort.
			if(data->isBadPos[i][j]) {
				index = j+1;
			}

			//take the max of the two indices
			if(minBeginDependInd[i][j] < index) {
				minBeginDependInd[i][j] = index;
			}
		}
	}
}

static
void precomputeBgscore(Markov *markov, Dataset *data, int span) {
	if(!markov->validSpan[span]) {
		printf("Error: proposed span is invalid.\n");
		exit(1);
	}

	//initialization
	double **score = markov->bgscore[span];

	for(int i = 0; i < data->numseqs; i++) {
		for(int j= 0; j < data->seqlen[i]; j++) {
			score[i][j] = NAN; //for out-of-bound (j >= seqlen[i] - span + 1)
		}
	}

	//precompute
	if(markov->bgmodel == BG_GIBBSMARKOV) {
		for(int i = 0; i < data->numseqs; i++) {
			for(int j= 0; j < data->seqlen[i] - span + 1; j++) {
				if(!data->isValidSite[span][i][j]) {
					score[i][j] = PINF;
				}
				else {
					score[i][j] = 1.0;
					for(int m = 0; m < span; m++) {
						score[i][j] *= getTransProb(markov,data, i, j+m, 0);
					}
					for(int k = 0; k < markov->maxorder; k++) { //exception from using "<= maxorder". See DEBUG0 check below.
						//If the motif is at the right edge of the forward-strand,
						//then the ratio will be 1.0 because j+span+k is at the reverse-strand.
						if(j+span+k < data->seqlen[i] && !data->isBadPos[i][j+span+k]) {
							score[i][j] *= getTransProb(markov,data, i, j+span+k, 0)
								/ getTransProb(markov, data, i, j+span+k, j+span);
						}
						else {
							break;
						}
					}

					if(DEBUG0) {
						int k = markov->maxorder;
						if(j+span+k < data->seqlen[i] && !data->isBadPos[i][j+span+k] ) {
							double ratio = getTransProb(markov,data, i, j+span+k, 0)
								/ getTransProb(markov, data, i, j+span+k, j+span);
							if(fabs(ratio - 1.0) > 0.000000001) {
								printf("Error: inconsistency in transition probability ratios\n");
								exit(1);
							}
						}
					}
				}
			}
		}
	}
	else if(markov->bgmodel == BG_BIOPRO) {
		for(int i = 0; i < data->numseqs; i++) {
			for(int j= 0; j < data->seqlen[i] - span + 1; j++) {
				if(!data->isValidSite[span][i][j]) {
					score[i][j] = PINF;
				}
				else {
					score[i][j] = 1.0;
					for(int m = 0; m < span; m++) {
						score[i][j] *= getTransProb(markov,data, i, j+m, j);
					}
				}
			}
		}
	}
	else if(markov->bgmodel == BG_MOTIFSAMPLER) {
		for(int i = 0; i < data->numseqs; i++) {
			for(int j= 0; j < data->seqlen[i] - span + 1; j++) {
				if(!data->isValidSite[span][i][j]) {
					score[i][j] = PINF;
				}
				else {
					score[i][j] = 1.0;
					for(int m = 0; m < span; m++) {
						score[i][j] *= getTransProb(markov,data, i, j+m, 0);
					}
				}
			}
		}
	}
	else if(markov->bgmodel == BG_GMEAN) {
		//similar to MOTIF_SAMPLER
		for(int i = 0; i < data->numseqs; i++) {
			for(int j= 0; j < data->seqlen[i] - span + 1; j++) {
				if(!data->isValidSite[span][i][j]) {
					score[i][j] = PINF;
				}
				else {
					score[i][j] = 1.0;
					for(int m = 0; m < span; m++) {
						double score1 = getTransProb(markov,data, i, j+m, 0);
						double score2 = getTransProb(markov, data, i, getConcatPosOfOppStrand(data, i, j+m, 1), 0);

						score[i][j] *= sqrt(score1 * score2);
					}
				}
			}
		}
	}
	else if(markov->bgmodel == BG_AMEAN) {
		//similar to MOTIF_SAMPLER
		for(int i = 0; i < data->numseqs; i++) {
			for(int j= 0; j < data->seqlen[i] - span + 1; j++) {
				if(!data->isValidSite[span][i][j]) {
					score[i][j] = PINF;
				}
				else {
					score[i][j] = 1.0;
					for(int m = 0; m < span; m++) {
						double score1 = getTransProb(markov,data, i, j+m, 0);
						double score2 = getTransProb(markov, data, i, getConcatPosOfOppStrand(data, i, j+m, 1), 0);

						score[i][j] *= (score1 + score2) / 2.0;
					}
				}
			}
		}
	}

	else {
		printf("Error: Invalid background type\n"); exit(1);
	}

}

//--------------------------------------------------------------------------------------------
// Core functions to use Markov probabilites
//
// Note: Do not add any function that requires headers other than dataset.h to avoid bloating!
//--------------------------------------------------------------------------------------------

//bgfile - dataset used to estimate markov probability
//data - dataset used to calculate precomputation of bgscore

Markov* initMarkov(Dataset *bgfile, Dataset *data, int maxorder, int minspan, int maxspan,
				   enum BackgroundModel bgmodel, enum BackgroundDimType bgdimtyp) {
	if(DEBUG0) {
		assert(maxorder <= MARKOV_ORDER_BOUND);
	}

	Markov *markov = (Markov*) malloc(sizeof(Markov));
	markov->bgmodel = bgmodel;
	markov->maxorder = maxorder;

	for(int i = 0; i < NUMALPHAS; i++) {
		markov->bgmodelFreq[i] =  bgfile->ntFreq[i];
	}

	//initialize validSpan
	for(int i = 0; i < SPAN_CAPACITY; i++) {
		markov->validSpan[i] = false;
	}
	for(int i = minspan; i <= maxspan; i++) {
		markov->validSpan[i] = true;
	}

	//compute transprob (using bgfile)
	TransProb **seqtransprob = NULL;
	if(bgdimtyp == BGDIM_AGGREGATE || bgdimtyp == BGDIM_NUMSEQS) {
		seqtransprob = calculateSeqTransProb(markov->maxorder, bgfile, bgdimtyp, data->numseqs);
		if(DEBUG1) {
			printSeqTransProb(stderr, seqtransprob, data->numseqs);
		}
	}
	//else if(bgdimtyp == BGDIM_WND) {
	//	seqtransprob = NULL;
	//}
	else {
		fprintf(stderr, "Invalid bgdimtyp");
		exit(1);
	}

	//construct bgscore dimension (which is the same dimension as data)
	markov->numseqs = data->numseqs;
	markov->seqlen = (int*) malloc(sizeof(int) * markov->numseqs);
	for(int i = 0; i < markov->numseqs; i++) {
		markov->seqlen[i] = data->seqlen[i];
	}

	markov->bgscore = (double***) malloc(sizeof(double**) * SPAN_CAPACITY);
	for(int m = minspan; m <= maxspan; m++) {
		markov->bgscore[m] = (double**) malloc(sizeof(double*) * data->numseqs);
		for(int i=0; i < data->numseqs; i++) {
			markov->bgscore[m][i] = (double*) malloc(sizeof(double) * data->seqlen[i]);
		}
	}
	markov->posTransprob = (double***) malloc(sizeof(double**) * (markov->maxorder+1));
	for(int ord = 0; ord <= markov->maxorder; ord++) {
		markov->posTransprob[ord] = (double**) malloc(sizeof(double*) * data->numseqs);
		for(int i=0; i < data->numseqs; i++) {
			markov->posTransprob[ord][i] = (double*) malloc(sizeof(double) * data->seqlen[i]);
		}
	}


	markov->minBeginDependInd = (int**) malloc(sizeof(int*) * data->numseqs);
	for(int i=0; i < data->numseqs; i++) {
		markov->minBeginDependInd[i] = (int*) malloc(sizeof(int) * data->seqlen[i]);
	}

	computeMinBeginDependInd(data, markov->minBeginDependInd);

	//print out debugging information
	if(DEBUG1) {
		fprintf(stderr, "numseqs %d\n", data->numseqs);

		for(int span = minspan; span <= maxspan; span++) {
			//printing isValidSite
			fprintf(stderr, "isValidSite for span %d\n", span);
			for(int i = 0; i < data->numseqs; i++) { //takeout
				fprintf(stderr, "seqlen[%d] %d\n", i, data->seqlen[i]);
				for(int j= 0; j < data->seqlen[i]; j++) {
					fprintf(stderr, "%d(%d) ", j, data->isValidSite[span][i][j]);
				}
				fprintf(stderr, "\n");
			}
			fprintf(stderr, "\n");
		}

		//printing isBadPos
		fprintf(stderr, "isBadPos\n");
		for(int i = 0; i < data->numseqs; i++) { //takeout
			fprintf(stderr, "seqlen[%d] %d\n", i, data->seqlen[i]);
			for(int j= 0; j < data->seqlen[i]; j++) {
				fprintf(stderr, "%d(%d) ", j, data->isBadPos[i][j]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");

		//printing minBeginInd
		fprintf(stderr, "minBeginDependInd\n");
		for(int i = 0; i < data->numseqs; i++) { //takeout
			fprintf(stderr, "seqlen[%d] %d\n", i, data->seqlen[i]);
			for(int j= 0; j < data->seqlen[i]; j++) {
				fprintf(stderr, "%d(%d) ", j, markov->minBeginDependInd[i][j]);
			}
			fprintf(stderr, "\n");
		}
		fprintf(stderr, "\n");
	}

	if(bgdimtyp == BGDIM_AGGREGATE || bgdimtyp == BGDIM_NUMSEQS) {
		precomputePosTransprob(markov, data, seqtransprob);
	}
	//else if(bgdimtyp == BGDIM_WND) {
	//	precomputePosTransprobInWindows(wndbgm, markov,data);
	//}
	else {
		fprintf(stderr, "Invalid bgdimtyp");
		exit(1);
	}
	if(DEBUG1) {
		//printing posTransprob
		for(int ord = 0; ord <= markov->maxorder; ord++) {
			fprintf(stderr, "posTransprob for markov-order %d\n", ord);
			for(int i = 0; i < data->numseqs; i++) { //takeout
				fprintf(stderr, "seqlen[%d] %d\n", i, data->seqlen[i]);
				for(int j= 0; j < data->seqlen[i]; j++) {
					fprintf(stderr, "%d(%.3lg) ", j, markov->posTransprob[ord][i][j]);
				}
				fprintf(stderr, "\n");
			}
			fprintf(stderr, "\n");
		}

	}

	//deallocate because it is only needed for posTransprob
	if(seqtransprob != NULL) {
		if(bgdimtyp == BGDIM_AGGREGATE) {
			free(seqtransprob[0]);
			//the other seqtransprob[i] are references
		}
		else if(bgdimtyp == BGDIM_NUMSEQS) {
			for(int i = 0; i < markov->numseqs; i++) {
				free(seqtransprob[i]);
			}
		}
		free(seqtransprob);
	}

	//precompute bgscore
	for(int m = minspan; m <= maxspan; m++) {
		if(DEBUG1) {
			fprintf(stderr, "Building markov background matrix for span %d\n", m);
		}

		precomputeBgscore(markov, data, m);

		if(DEBUG1) {
			int span = m;
			//printing bgscore
			fprintf(stderr, "bgscore\n");
			for(int i = 0; i < data->numseqs; i++) { //takeout
				fprintf(stderr, "seqlen[%d] %d\n", i, data->seqlen[i]);
				for(int j= 0; j < data->seqlen[i]; j++) {
					fprintf(stderr, "%d(%.2lg) ", j, markov->bgscore[span][i][j]);
				}
				fprintf(stderr, "\n");
			}
		}
	}

	return markov;
}

extern void nilMarkov(Markov *markov) {
	//do not deallocate Dataset *data
	for(int m = 0; m < SPAN_CAPACITY; m++) {
		if(markov->validSpan[m]) {
			for(int i=0; i < markov->numseqs; i++) {
				free(markov->bgscore[m][i]);
			}
			free(markov->bgscore[m]);
		}
	}
	for(int i=0; i < markov->numseqs; i++) {
		free(markov->minBeginDependInd[i]);
	}
	free(markov->minBeginDependInd);

	for(int ord = 0; ord <= markov->maxorder; ord++) {
		for(int i = 0; i < markov->numseqs; i++) {
			free(markov->posTransprob[ord][i]);
		}
		free(markov->posTransprob[ord]);
	}
	free(markov->posTransprob);

	free(markov->bgscore);
	free(markov->seqlen);
	free(markov);
}

