#include "em_alg.h"

#define THRSHLD 0.01

			
//--------------------------------------------------------------------------------------------
// EM algorithm (Markov variant)
//--------------------------------------------------------------------------------------------

static
double EM_step(Markov *markov, Profile *freqmat, Dataset *data, double *posScore, double pseudocount[NUMALPHAS]) {
	double ilrval = 0.0;
	double scoremat[SPAN_CAPACITY][NUMCHARS]; //NUMCHARS include GAP_CHAR
	double newmat[SPAN_CAPACITY][NUMCHARS]; //NUMCHARS include GAP_CHAR

	for(int x = 0; x < freqmat->span; x++) {
		for(int y = 0; y < NUMALPHAS; y++) {
			scoremat[x][y] = freqmat->mat[x][y];
			newmat[x][y] = 0.0;
		}
	}
	for(int x = 0; x < freqmat->span; x++) {
		scoremat[x][GAP_CHAR] = 0.0;
		newmat[x][GAP_CHAR] = 0.0;
	}

	//setProfile(freqmat, 0.0);
	for(int i=0; i< data->numseqs; i++) {
		//E-step
		double sum = 0.0;
		for(int pos=0; pos < data->seqlen[i] - freqmat->span + 1; pos++)  {
			posScore[pos] = (double) data->isValidSite[freqmat->span][i][pos];
			for(int x=0; x< freqmat->span; x++) {
				posScore[pos] *= scoremat[x][ data->seqs[i][pos+x] ];
			}
			posScore[pos] /= markov->bgscore[freqmat->span][i][pos];
			sum += posScore[pos]; 
		}

		//M-step
		for(int pos=0; pos < data->seqlen[i] - freqmat->span+1; pos++)  {
			//normalize positional score
			double norm = posScore[pos] / sum;

			for(int x=0; x < freqmat->span; x++) {
				newmat[x][ data->seqs[i][pos+x] ] += norm;
				//pswm is constructed from the normalized probability value of
				//each nucleotide in the words that are at good positions 
			}
		}

		//ilrval is not part of the EM, but can be easily computed from EM_step()
		ilrval += log((double)sum) - log((double)(data->seqlen[i] - freqmat->span + 1));
	}

	for(int x = 0; x < freqmat->span; x++) {
		for(int y = 0; y < NUMALPHAS; y++) {
			freqmat->mat[x][y] = newmat[x][y];
		}
	}

	//normalize to frequency matrix
	for(int x=0; x< freqmat->span; x++)  {
		if(!freqmat->isgap[x]) {
			double sumthisx = 0.0;
			for(int y=0; y < NUMALPHAS; y++) {
				sumthisx += freqmat->mat[x][y];
			}
			for(int y=0; y< NUMALPHAS; y++)   {
				freqmat->mat[x][y] /= sumthisx; 
			}

			//add pseudocount and then normalize again
			sumthisx = 0.0;
			for(int y=0; y < NUMALPHAS; y++) {
				//divide by numseqs because pseudocount are "counts"
				freqmat->mat[x][y] += (pseudocount[y] / data->numseqs); 
				sumthisx += freqmat->mat[x][y];
			}
			for(int y=0; y< NUMALPHAS; y++)   {
				freqmat->mat[x][y] /= sumthisx; 
			}
		}
		else {
			for(int y=0; y < NUMALPHAS; y++) { 
				freqmat->mat[x][y] = 1.0;
			}
		}
	}
	return ilrval;
}

void runEmSteps(Markov *markov, Profile *freqmat, Dataset *data, int steps, boolean fastMode) {
	double pseudocount[NUMALPHAS] = {0.0, 0.0, 0.0, 0.0};
	runEmSteps(markov, freqmat, data, steps, fastMode, pseudocount);
}

void runEmSteps(Markov *markov, Profile *freqmat, Dataset *data, int steps, boolean fastMode, double pseudocount[NUMALPHAS]) 
{
	if(steps == 0) {
		//pswm remains unchanged
		return;
	}

	int i;
	double score = 0.0, origscore = 0.0, bestscore = -DBL_MAX;
	double thrshld = THRSHLD;

	double *posScore = (double*) malloc(data->maxseqlen * sizeof(double));

	int stepcount = 0; 
	int plateau = 0;
	for( i = 0; i < steps; i++) {
		score = EM_step(markov, freqmat, data, posScore, pseudocount);

		if(i == 0) {
			origscore = score;
		}
		
		if(fastMode) {
			if( score  > thrshld + bestscore ) { //thrshld > 0
				bestscore = score;
				plateau = 0;
			}
			else {
				plateau++;
			}
			if(plateau > 2) {
				break;
			}
		}

		stepcount++;
	}
	free(posScore);

	if(DEBUG1) {
		fprintf(stderr, "Before EM: %.3lf. After EM: %.3lf. Steps: %d\n", 
			origscore, score, stepcount);
	}

}

