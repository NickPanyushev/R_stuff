#include "psprior.hpp"

#define print_error(str) for(fprintf(stderr,"%s\n",str); TRUE; exit(1))


static
FILE *openFile(char *filename)
{
	FILE *fptr;
	if((fptr = fopen(filename,"r")) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n",filename);
		print_error("File does not exist!\n");
	}
	return fptr;
}


static
void printPriors(FILE *fptr, Psprior *psp) {
	for(int i = 0; i < psp->numseqs; i++) {
		for(int j = 0; j < psp->seqlen[i]; j++) {
			fprintf(fptr, "%.4lg ", psp->prior[i][j]);
		}
		fprintf(fptr, "\n");
	}
}

static
bool isFloatingPoint(char c) {
	if(c >= '0' && c <= '9') {
		return true;
	}
	else if(c == 'E' || c == 'e' || c == '-' || c == '.') {
		return true;
	}
	else {
		return false;
	}
}

static
bool isSpaceChar(char c) {
	return (c == ' ' || c == '\n' || c == '\r' || c == '\t' );
}

//---------------------------------------------------------------------------------
// Converting to Reverse-complementary (double strand)
//---------------------------------------------------------------------------------
static
void augumentPspriorWithRevcompl(Psprior *psp) {
	if(!psp->useRevcompl) {
		printf("Error: invalid usage of revcompl for psp\n");
		abort();
	}

	//concatenate reverse complement to each strand
	for(int i = 0; i < psp->numseqs; i++) {
		//fprintf(stderr, "seqlen[%d]=%d\n", i, psp->seqlen[i]);
		//should have already allocated enough space
		for(int j = psp->seqlen[i]; j < 2*psp->seqlen[i] - psp->span + 1; j++) {
			psp->prior[i][j] = psp->prior[i][ 2*psp->seqlen[i] - j - psp->span];
		}
		for(int j = 2*psp->seqlen[i] - psp->span + 1; j < 2*psp->seqlen[i]; j++) {
			psp->prior[i][j] = 0.0;
		}
	}

	//modify length of sequence (2 * current length)
	for(int i = 0; i < psp->numseqs; i++) {
		psp->seqlen[i] = psp->seqlen[i] * 2;
	}

}

//
//---------------------------------------------------------------------------------
// Reading sequences from file input
//---------------------------------------------------------------------------------

//returns the length of the sequence
//This function can easily leads to buffer overflow, but I don't care about security now.
static
int readSinglePspSeq(FILE *fptr, Psprior *psp, int seqInd) {
	//int siteInd = 0;
	double flt = NAN;
	int c;
	for(int i = 0; i < psp->seqlen[seqInd]; i++) {
		fscanf (fptr, "%lf", &flt);
		psp->prior[seqInd][i] = flt;
	}
	while((c = fgetc(fptr)) != '>' && c!= EOF) {
		if( isSpaceChar(c)) {
			continue;
		}
		else {
			fprintf(stderr, "ERROR: invalid character found at end of sequence %c", c);
		}
	}
	ungetc(c, fptr);

	return psp->seqlen[seqInd];
}

static
void readPspSeqs(FILE *fptr, Psprior *psp) {
	char c;
	char prevC = 0; //null
	int seqInd = 0;
	boolean inHeader = FALSE;
	while ((c = fgetc(fptr)) != EOF) {
		if( c == '>' ) {
			//Allows for >> headers
			if(inHeader && prevC != '>') {
				print_error("Error: There is non-consecutive > in header (readSeqs).\n");
			}
			inHeader = TRUE;
			if(seqInd >= psp->numseqs) {
				print_error("ERROR in readSeqs: number of seqs exceed max.");
			}
		}
		else if ( inHeader && c == '\n') {
			inHeader = FALSE;
			int seqlen = readSinglePspSeq(fptr, psp, seqInd);
			if(seqlen != psp->seqlen[seqInd]) {
				print_error("ERROR: length of seqs does not match");
			}
			seqInd++;
		}
		prevC = c;
	}
	if(seqInd != psp->numseqs) {
		print_error("ERROR in readSeqs: Inconsistent number of sequences.");
	}
}


static
int getNumOfPspSeqs(FILE *fptr) {
	char c;
	char prevC = 0;
	int seqInd = 0;
	while ((c = fgetc(fptr)) != EOF) {
		if( c == '>' && prevC != '>') {
			seqInd++;
		}
		prevC = c;
	}
	return seqInd;
}

static
void getLenOfPspSeqs(FILE *fptr, Psprior *psp) {
	char c;
	char prevC = 0;
	int seqInd = -1;
	boolean inHeader = FALSE;
	boolean inSeq = FALSE;
	while ((c = fgetc(fptr)) != EOF) {
		if( c == '>' ) {
			inHeader = TRUE;
			inSeq = FALSE;
			if(prevC != '>') {
				seqInd++; //do not increment if >>  in header
			}
			if(seqInd >= psp->numseqs) {
				print_error("ERROR in getLenOfSeqs: number of seqs exceed max.");
			}
			psp->seqlen[seqInd] = 0;
		}
		else if ( inHeader && c == '\n') {
			inHeader = FALSE;
			inSeq = TRUE;
		}
		else if (inSeq) {
			if(isFloatingPoint(prevC) && isSpaceChar(c)) {
				psp->seqlen[seqInd]++;
				//fprintf(stderr, "case1\n");
			}
			else if(isSpaceChar(prevC) && isSpaceChar(c)) {
				//do nothing
				//fprintf(stderr, "case2\n");
			}
			else if(isSpaceChar(prevC) && isFloatingPoint(c)) {
				//do nothing
				//fprintf(stderr, "case3\n");
			}
			else if(isFloatingPoint(prevC) && isFloatingPoint(c)) {
				//do nothing
				//fprintf(stderr, "case4\n");
			}
			else {
				fprintf(stderr, "character: %d %d\n", (int)prevC, (int)c);
				fprintf(stderr, "seqind: %d\n", seqInd);
				fprintf(stderr, "Invalid format found: %c%c\n", prevC, c);
				exit(1);
			}
		}
		prevC = c;
	}
}


static
Psprior* _readPsprior(char *filename, bool useRevcompl, int span) {
	Psprior *psp;
	FILE *fptr;

	if(DEBUG1) {
		fprintf(stderr, "PSP open: %s\n", filename);
	}

	fptr = openFile(filename);
	psp = (Psprior*) malloc(sizeof(Psprior));
	psp->useRevcompl = useRevcompl;
	psp->span = span;


	//first character is '>', which is a FASTA format style
	if(fgetc(fptr) != '>') {
		fprintf(stderr, "Error: %s does not start with >\n", filename);
		exit(1);
	}

	//number of sequence
	rewind(fptr);
	psp->numseqs = getNumOfPspSeqs(fptr);

	psp->seqlen = (int*) malloc( psp->numseqs * sizeof(int));
	psp->prior = (double**) malloc(psp->numseqs * sizeof(double*));

	//len per seq
	rewind(fptr);
	getLenOfPspSeqs(fptr, psp);

	for(int i = 0; i < psp->numseqs; i++) {
		int len;

		//reverse-complementary is twice as long after concatenating
		//note that we do not change seqlen[i] here (this is done in augumentPspriorWithRevcompl())
		if(psp->useRevcompl) {
			len = psp->seqlen[i] * 2;
		}
		else {
			len = psp->seqlen[i];
		}
		psp->prior[i] = (double*) malloc(len * sizeof(double));
	}

	//read sequences (core of the function)
	rewind(fptr);
	readPspSeqs(fptr, psp);
	fclose(fptr);

	if(DEBUG1) {
		fprintf(stderr, "FASTA closed\n");
	}

	if(useRevcompl) {
		if(DEBUG1) {
			fprintf(stderr, "PSP with revcompl...\n");
		}
		augumentPspriorWithRevcompl(psp);
	}


	if(DEBUG1) {
		fprintf(stderr, "non-normalized\n");
		printPriors(stderr, psp);
	}


	for(int i = 0; i < psp->numseqs; i++) {
		double sum = 0.0;
		for(int j = 0; j < psp->seqlen[i]; j++) {
			sum += psp->prior[i][j];
		}
		for(int j = 0; j < psp->seqlen[i]; j++) {
			psp->prior[i][j] /= sum ;
		}
	}
	if(DEBUG1) {
		fprintf(stderr, "normalized\n");
		printPriors(stderr, psp);
	}

	//Find max/min of sequence length
	psp->maxseqlen = 0;
	for(int i = 0; i < psp->numseqs; i++) {
		if(psp->maxseqlen < psp->seqlen[i]) {
			psp->maxseqlen = psp->seqlen[i];
		}
	}
	psp->minseqlen = INT_MAX;
	for(int i = 0; i < psp->numseqs; i++) {
		if(psp->minseqlen > psp->seqlen[i]) {
			psp->minseqlen = psp->seqlen[i];
		}
	}

	return psp;
}

//---------------------------------------------------------------------------------
// Main functions
//---------------------------------------------------------------------------------

Psprior *openPsprior(char *filename, bool useRevcompl, int span) {
	//This only reserves the space for revcompl (does not actually set the data or
	//increase the seqlen)
	Psprior *psp = _readPsprior(filename, useRevcompl, span);

	return psp;
}


void nilPsprior(Psprior *psp) {
	free(psp->seqlen);

	for(int i = 0; i< psp->numseqs; i++) {
		free(psp->prior[i]);
	}
	free(psp->prior);
	free(psp);
}

void incorporatePspriorToBgscore(Markov *markov, Dataset *data, Psprior *psp){
	//sanity check
	assert(markov->numseqs == psp->numseqs);
	for(int i = 0; i < markov->numseqs; i++) {
		assert(markov->seqlen[i] == psp->seqlen[i]);
	}
	assert(markov->validSpan[psp->span]);

	for(int i = 0; i < markov->numseqs; i++) {
		for(int j = 0; j < markov->seqlen[i]; j++) {
			//Divide the background score by psp because bgscore is the denominator in the LLR score.
			//The numValidSites is present for backward compatibility. The original ZOOPS model divides
			//the ZOOPS sampling score by numValidSites (see zeroOccurScore in entsamp.c as well as computeEntropy).

			markov->bgscore[psp->span][i][j] /= (psp->prior[i][j] * data->numValidSites[psp->span][i]);
		}
	}

}



