#ifndef _PSPRIOR_H
#define _PSPRIOR_H

#include "stdinc.h"
#include "markov.h"
#include "dataset.h"

//static const int FASTA_HEADER_CAPACITY = 80;


typedef struct {
	//sequence and site index starts at 0.
	//the prior[][] and seqlen is after revcompl
	double  **prior; //seqs[numseqs][seqlen[i]] = {0, 1, 2, 3}
	//double *sumPrior;
	int	numseqs; //number of sequences
	int *seqlen; //length of each sequence
	int maxseqlen;
	int minseqlen;
	int span;

	//reverse-complementary
	bool useRevcompl;
} Psprior;

extern Psprior *openPsprior(char *filename, bool useRevcompl, int span);
extern void nilPsprior(Psprior *);

//incorporate bgscore with psp (also perform sanity check)
extern void incorporatePspriorToBgscore(Markov *markov, Dataset *data, Psprior *psp);

//use functions from Dataset for concat-pos / double-strand-pos lookup


#endif
