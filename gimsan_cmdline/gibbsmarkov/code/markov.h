#ifndef _MARKOV_H
#define _MARKOV_H

#include "stdinc.h"
#include "dataset.h"
#include "transprob.h"


//BG_BIOPRO - BioProspector
//BG_MOTIFSAMPLER - Thijs, et al. "A higher-order background model ..."
enum BackgroundModel { BG_INVALID, BG_GIBBSMARKOV, BG_BIOPRO, BG_MOTIFSAMPLER, BG_GMEAN, BG_AMEAN};
enum BackgroundDimType { BGDIM_AGGREGATE, BGDIM_NUMSEQS, BGDIM_WND};


typedef struct {
	int maxorder; //max order counted/calculated
	bool validSpan[SPAN_CAPACITY];

	//INFINITY for invalid sites (non-ACGT), NAN for out of bound
	double ***bgscore; //dimension: [SPAN_CAPACITY][numseqs][seqlen]

	//positional transition-probability [markov-order][numseqs][seqlen]
	//For invalid markov-order, the transprob with largest valid markov-order is used
	//For bad positions (GAP_CHAR), the transprob is 0.0.
	double ***posTransprob;

	int **minBeginDependInd; //dimension: [numseqs][seqlen]
	int numseqs; //dim of bgscore
	int *seqlen; //dim of bgscore

	double bgmodelFreq[NUMALPHAS];

	enum BackgroundModel bgmodel;
} Markov;

extern Markov* initMarkov(Dataset *bgfile, Dataset *data, int maxorder, int minspan, int maxspan,
						  enum BackgroundModel bgmodel, enum BackgroundDimType bgdimtyp);
extern void nilMarkov(Markov *markov);



#endif
