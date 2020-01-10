#ifndef _MARKOV_SCORE_H
#define _MARKOV_SCORE_H

#include "markov.h"
#include "zoops.h"

//--------------------------------------------------------------------------------------------
// Supplementary functions 
//
// Note: Do not add any function that requires headers other than dataset.h to avoid bloating!
//--------------------------------------------------------------------------------------------

// Possible generalization of "entropy score" under k-th order Markov background.
extern double computeEntropy(Markov *markov, Dataset *data, Zoops *zoops, int *sites, int span);
//extern double computeEntropy(Markov *markov, Dataset *data, Zoops *zoops, int *sites, int span, double pseudocount[NUMALPHAS]);

double computeEntropyPerSite(Markov *markov, Dataset *data, Zoops *zoops, int *sites, int seqind, int span);

extern double computeIlr(Markov *markov, Dataset*, Zoops *zoops, int *sites, int span);
extern double computeIlr(Markov *markov, Dataset*, double** pwm, int span);

#endif
