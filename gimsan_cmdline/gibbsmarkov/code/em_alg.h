#ifndef _EM_ALG_H
#define _EM_ALG_H

#include "profile.h"
#include "dataset.h"
#include "gibbs_util.h"

extern void runEmSteps(Markov *markov, Profile *freqmat, Dataset *data, int steps, boolean fastMode);
extern void runEmSteps(Markov *markov, Profile *freqmat, Dataset *data, int steps, boolean fastMode, double pseudocount[NUMALPHAS]);

#endif
