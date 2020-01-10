#ifndef _PRINT_RESULTS_H
#define _PRINT_RESULTS_H

#include "gibbs_util.h"

extern void printRunNode(RunNode *node, enum ScoreMetric metric, Markov *markov, Dataset *data, Zoops *zoops);
extern void printTopRankRunNodes(int numTop, RunSet*, enum ScoreMetric metric, Markov*, Dataset*, Zoops*, double pvalCutoff);


#endif
