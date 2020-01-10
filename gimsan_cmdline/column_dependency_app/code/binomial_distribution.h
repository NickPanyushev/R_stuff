#ifndef _BINOMIAL_DISTRIBUTION_H
#define _BINOMIAL_DISTRIBUTION_H

#include "col_depend.h"


extern double binom_distrib_find_upper_bound_on_p(int trials, int successes, double alpha);
extern double binom_distrib_find_lower_bound_on_p(int trials, int successes, double alpha);

#endif
