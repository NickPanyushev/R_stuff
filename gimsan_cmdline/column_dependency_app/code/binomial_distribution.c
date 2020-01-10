#include "binomial_distribution.h"
#include "dcdflib.h"
//#include <boost/math/distributions/normal.hpp>

//distribution<>::find_lower_bound_on_p(numperms, colpair->numObsEntGe, (1.0-coldep->beta)/2);
//distribution<>::find_upper_bound_on_p(numperms, colpair->numObsEntGe, (1.0-coldep->beta)/2);
//
//

static
double my_quantile(double p) {
	//boost::math::normal s;
	//return boost::math::quantile(s, 1.0 - alpha);
	//
	
	double q = 1.0 - p;
	return dinvnr(&p, &q);
	
}

static
double compute_kappa(double alpha) {
	return my_quantile(1.0 - alpha);
}

static
double binom_distrib_interval(int trials, int successes, double alpha, int sign) {
	//sign should either be +1 or -1
	//Implements the Agresti-Coull confidence interval
	double kappa = compute_kappa(alpha);
	double x_tilde = successes + kappa * kappa / 2.0;
	double n_tilde = trials + kappa * kappa;
	double p_tilde = x_tilde / n_tilde;
	double value = p_tilde + sign * fabs(kappa * sqrt(p_tilde * (1-p_tilde)) / sqrt(n_tilde)) ;
	if(value < 0.0) {
		value = 0.0;
	}
	if(value > 1.0) {
		value = 1.0;
	}
	return value;
}


double binom_distrib_find_lower_bound_on_p(int trials, int successes, double alpha) {
	return binom_distrib_interval(trials, successes, alpha, -1) ;
}

double binom_distrib_find_upper_bound_on_p(int trials, int successes, double alpha) {
	return binom_distrib_interval(trials, successes, alpha, +1);
}

