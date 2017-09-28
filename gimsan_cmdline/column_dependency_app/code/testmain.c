#include <stdio.h>
#include <stdlib.h>
//#include <boost/math/distributions/normal.hpp>
//#include <boost/math/distributions/binomial.hpp>
#include "binomial_distribution.h"
#include "random.h"

//using boost::math::normal;
//----------------------------------------------------------------------
// Main
//----------------------------------------------------------------------
int main(int argc, char *argv[]) {
	

	for(int i = 0; i < 200; i++) {
		int trials = (int) (Random() * 1000);
		int successes = (int) (Random() * trials);
		double alpha = Random() * 1e-5;
		//double boost_lower = boost::math::binomial::find_lower_bound_on_p(trials, successes, alpha);
		//double boost_upper = boost::math::binomial::find_upper_bound_on_p(trials, successes, alpha);
		double ac_lower = binom_distrib_find_lower_bound_on_p(trials, successes, alpha);
		double ac_upper = binom_distrib_find_upper_bound_on_p(trials, successes, alpha);

		//printf("trials=%d, successes=%d, alpha=%.4lf, boost-lower=%.5lf, boost-upper=%.5lf\n", 
		//		trials, successes, alpha, boost_lower, boost_upper
		//		);
		printf("trials=%d, successes=%d, alpha=%.4lf,    ac-lower=%.5lf,    ac-upper=%.5lf\n", 
				trials, successes, alpha, ac_lower, ac_upper
				);
	}
				



	return(0);
}

