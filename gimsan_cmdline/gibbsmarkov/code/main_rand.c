#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include "random.h"

//----------------------------------------------------------------------
// Print-screen parameters and errors
//----------------------------------------------------------------------

static
void printerr() {
	printf("usage: rand -s <int> -n <int>\n");
	printf("\n");
	printf("Options:\n\n");
	printf("-s <unsigned long>   set random seed\n");
	printf("-n                   number of random numbers to generate\n");
	printf("\n");
	printf("ULONG_MAX = %lu\n", ULONG_MAX);
	printf("\n");

	exit(1);
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		printerr();
	}

	unsigned long randseed = (unsigned long)time(NULL); //random up to seconds;
	int num_rands = 1;

	int i = 1;
	int error;
	while(i < argc) {
		if (!strcmp(argv[i],"-s")) {
			i++;
			error = sscanf(argv[i], "%lu", &(randseed)); 
			if(error<1) { fprintf(stderr, "Error in specified randoom seed.\n"); exit(1);}
		}
		else if (!strcmp(argv[i],"-n")) {
			i++;
			error = sscanf(argv[i], "%d", &(num_rands));
			if(error<1) { fprintf(stderr, "Error in -n value\n"); exit(1);}
		}
		else {
			fprintf(stderr, "Unknown command: %s\n", argv[i]);
			printerr();
		}
		i++;
	}

	fprintf(stderr, "Random seed: %lu\n", randseed);
	sRandom(randseed);

	for(int i = 0; i < num_rands; i++) {
		printf("%lf\n", Random());
	}

	return 0;
}

