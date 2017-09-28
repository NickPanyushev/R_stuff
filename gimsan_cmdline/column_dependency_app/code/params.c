#include "params.h"

void printUsage() {
	printf("usage: column_dependency.out [arguments]\n");
	printf("Maximum span allowed: %d\n\n", SPAN_CAPACITY - 1);
	printf("Required arguments:\n");
	printf("-fsa <file>       FASTA file\n");
	printf("-loc <file>       Motif location file (default: each sequence has a site at position 1)\n");
	printf("                  Sequence index ranges from 0 to n-1, where n is the number of sequences.\n");
	printf("                  Location index ranges from 1 to m or from -1 to -m, where m is the sequence length.\n");
	printf("\n");
	printf("Optional arguments:\n");
	printf("-l <int>          Motif span (default: sequence length)\n");
	printf("-s <int>          Set random seed (default: different up to seconds)\n");
	printf("-ic-min <flt>     Information-content lower bound. Range: [0.0, 2.0] (default: 0.25)\n");
	printf("-ic-max <flt>     Information-content upper bound. Range: [0.0, 2.0] (default: 1.75)\n");
	printf("-alpha <flt>      Significance level (default: 0.05)\n");
	printf("-beta <flt>       Binomial confidence interval level (default: 0.95)\n");
	printf("\n");
	printf("-alpha_mode       Stopping criterion using alpha threshold\n");
	printf("-stdev_mode       Stopping criterion using standard deviation threshold\n");
	printf("\n");

	exit(1);
}

static
void printParams(FILE *fptr, Params *params, int argc, char *argv[]) {
	
	for(int i = 0; i < argc; i++) {
		fprintf(fptr, "%s ", argv[i]);
	}
	fprintf(fptr, "\n");

	fprintf(fptr, "FASTA Filename: %s\n", params->fastaFilename );
	fprintf(fptr, "Number of sequences in FASTA: %d\n", params->fasta->numseqs);
	fprintf(fptr, "Minimum sequence length in FASTA: %d\n", params->fasta->minseqlen);
	fprintf(fptr, "Maximum sequence length in FASTA: %d\n", params->fasta->maxseqlen);
	fprintf(fptr, "\n");
	fprintf(fptr, "SiteLoc Filename: %s\n", params->sitelocFilename );
	fprintf(fptr, "Number of sites: %d\n", params->siteloc->numsites);
	fprintf(fptr, "\n");
	fprintf(fptr, "Random seed: %u\n", params->randomSeed);
	fprintf(fptr, "Motif span: %d\n", params->motifspan);
	fprintf(fptr, "Alpha (significance level): %lf\n", params->alpha);
	fprintf(fptr, "Beta (binomial confidence interval level): %lf\n", params->beta);
	if(params->stopMode == ALPHA_MODE) {
		fprintf(fptr, "Stopping mode: alpha-mode");
	}
	else if(params->stopMode == STDEV_MODE) {
		fprintf(fptr, "Stopping mode: stdev-mode");
	}
	else {
		abort();
	}
	fprintf(fptr, "\n");
	fprintf(fptr, "Initial number of permutations per column pair: %d\n", params->initNumPerms);
	fprintf(fptr, "Max number of permutations per column pair: %d\n", params->capacityNumPerms);
	fprintf(fptr, "\n");
	fprintf(fptr, "Column information content lower bound: %0.2lf\n", params->lowerColInfoContent);
	fprintf(fptr, "Column information content upper bound: %0.2lf\n", params->upperColInfoContent);
	fprintf(fptr, "\n");

	//TODO: add all parameters value and print out siteloc input
}


void setColDependParams(Params *params, ColDepend *coldep) {
	coldep->initNumPerms = params->initNumPerms;
	coldep->capacityNumPerms = params->capacityNumPerms; //maximum number of possible permutations allowed

	coldep->motifspan = params->motifspan;
	coldep->alpha = params->alpha;
	coldep->beta = params->beta;
	coldep->stopMode = params->stopMode;
	coldep->lowerColInfoContent = params->lowerColInfoContent;
	coldep->upperColInfoContent = params->upperColInfoContent;

	//set k-mer set
	coldep->kmerset = constructIntMat(params->siteloc->numsites, params->motifspan);	
	for(int i = 0; i < params->siteloc->numsites; i++) {
		int seqind = params->siteloc->seqInd[i];
		int concatpos = getDoubleStrand2ConcatPos(params->fasta, seqind, params->siteloc->dsPos[i]);
		if(seqind < 0 || seqind >= params->fasta->numseqs) {
			fprintf(stderr, "Error: invalid seqind found.\n");
			fprintf(stderr, "Sequence index %d when FASTA only has %d number of seqs.\n", seqind, params->fasta->numseqs);
			abort();
		}
		if(isValidConcatPos(params->fasta, seqind, concatpos, params->motifspan)) {
			for(int j = 0; j < params->motifspan; j++) {
				coldep->kmerset->vals[i][j] = params->fasta->seqs[seqind][concatpos + j];
			}
		}
		else {
			fprintf(stderr, "Error: invalid motif site at seq %d, pos %d (concatpos %d).\n", 
				seqind, params->siteloc->dsPos[i], concatpos);
			abort();
		}
	}
	
	//set countmat and freqmat
	coldep->countmat = constructIntMat(params->motifspan, NUMALPHAS);
	coldep->freqmat = constructDblMat(params->motifspan, NUMALPHAS);
	for(int i = 0; i < coldep->countmat->numrows; i++) {
		for(int j = 0; j < coldep->countmat->numcols; j++) {
			coldep->countmat->vals[i][j] = 0;
		}
	}
	for(int i = 0; i < params->siteloc->numsites; i++) {
		for(int j = 0; j < params->motifspan; j++) {
			coldep->countmat->vals[j][coldep->kmerset->vals[i][j]]++;
		}
	}
	for(int i = 0; i < coldep->countmat->numrows; i++) {
		double sum = 0.0;
		for(int j = 0; j < coldep->countmat->numcols; j++) {
			sum += coldep->countmat->vals[i][j];
		}
		for(int j = 0; j < coldep->countmat->numcols; j++) {
			coldep->freqmat->vals[i][j] = coldep->countmat->vals[i][j] / ((double)sum);
		}
	}
	if(DEBUG0) {
		for(int i = 0; i < coldep->countmat->numrows; i++) {
			double sum = 0.0;
			for(int j = 0; j < coldep->countmat->numcols; j++) {
				sum += coldep->countmat->vals[i][j];
			}
			assert(fabs(sum - params->siteloc->numsites) < 0.000001);
		}
	}

	coldep->pairs = NULL;
	coldep->numPairs = -1;
}


//open FASTA and set motif span
void openFastaHelperAndSetMotifSpan(Params *params) {
	if(params->motifspan < SPAN_MIN_CAPACITY) { //if it is undefined by user
		//this code-block determines the motif span
		Dataset *fsa = openBackgroundData(params->fastaFilename, false); //don't use revcompl here

		//check whether sequence lengths are equal
		if(fsa->minseqlen != fsa->maxseqlen) {
			fprintf(stderr, "Error: sequence lengths are not equal and motif span is not specified.\n");
			fprintf(stderr, "Filename: %s\n", params->fastaFilename);
			abort();
		}

		//set motif span to min/max seqlen
		params->motifspan = fsa->minseqlen;
		nilDataset(fsa); //destruct it (it's not very efficient but it's clean)
	}
	if(DEBUG1) {
		fprintf(stderr, "Motif span determined at helper function: %d\n", params->motifspan);
	}

	params->fasta = openDataset(params->fastaFilename, true, params->motifspan, params->motifspan); //use revcompl here
}

void openSitelocHelper(Params *params) {
	if(params->sitelocFilename == NULL) {
		params->siteloc = constructDegenerateSiteloc(params->fasta->numseqs, 1);
	}
	else {
		params->siteloc = openSiteloc(params->sitelocFilename);
	}
}



void setParamsFromProgArg(Params *params, int argc, char *argv[]) {
	if (argc <= 1) {
		printUsage();
	}

	int i = 1;
	while(i < argc) {
		if (!strcmp(argv[i],"-s")) {
			i++;
			int err = sscanf(argv[i], "%u", &(params->randomSeed)); 
			if(err<1) printUsage();
		}
		else if (!strcmp(argv[i],"-fsa")) {
			i++;
			params->fastaFilename = argv[i];
		}
		else if (!strcmp(argv[i],"-loc")) {
			i++;
			params->sitelocFilename = argv[i];
		}
		else if (!strcmp(argv[i],"-l")) {
			i++;
			int err = sscanf(argv[i], "%d", &(params->motifspan));
			if(err<1) printUsage();
		}
		else if (!strcmp(argv[i],"-ic-min")) {
			i++;
			int err = sscanf(argv[i], "%lf", &(params->lowerColInfoContent));
			if(err<1) printUsage();
		}
		else if (!strcmp(argv[i],"-ic-max")) {
			i++;
			int err = sscanf(argv[i], "%lf", &(params->upperColInfoContent));
			if(err<1) printUsage();
		}
		else if (!strcmp(argv[i],"-alpha")) {
			i++;
			int err = sscanf(argv[i], "%lf", &(params->alpha));
			if(err<1) printUsage();
		}
		else if (!strcmp(argv[i],"-beta")) {
			i++;
			int err = sscanf(argv[i], "%lf", &(params->beta));
			if(err<1) printUsage();
		}
		else if (!strcmp(argv[i],"-alpha_mode")) {
			params->stopMode = ALPHA_MODE;
		}
		else if (!strcmp(argv[i],"-stdev_mode")) {
			params->stopMode = STDEV_MODE;
		}
		else {
			printf("Unknown command: %s\n", argv[i]);
			printUsage();
		}
		i++;
	}

}
void setDefaultParams(Params *params) {
	params->fastaFilename = NULL;
	params->fasta = NULL;

	params->sitelocFilename = NULL;
	params->siteloc = NULL;

	params->randomSeed = (unsigned int)time(NULL); //random up to seconds

	params->fptr = stdout;

	//ColDepend params
	params->initNumPerms = 1000;
	params->capacityNumPerms = 1000000; //maximum number of possible permutations allowed
	params->motifspan = -1; //later reset to be the min/max sequence length
	params->alpha = 0.05;
	params->beta = 0.95;
	params->lowerColInfoContent = 0.25;
	params->upperColInfoContent = 1.75;
	params->stopMode = ALPHA_MODE;
}

void boundCheckParams(Params *params) {
	if(params->fastaFilename == NULL) {
		fprintf(stderr, "Error: FASTA must be provided.\n");
		abort();
	}
	else if(params->lowerColInfoContent >= params->upperColInfoContent 
		|| params->lowerColInfoContent < 0.0
		|| params->upperColInfoContent > 2.0) 
	{
		fprintf(stderr, "Error: inconsistent information content. Range must be from 0.0 to 2.0.\n");
		abort();
	}
	else if(params->motifspan < SPAN_MIN_CAPACITY || params->motifspan >= SPAN_CAPACITY) {
		fprintf(stderr, "Error: span must be in range [%d,%d]\n", SPAN_MIN_CAPACITY, SPAN_CAPACITY-1);
		abort();
	}
	//else if(params->alpha <= 0.0 || params->alpha >= 1.0) {
	//	fprintf(stderr, "Error: alpha must be in range (0.0, 1.0)\n");
	//	abort();
	//}
	else if(params->beta <= 0.0 || params->beta >= 1.0) {
		fprintf(stderr, "Error: beta must be in range (0.0, 1.0)\n");
		abort();
	}
}

Params* constructParams(int argc, char *argv[]) {
	Params *params = (Params*) malloc(sizeof(Params));
	setDefaultParams(params);
	setParamsFromProgArg(params, argc, argv);

	openFastaHelperAndSetMotifSpan(params);
	openSitelocHelper(params);

	boundCheckParams(params);
	printParams(params->fptr, params, argc, argv);

	//use parameters
	if(params->randomSeed == 0) {
		params->randomSeed = 1; //cannot start with seed 0
	}
	sRandom(params->randomSeed);
	if(DEBUG1) {
		fprintf(stderr, "Random seed: %u\n", params->randomSeed);
	}

	return params;

}

void destructParams(Params *params) {
	destructSiteloc(params->siteloc);
	nilDataset(params->fasta);
	free(params);
}

