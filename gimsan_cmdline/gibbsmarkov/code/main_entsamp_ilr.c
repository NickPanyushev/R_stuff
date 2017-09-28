/**
*
* This file contains the main() function for EntropySamplerILR. It contains functions to
* parse arguments from the executable. It is also the engine to run the
* core functions of the motif finder.
*
**/

#include "entsamp_ilr.h"
#include "print_results.h"

//----------------------------------------------------------------------
// Print-screen parameters and errors
//----------------------------------------------------------------------

static
void printerr() {
	printf("usage: gibbsmarkov <filename> -l <span> [options]\n");
	printf("Maximum span allowed: %d\n\n", SPAN_CAPACITY - 1);
	printf("Options:\n\n");
	printf("-s <int>           set random seed (default: different up to seconds)\n");
	printf("-t <int>           number of runs\n");
	printf("-cput <int>        cpu time as stopping criterion (in seconds)\n");
	printf("-L <int>           set rapid convergence limit\n");
	printf("-F <int>           fixed number of iterations\n");
	printf("\n");
	printf("-p <flt>           pseudo-weight [-gibbsamp only] (default: 0.05)\n");
	//printf("-opt_pw          choose best based on MCLR/MILR with pseudoweight [-gibbsamp only]\n");
	//printf("                 For -best_clr, pseudoweight is used for both best iters and runs.\n");
	//printf("                 For -best_ilr, pseudoweight is used for best iter. Then pseudoweight\n");
	//printf("                 is added to the PWM at each iteration of M-step of EM. \n");
	printf("\n");
	//printf("-r <int>           print multiple best motifs (default: 1)\n");
	printf("-ds                consider double-strand, i.e. reverse-complement (default: off)\n");
	printf("-ps <flt>          phase shift frequency (default: 0.4)\n");
	//printf("-em <int>          number of iterations to run EM (default: 0)\n");
	printf("-zoops <flt>       Bayesian ZOOPS mode with pseudoweight (recommend: 0.2)\n");
	printf("-print_runs        print scores of each run\n");
	printf("\n");
	//printf("Sampling technique (default: -gibbsamp):\n");
	//printf("   -gibbsamp       use odds-ratio (original gibbs) to sample\n");
	//printf("   -clrsamp        use CLR/entropy to sample\n");
	////printf("   -ilrsamp      use ILR to sample\n");
	//printf("\n");
	printf("Scoring scheme (default: -best_clr):\n");
	printf("   -best_clr       choose best by CLR/entropy score\n");
	printf("   -best_ilr       choose best by ILR score\n");
	printf("\n");
	printf("-markov <int>      order of Markov background (default: 3)\n");
	printf("-bfile <name>      fasta file to estimate background probabilities\n");
	printf("-psp <name>        filename of position-specific prior\n");
	//printf("-seqspecbgm        sequence-specific background model (same as -bgdim_numseqs)\n");
	//printf("\n");
	//printf("-wndspecbgm <int>  window-specific background model with window-size\n");
	//printf("-markovfile <name> filename that contains the Markov chain for wndspecbgm\n");
	printf("\n");
	printf("Background model (default: bg_gm):\n");
	printf("   -bg_gm           use GibbsMarkov model\n");
	printf("   -gmean_strands   use geometric mean of the two strands\n");
	//printf("   -bg_biopro    use BioProspector model\n");
	//printf("   -bg_ms        use MotifSampler model\n");

	exit(1);
}

static
void printparams(GibbsIlr *gilr, int argc, char **argv) {
	Gibbs *gibbs = gilr->gibbs;

	for(int i = 0; i < argc; i++) {
		fprintf(gibbs->fptr, "%s ", argv[i]);
	}
	fprintf(gibbs->fptr, "\n\n");

	fprintf(gibbs->fptr, "Filename: %s\n", gibbs->fastafile);
	fprintf(gibbs->fptr, "Number of sequences: %d\n", gibbs->data->numseqs);
	fprintf(gibbs->fptr, "Total residues: %d\n", (gibbs->useRevcompl ? gibbs->data->totalCount / 2 : gibbs->data->totalCount));
	fprintf(gibbs->fptr, "Average sequence length: %.1lf\n", (gibbs->useRevcompl ? gibbs->data->avgseqlen/ 2 : gibbs->data->avgseqlen));
	fprintf(gibbs->fptr, "Minimum sequence length: %d\n", (gibbs->useRevcompl ? gibbs->data->minseqlen/ 2 : gibbs->data->minseqlen));
	fprintf(gibbs->fptr, "Maximum sequence length: %d\n", (gibbs->useRevcompl ? gibbs->data->maxseqlen / 2 : gibbs->data->maxseqlen));
	fprintf(gibbs->fptr, "Include reverse-complement: %s\n", (gibbs->useRevcompl ? "yes" : "no") );
	fprintf(gibbs->fptr, "Span: %d\n", gibbs->span);
	fprintf(gibbs->fptr, "Random seed: %u\n", gibbs->randseed);
	fprintf(gibbs->fptr, "\n");

	fprintf(gibbs->fptr, "Number of motifs to be printed: %d\n", gibbs->numOfTopMotifs);
	if(gibbs->numOfTopMotifs > 1) {
		fprintf(gibbs->fptr, "Top motifs p-value cutoff: %lf\n", gibbs->topMotifsPvalCutoff);
	}
	fprintf(gibbs->fptr, "\n");
	if(gibbs->runThrshldTyp == NUMRUNS) {
		fprintf(gibbs->fptr, "Stopping criterion: Fixed number of runs\n");
		fprintf(gibbs->fptr, "Number of runs: %d\n", gibbs->numruns);
	}
	else {
		fprintf(gibbs->fptr, "Stopping criterion: CPU-time\n");
		fprintf(gibbs->fptr, "CPU time in seconds: %d\n", (int)gibbs->cpuSecThrshld);
	}

	fprintf(gibbs->fptr, "\n");
	fprintf(gibbs->fptr, "Rapid convergence: %s\n", (gibbs->useRapidConv ? "on" : "off"));
	if(gibbs->useRapidConv) {
		fprintf(gibbs->fptr, "Iterations plateau (-L): %d\n", gibbs->iterPlateauLen);
	}
	else {
		fprintf(gibbs->fptr, "Fixed iterations (-F): %d\n", gibbs->numFixIters);
	}

	fprintf(gibbs->fptr, "Number of trial iterations: %d\n", gibbs->trialIters);
	fprintf(gibbs->fptr, "Phase-shift frequency (-ps): %.2lf\n", gibbs->phaseShiftFreq);
	fprintf(gibbs->fptr, "\n");

	fprintf(gibbs->fptr, "Occurrences per sequence: %s\n", (gibbs->isZoopsMode ? "Bayesian ZOOPS" : "OOPS"));
	if(gibbs->isZoopsMode) {
		fprintf(gibbs->fptr, "Bayesian ZOOPS pseudoweight: %.2lf\n", gibbs->zoopsPseudoweight);
	}
	fprintf(gibbs->fptr, "\n");

	fprintf(gibbs->fptr, "Pseudo-weight (-p): %.2lf\n", gibbs->pseudoweight);
	if(gibbs->samptyp == CLR_SAMP) {
		fprintf(gibbs->fptr, "Sampling Technique: CLR Sampler\n");
	}
	else if(gibbs->samptyp == ODDSRATIO_SAMP){
		fprintf(gibbs->fptr, "Sampling Technique: Gibbs Sampler\n");
	}
	else if(gibbs->samptyp == ILR_SAMP){
		fprintf(gibbs->fptr, "Sampling Technique: ILR Sampler\n");
	}
	else {
		exit(1);
	}

	fprintf(gibbs->fptr, "\n");
	fprintf(gibbs->fptr, "Order of Markov background (-markov): %d\n", gibbs->markovOrder);
	if(gibbs->bgmodel == BG_GIBBSMARKOV) {
		fprintf(gibbs->fptr, "Markov background model: %s\n", "GibbsMarkov");
	}
	else if(gibbs->bgmodel == BG_BIOPRO) {
		fprintf(gibbs->fptr, "Markov background model: %s\n", "BioProspector");
	}
	else if(gibbs->bgmodel == BG_MOTIFSAMPLER) {
		fprintf(gibbs->fptr, "Markov background model: %s\n", "MotifSampler");
	}
	else if(gibbs->bgmodel == BG_GMEAN) {
		fprintf(gibbs->fptr, "Markov background model: %s\n", "Two-strand geometric mean");
	}
	else if(gibbs->bgmodel == BG_AMEAN) {
		fprintf(gibbs->fptr, "Markov background model: %s\n", "Two-strand arithmetic mean");
	}


	if(gibbs->bgdimtyp == BGDIM_AGGREGATE) {
		fprintf(gibbs->fptr, "Background dimension: Aggregate\n");
		fprintf(gibbs->fptr, "Background file (-bfile): %s\n", gibbs->bgfilename);
	}
	else if(gibbs->bgdimtyp == BGDIM_NUMSEQS) {
		fprintf(gibbs->fptr, "Background dimension: Number of sequences\n");
		fprintf(gibbs->fptr, "Background file (-bfile): %s\n", gibbs->bgfilename);
	}
	else if(gibbs->bgdimtyp == BGDIM_WND) {
		fprintf(gibbs->fptr, "Background dimension: Windows\n");
		if(gibbs->markovModelsFilename != NULL) {
			fprintf(gibbs->fptr, "Markov chain file: %s\n", gibbs->markovModelsFilename);
		}
		fprintf(gibbs->fptr, "Number of Markov models for -wndspecbgm: %d\n", gibbs->numModelsForWndbgm);
		fprintf(gibbs->fptr, "Background window size: %d\n", gibbs->bgWndSize);
	}

	if(gibbs->pspfilename == NULL) {
		fprintf(gibbs->fptr, "Position-specific prior file (-psp): none\n");
	}
	else {
		fprintf(gibbs->fptr, "Position-specific prior file (-psp): %s\n", gibbs->pspfilename);
	}
	fprintf(gibbs->fptr, "\n");

	fprintf(gibbs->fptr, "Scoring scheme: %s\n", scoreMetricToStr(gilr->metric));
	fprintf(gibbs->fptr, "EM iterations (-em): %d\n", gilr->emStep);
	fprintf(gibbs->fptr, "\n");

	//input FASTA freq
	int total = 0;
	for(int i = 0; i < NUMALPHAS; i++) {
		total+= gibbs->data->count[i];
	}
	for(int i = 0; i < NUMALPHAS; i++) {
		fprintf(gibbs->fptr, "Input FASTA freq of %c: %.4lf\n", numToChar(i), ((double)gibbs->data->count[i])/total);
	}
	fprintf(gibbs->fptr, "\n");

	//background model
	for(int i = 0; i < NUMALPHAS; i++) {
		fprintf(gibbs->fptr, "Background model freq of %c: %.4lf\n", numToChar(i), gibbs->markov->bgmodelFreq[i]);
	}
	fprintf(gibbs->fptr, "\n");

	//print pseudocount
	//for(int i = 0; i < NUMALPHAS; i++) {
	//	fprintf(gibbs->fptr, "Pseudocount of %c: %.4lf\n", numToChar(i),
	//		gibbs->pseudocount[i]);
	//}
	fprintf(gibbs->fptr, "\n");

}


//----------------------------------------------------------------------
// Params and Initialization
//----------------------------------------------------------------------


void setDefaultGibbsParams(GibbsIlr *gilr) {
	Gibbs *gibbs = gilr->gibbs;

	//defaults
	gibbs->numalphas = NUMALPHAS;
	gibbs->span = -1;
	gibbs->randseed = (unsigned int)time(NULL); //random up to seconds

	//runs
	gibbs->runThrshldTyp = INVALID_RUN_THRSHLD_TYP;
	gibbs->cpuSecThrshld = -1.0;
	gibbs->numruns = -1;

	gibbs->trialIters = 5; //no parameters to set this explicitly, but -ilrsamp is set differently

	//iters
	gibbs->useRapidConv = TRUE;
	gibbs->iterPlateauLen = -1;
	gibbs->numFixIters = -1;

	gibbs->pseudoweight = 0.05;

	gibbs->phaseShiftFreq = 0.4;
	gibbs->markovOrder = 3;
	gibbs->useRevcompl = false;

	gibbs->numOfTopMotifs = 1;
	gibbs->topMotifsPvalCutoff = 0.001;

	//gibbs->samptyp = INVALID_SAMP;
	gibbs->samptyp = ODDSRATIO_SAMP; //-gibbsamp
	gibbs->bgmodel = BG_GIBBSMARKOV;
	gibbs->bgdimtyp = BGDIM_AGGREGATE;
	gibbs->bgWndSize = -1;

	gibbs->fptr = stdout;
	gibbs->printScorePerRun = false;

	gibbs->useBgfileOption = false;
	gibbs->bgfilename = NULL;
	gibbs->pspfilename = NULL;

	gibbs->markovModelsFilename = NULL;
	gibbs->numModelsForWndbgm = 0;

	//ZOOPS
	gibbs->isZoopsMode = false;
	gibbs->zoopsPseudoweight = NAN;

	//gilr
	gilr->emStep = 0;
	//gilr->metric = NO_SCORE;
	gilr->metric = CLR; //-best_clr

	//tracker - for output only
	gibbs->totalIters = 0; //no params to set this
}

static
boolean boundCheckGibbsParams(GibbsIlr *gilr, int argc, char** argv) {
	Gibbs *gibbs = gilr->gibbs;

	if(gibbs->span < 2 || gibbs->span >= SPAN_CAPACITY) {
		printf("Span out of bound\n");
		return FALSE;
	}
	if(gibbs->iterPlateauLen <= 0 && gibbs->numFixIters <= 0) {
		printf("Number of iterations (-L or -F) not specified.\n");
		return FALSE;
	}
	if(gibbs->iterPlateauLen >= 0 && gibbs->numFixIters >= 0) {
		printf("Error: -L and -F cannot be used together.\n");
		return FALSE;
	}

	if(gibbs->runThrshldTyp == NUMRUNS && gibbs->numruns <= 0 ) {
		printf("Invalid number of runs.\n");
		return FALSE;
	}

	if(	gibbs->runThrshldTyp == INVALID_RUN_THRSHLD_TYP) {
		printf("Error: the options -t or -cput must be specified.\n");
		return FALSE;
	}
	if(gibbs->cpuSecThrshld > 0.0 && gibbs->numruns > 0) {
		printf("Error: the options -t and -cput cannot be used together.\n");
		return FALSE;
	}

	if(gibbs->pseudoweight < 0.0 || gibbs->pseudoweight > 1.0) {
		printf("Pseudo weight out of bound\n");
		return FALSE;
	}
	if(gibbs->phaseShiftFreq < 0.0 || gibbs->phaseShiftFreq > 1.0) {
		printf("Column shift frequency out of bound\n");
		return FALSE;
	}
	if(gibbs->samptyp == INVALID_SAMP) {
		printf("Sampling Technique must be specified.\n");
		return FALSE;
	}
	if(gibbs->markovOrder < 0 || gibbs->markovOrder > MARKOV_ORDER_BOUND) {
		printf("The order of Markov background exceeds bound.\n");
		return FALSE;
	}
	if(gibbs->numOfTopMotifs <= 0) {
		printf("Invalid number of motifs to be printed.\n");
		return FALSE;
	}

	//if(gibbs->markovModelsFilename != NULL && gibbs->bgdimtyp != BGDIM_WND) {
	//	printf("Markovfile can only be specified with -wndspecbgm.\n");
	//	return FALSE;
	//}

	if(gibbs->isZoopsMode && gibbs->zoopsPseudoweight < 0.0) {
		printf("ZOOPS pseudoweight must be greater than 0.0\n");
		return FALSE;
	}
	if(gibbs->isZoopsMode && gilr->emStep > 0) {
		printf("Error: EM for ZOOPS mode has not been implemented.\n");
		return FALSE;
	}
	//if(gilr->emStep > 0) {
	//	printf("EM algorithm is deprecated.\n");
	//	return FALSE;
	//}

	//checking whether -F and -L are used at the same time
	boolean seenL = FALSE;
	boolean seenF = FALSE;
	for(int i = 0; i < argc; i++) {
		if(!strcmp(argv[i],"-L")) {
			seenL = TRUE;
		}
		if(!strcmp(argv[i],"-F")) {
			seenF = TRUE;
		}
	}
	if(seenL && seenF) {
		printf("Fixed iterations and rapid convergence cannot be used together.\n");
		return FALSE;
	}

	if( gibbs->samptyp == CLR_SAMP) {
		for(int i = 0; i < argc; i++) {
			if(!strcmp(argv[i],"-p")) {
				printf("Pseudocount cannot used when using entropy to sample\n");
				return FALSE;
			}
		}
	}
	if(gilr->metric == NO_SCORE) {
		printf("Scoring scheme must be specified.\n");
		return FALSE;
	}
	if(gilr->emStep < 0) {
		printf("Number of EM step must be >= 0\n");
		return FALSE;
	}
	if(gibbs->bgdimtyp == BGDIM_WND) {
		if(gibbs->useBgfileOption) {
			printf("Cannot use -bfile and -wnd_bgm at the same time\n");
			return FALSE;
		}
	}
	if(gibbs->numOfTopMotifs != 1) {
		printf("Number of top motifs > 1 is not implemented yet.\n");
		return FALSE;
	}

	return TRUE;
}

static
GibbsIlr *parseArgs(int argc, char **argv) {
	Gibbs *gibbs = (Gibbs*) malloc(sizeof(Gibbs));
	GibbsIlr *gilr = (GibbsIlr*) malloc(sizeof(GibbsIlr));
	gilr->gibbs = gibbs;

	if (argc < 3) {
		printerr();
	}
	gibbs->fastafile = argv[1];
	setDefaultGibbsParams(gilr);

	int i = 2;
	int error;
	while(i < argc) {
		if (!strcmp(argv[i],"-s")) {
			i++;
			error = sscanf(argv[i], "%u", &(gibbs->randseed));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-l")) {
			i++;
			error = sscanf(argv[i], "%d", &(gibbs->span));
			if(error<1) printerr();
		}

		else if (!strcmp(argv[i],"-t")) {
			i++;
			gibbs->runThrshldTyp = NUMRUNS;
			error = sscanf(argv[i], "%d", &(gibbs->numruns));
			if(error < 1) printerr();
		}
		else if (!strcmp(argv[i],"-cput")) {
			i++;
			gibbs->runThrshldTyp = CPUTIME;
			int int_cpusec = 0;
			error = sscanf(argv[i], "%d", &(int_cpusec));
			gibbs->cpuSecThrshld = (double) int_cpusec;
			if(error < 1) printerr();
		}

		else if (!strcmp(argv[i],"-L")) {
			i++;
			error = sscanf(argv[i], "%d", &(gibbs->iterPlateauLen));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-F")) {
			i++;
			error = sscanf(argv[i], "%d", &(gibbs->numFixIters));
			gibbs->useRapidConv = FALSE;
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-ds")) {
			gibbs->useRevcompl = true;
		}
		else if (!strcmp(argv[i],"-ps")) {
			i++;
			error = sscanf(argv[i], "%lf", &(gibbs->phaseShiftFreq));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-r")) {
			i++;
			error = sscanf(argv[i], "%d", &(gibbs->numOfTopMotifs));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-zoops")) {
			gibbs->isZoopsMode = true;
			i++;
			error = sscanf(argv[i], "%lf", &(gibbs->zoopsPseudoweight));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-markov")) {
			i++;
			error = sscanf(argv[i], "%d", &(gibbs->markovOrder));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-bg_gm")) {
			gibbs->bgmodel = BG_GIBBSMARKOV;
		}
		else if (!strcmp(argv[i],"-gmean_strands")) {
			gibbs->bgmodel = BG_GMEAN;
		}
		else if (!strcmp(argv[i],"-amean_strands")) {
			gibbs->bgmodel = BG_AMEAN;
		}
		//else if (!strcmp(argv[i],"-bg_biopro")) {
		//	gibbs->bgmodel = BG_BIOPRO;
		//}
		//else if (!strcmp(argv[i],"-bg_ms")) {
		//	gibbs->bgmodel = BG_MOTIFSAMPLER;
		//}
		//else if (!strcmp(argv[i],"-bgdim_numseqs") || !strcmp(argv[i],"-seqspecbgm")) {
		//	gibbs->bgdimtyp = BGDIM_NUMSEQS;
		//}
		//else if (!strcmp(argv[i],"-wndspecbgm")) {
		//	gibbs->bgdimtyp = BGDIM_WND;
		//	i++;
		//	error = sscanf(argv[i], "%d", &(gibbs->bgWndSize));
		//	if(error<1) printerr();
		//}
		//else if(!strcmp(argv[i], "-markovfile")) {
		//	i++;
		//	gibbs->markovModelsFilename = argv[i];
		//}
		else if (!strcmp(argv[i],"-p")) {
			i++;
			error = sscanf(argv[i], "%lf", &(gibbs->pseudoweight));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-gibbsamp")) {
			gibbs->samptyp = ODDSRATIO_SAMP;
		}
		else if (!strcmp(argv[i],"-clrsamp") || !strcmp(argv[i],"-entsamp")) {
			gibbs->samptyp = CLR_SAMP;
		}
		//else if (!strcmp(argv[i],"-ilrsamp")) {
		//	gibbs->samptyp = ILR_SAMP;
		//}
		else if(!strcmp(argv[i], "-bfile")) {
			i++;
			gibbs->bgfilename = argv[i];
			gibbs->useBgfileOption = true;
		}
		else if(!strcmp(argv[i], "-psp")) {
			//default is NULL (when uninitialized)
			i++;
			gibbs->pspfilename = argv[i];
		}
		else if (!strcmp(argv[i],"-em")) {
			i++;
			error = sscanf(argv[i], "%d", &(gilr->emStep));
			if(error<1) printerr();
		}
		else if (!strcmp(argv[i],"-best_clr") || !strcmp(argv[i],"-best_ent")) {
			gilr->metric = CLR;
		}
		else if (!strcmp(argv[i],"-best_ilr")) {
			gilr->metric = ILR;
		}
		else if (!strcmp(argv[i],"-print_runs")) {
			gilr->gibbs->printScorePerRun = true;
		}


		else if (!strcmp(argv[i],"-n")) {
			//nothing, doesn't apply for this gibbs finder
		}
		else {
			printf("Unknown command: %s\n", argv[i]);
			printerr();
		}
		i++;
	}

	boolean isValid = boundCheckGibbsParams(gilr, argc, argv);

	if(!isValid) {
		//printerr();
		exit(1);
	}

	return gilr;
}

static
void printElapsed(Gibbs *gibbs) {
	double elapsed = double(clock())/CLOCKS_PER_SEC  - double(gibbs->cpuClockStart)/CLOCKS_PER_SEC;
	double preprocess = double(gibbs->cpuClockPreprocessEnd)/CLOCKS_PER_SEC  - double(gibbs->cpuClockStart)/CLOCKS_PER_SEC;
	double sampling = double(gibbs->cpuClockSamplingEnd)/CLOCKS_PER_SEC  - double(gibbs->cpuClockPreprocessEnd)/CLOCKS_PER_SEC;
	double postprocess = fabs(double(clock())/CLOCKS_PER_SEC  - double(gibbs->cpuClockSamplingEnd)/CLOCKS_PER_SEC);

	printf("Number of runs: %d\n", gibbs->runset->len);
	printf("Preprocessing CPU time (in seconds): %.2lf\n", preprocess);
	printf("Sampling-step CPU time (in seconds): %.2lf\n", sampling);
	printf("Postprocessing CPU time (in seconds): %.2lf\n", postprocess);
	printf("Total elapsed CPU time (in seconds): %.2lf\n", elapsed);
	printf("\n");
}

//A simple loop will be printing the scores in backward order
static
void recursion_printScorePerRun(RunNode *node) {
	if(node != NULL) {
		recursion_printScorePerRun(node->next);
		printf("run %4d: %7.2lf\n", node->runId, node->score);
	}
}

static
void printScorePerRun(RunSet *runset, enum ScoreMetric metric) {
	if(metric == ILR) {
		printf("Log Markovian-ILR of PWM (after EM) for all runs:\n");
	}
	else {
		printf("Log Markovian-CLR for all runs:\n");
	}
	recursion_printScorePerRun(runset->head);
	printf("\n");
}




//----------------------------------------------------------------------
// Main
//----------------------------------------------------------------------
int main(int argc, char *argv[]) {
	if(DEBUG0) {
		printf("WARNING: This is currently running under DEBUG mode.\n");
	}
	if(DEBUG1) {
		printf("Verbose mode.\n");
	}
	printf("Compiled on " __DATE__ " " __TIME__ "\n");
	printf("\n");

	printf("ChangeLog\n");
	printf("Version 2.2.0\n");
	printf("20100609. PSP mode\n");
	printf("20100527. Allow for multiple consecutive > in the header\n");
	printf("20081128. ZOOPS mode is allow to have sequences that are shorter than span\n");
	printf("20081115. Modified the chooseBestSite routine to use P(Y_i > 0)\n");
	printf("20080807. Code cleanup for GIMSAN\n");
	printf("20080625. Output of background model frequency\n");
	printf("20080625. Out of range sequence are counted toward char 'X'\n");
	printf("20080408. Increased span-capacity and fixed bug where span may be larger than minimum seqlen.\n");
	printf("\n");

	RunNode *bestnode;

	GibbsIlr *gilr = parseArgs(argc, argv);
	Gibbs *gibbs = gilr->gibbs;

	gilr->gibbs->cpuClockStart = clock();

	initGibbsStructs(gibbs);
	printparams(gilr, argc, argv);

	gilr->gibbs->cpuClockPreprocessEnd = clock();

	bestnode = runGibbsIlr(gilr);

	gilr->gibbs->cpuClockSamplingEnd = clock();

	if(gibbs->numOfTopMotifs == 1) {
		printRunNode(bestnode, gilr->metric, gibbs->markov, gibbs->data, gibbs->zoops);
	}
	else {
		//printTopRankRunNodes(gibbs->numOfTopMotifs, gibbs->runset, gilr->metric,
		//	gibbs->markov, gibbs->data, gibbs->zoops, gibbs->topMotifsPvalCutoff);
	}

	printf("Average number of iterations per run: %.2lf\n",
		((double)gibbs->totalIters) / gibbs->runset->len);
	printElapsed(gibbs);

	if(gibbs->printScorePerRun) {
		printScorePerRun(gibbs->runset, gilr->metric);
	}

	nilGibbsIlr(gilr);
	return(0);
}

