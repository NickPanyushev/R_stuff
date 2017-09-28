/***
*
* This file contains the core functions for the Gibbs-sampling algorithm. 
*
***/

#ifndef _ENTSAMP_H
#define _ENTSAMP_H

#include "gibbs_util.h"


enum Sampler {INVALID_SAMP, CLR_SAMP, ODDSRATIO_SAMP, ILR_SAMP};
enum RunThrshldTyp {INVALID_RUN_THRSHLD_TYP, NUMRUNS, CPUTIME};

//------------------------------------------------------------
// Structs
//------------------------------------------------------------
typedef struct {
	double **scoremat; //lookup table [numalphas][0:maxcount]
	int numalphas;
	int maxcount;
} EntropySampTable; //lookup table for entropy sampling

//Struct for bookkeeping
typedef struct {
	char *fastafile; //filename of the sequence file
	int span; //given span of the motif (-l param)
	int numalphas; //set to 4 for nucleotides {A,C,G,T}

	//number of runs 
	enum RunThrshldTyp runThrshldTyp;
	int numruns; // number of sampling runs (-t param)
	double cpuSecThrshld; //stopping criteria in terms of cpu time in seconds
	clock_t cpuClockStart; //starting position of clock ticks
	clock_t cpuClockPreprocessEnd;
	clock_t cpuClockSamplingEnd;

	//Other options
	int trialIters; //number of trial updates iteration (less costly)
	boolean useRapidConv;
	int iterPlateauLen; // rapid convergence (-L param)
	int numFixIters; //fixed number of iterations (-F param)
	unsigned int randseed; //random seed

	bool useRevcompl;

	//Background
	//pseudoweight is in [0,1]; For entropy-sampling, 
	//this is only used when creating profile before EM.
	double pseudoweight; 

	//frequency of phase shifting [0,1]
	//setting phase shift to 1.0 will have column shift in every iteration
	double phaseShiftFreq; 

	//dim = [max-length of sequences]
	//double *posScore; 
	double *cumsumArr;

	//combination of isValid and bgscore
	double ***validBgscore;

	enum Sampler samptyp;

	EntropySampTable *entsampTable;
	Dataset *data; 
	RunSet *runset; //stores results of all runs
	FILE *fptr;

	bool isZoopsMode;
	double zoopsPseudoweight;
	Zoops *zoops;

	Markov *markov;
	int markovOrder;
	enum BackgroundModel bgmodel;
	enum BackgroundDimType bgdimtyp;

	//wndbgm option
	char *markovModelsFilename;
	int numModelsForWndbgm;
	int bgWndSize; 

	char *bgfilename;
	bool useBgfileOption;
	char *pspfilename;

	int numOfTopMotifs; //number of top motifs to display (may overlap)
	double topMotifsPvalCutoff;

	//trackers - for output only
	int totalIters; //total number of iterations for all sampling runs
	bool printScorePerRun;
} Gibbs;

//---------------------------------------
// Functions
//---------------------------------------

extern void setDefaultGibbsParams(Gibbs *gibbs);
extern void printGibbsParams(Gibbs *gibbs, int argc, char **argv);
extern boolean boundCheckGibbsParams(Gibbs *gibbs, int argc, char** argv);

//Main engine of Gibbs-sampling; returns the run-node of the best run
//extern void runGibbs(Gibbs *gibbs);
extern void runIters(RunNode *rnode, int initspan, Gibbs *gibbs);


extern void initGibbsStructs(Gibbs *gibbs);
extern void nilGibbs(Gibbs *gibbs );


#endif

