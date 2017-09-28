/***
*
* This file contains the helper functions for the Gibbs-sampling algorithm. 
*
***/

#ifndef _GIBBS_UTIL_H
#define _GIBBS_UTIL_H

#include "stdinc.h"
#include "profile.h"
#include "dataset.h"
#include "random.h"
#include "markov_score.h"

#define EMPTY_SITE -1

enum ScoreMetric { NO_SCORE, ILR, CLR };

//------------------------------------------------------------
// Structs
//------------------------------------------------------------
typedef struct RunNode_el {
	struct RunNode_el *next;
	double score; 
	int runId; //run ID
	int *sites; //set of sites (one for each position)

	//These profiles may not be updated to be consistent
	//with the current set of sites. 
	//Hence, use with caution!!
	Profile *countmat; //kept to figure out span and columns
	Profile *pswm; //position specific weight matrix
} RunNode;


//List of RunNode
typedef struct {
	RunNode *head;
	int len;
} RunSet;


//---------------------------------------
// RunNode/RunSet functions
//---------------------------------------
extern RunNode* createRunNode(int runId, Dataset *data, int initspan, int maxspan);
extern void nilRunNode(RunNode *rnode);

extern RunSet* createRunSet();
extern void nilRunSet(RunSet *rset);

//---------------------------------------
// Scoring functions
//---------------------------------------

//Convert metric to string
extern char* scoreMetricToStr(enum ScoreMetric metric);

//---------------------------------------
// Matrix/profile updates functions
//---------------------------------------
extern void updateCountmatFromSites(Profile *countmat, int *sites, Dataset *data);
extern boolean validCountmatWithSites(Profile *countmat, int *sites, Dataset *data);


//----------------------------------------------------------------------
// add/remove/set sites functions
//----------------------------------------------------------------------
//set "sites" to have a new set of random starting positions
extern void setRandomSites(int *sites, int initspan, Dataset *data, Zoops *zoops);

extern int addSite(int newsite, int seqind, int *sites, int numsites, Profile *countmat, Dataset *data);
extern int removeSite(int seqind, int *sites, int numsites, Profile *countmat, Dataset *data);

//extern void findBestSites(Markov *markov, Profile *freqmat, int *sites, Dataset *data);
extern void findBestSitesByAvgLR(Markov *markov, Profile *freqmat, int *sites, Dataset *data);

//----------------------------------------------------------------------
// phase-shift
//----------------------------------------------------------------------
extern void attemptPhaseShift(Markov *markov, Dataset*, Profile *countmat, int *sites);

//----------------------------------------------------------------------
// Combination of bgscore and isValidSite
//----------------------------------------------------------------------
extern double*** buildValidBgscore(int **isValidSite[SPAN_CAPACITY], double ***bgscore, int minspan, int maxspan, Dataset *data);
extern void nilValidBgscore(double ***validBgscore, int minspan, int maxspan, int numseqs);

#endif

