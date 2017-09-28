/***
*
* This file contains the core functions for the Gibbs-sampling algorithm. 
*
***/

#ifndef _ENTSAMP_ILR_H
#define _ENTSAMP_ILR_H

#include "entsamp.h"
#include "em_alg.h"

//------------------------------------------------------------
// Structs
//------------------------------------------------------------

//Struct for bookkeeping
typedef struct {
	Gibbs *gibbs;

	int emStep; //number of EM steps

	//metric to choose the best run among all Gibbs-runs
	enum ScoreMetric metric; 
} GibbsIlr;

//---------------------------------------
// Functions
//---------------------------------------

//Main engine of GibbsIlr; returns the run-node of the best run
extern RunNode* runGibbsIlr(GibbsIlr *gilr);

//Deallocate the struct for bookkeeping and its children
extern void nilGibbsIlr(GibbsIlr *gilr );

#endif

