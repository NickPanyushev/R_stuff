#ifndef _STDINC_H
#define _STDINC_H

#ifdef WIN32
	#include "stdafx.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <time.h>

#if DEBUG
	#define DEBUG0 1 //sanity check
#else
	#define DEBUG0 0
#endif

#if VERBOSE
	#define DEBUG1 1 //verbose mode
#else
	#define DEBUG1 0
#endif

#define DEBUG2 0

#define TRUE true
#define FALSE false

#define NUMALPHAS 4 
#define NUMCHARS 5 //includes extra gap characters
#define GAP_CHAR 4
#define SPAN_MIN_CAPACITY 2
#define SPAN_CAPACITY 1000
#define LOGZERO -1e100

typedef bool boolean;

#ifdef WIN32
	// Add definitions of NAN, NINF, and PINF
	extern double NAN;
	extern double NINF;
	extern double PINF;
#else
	//INFINITY should not be used for compatability with WIN32
	#define NINF -INFINITY
	#define PINF INFINITY
	#define _snprintf snprintf 
#endif

#endif

