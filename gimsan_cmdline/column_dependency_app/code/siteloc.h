#ifndef _SITELOC_H
#define _SITELOC_H

#include "stdinc.h"

typedef struct {
	int numsites;
	int *seqInd;
	int *dsPos;
} Siteloc; 

extern Siteloc* constructDegenerateSiteloc(int numsites, int default_dspos);
extern Siteloc* openSiteloc(char *filename);
extern void destructSiteloc(Siteloc*);



#endif

