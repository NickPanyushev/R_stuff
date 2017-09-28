#ifndef _PROFILE_H
#define _PROFILE_H

#include "stdinc.h"

//#define getind(profile, pos) (pos)
//#define isgap(profile, pos) ((profile)->isgap[(pos)])

typedef struct {
	int cols; //number of motif columns (gapless)
	int span; //max len including gaps
	int maxspan;

	double **mat;

	bool *isgap; //length maxspan; default is all zero 
} Profile; 


extern Profile *initProfile(int initspan, int maxspan);
extern void nilProfile(Profile *profile);

//no gaps and set all value of pswm to zero
extern void resetProfile(Profile *profile, int initspan);

//keep gaps configuration of "profile" but set all values of profile to val
extern void setProfile(Profile *profile, double val);

//in-place reverse-complement (this will modify the profile)
extern void revcomplProfile(Profile *profile);

extern void copyProfile(Profile *dest, Profile *src);

extern void addFrontCol(Profile *profile, double *initval);
extern void addBackCol(Profile *profile, double *initval);
extern void addProfileCol(Profile *profile, int pos, double *initval);
extern void rmProfileCol(Profile *profile, int pos);

extern void shiftProfileLeft(Profile *profile, double *initval);
extern void shiftProfileRight(Profile *profile, double *initval);

extern void printProfile(FILE *fptr, Profile *profile);
extern void printCountmat(FILE *fptr, Profile *countmat);

#endif
