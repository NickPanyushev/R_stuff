
#ifndef _RANDOM_H
#define _RANDOM_H

#include "mt19937ar.h"

extern void sRandom(unsigned int seed);

//extern double Random();
inline double Random() { return genrand_real2(); }

#endif

