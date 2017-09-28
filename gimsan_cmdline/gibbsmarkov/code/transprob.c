#include "transprob.h"


static
bool isValidAlpha(int alpha) {
	bool flag = alpha >=0 && alpha < NUMALPHAS;
	assert(flag);
	return flag;
}

//static double getTransProb(TransProb *transprob, int order, int alpha0, int alpha1 = -1, int alpha2 = -1, 
//		 int alpha3 = -1, int alpha4 = -1); //for default values

//valid for 4th order markov chain
//There are default values in markov.h
double getTransProb(TransProb *transprob, int order, int alpha0, int alpha1, int alpha2, int alpha3, int alpha4, int alpha5) {
	if(DEBUG0) {
		//within bound
		assert(order >= 0);
		assert(order <= transprob->maxorder);

		if(order == 0) {
			isValidAlpha(alpha0);
			assert(alpha1 == -1);
			assert(alpha2 == -1);
			assert(alpha3 == -1);
			assert(alpha4 == -1);
			assert(alpha5 == -1);
		}
		else if(order == 1) {
			isValidAlpha(alpha0);
			isValidAlpha(alpha1);
			assert(alpha2 == -1);
			assert(alpha3 == -1);
			assert(alpha4 == -1);
			assert(alpha5 == -1);
		}
		else if(order == 2) {
			isValidAlpha(alpha0);
			isValidAlpha(alpha1);
			isValidAlpha(alpha2);
			assert(alpha3 == -1);
			assert(alpha4 == -1);
			assert(alpha5 == -1);
		}
		else if(order == 3) {
			isValidAlpha(alpha0);
			isValidAlpha(alpha1);
			isValidAlpha(alpha2);
			isValidAlpha(alpha3);
			assert(alpha4 == -1);
			assert(alpha5 == -1);
		}
		else if(order == 4) {
			isValidAlpha(alpha0);
			isValidAlpha(alpha1);
			isValidAlpha(alpha2);
			isValidAlpha(alpha3);
			isValidAlpha(alpha4);
			assert(alpha5 == -1);
		}
		else if(order == 5) {
			isValidAlpha(alpha0);
			isValidAlpha(alpha1);
			isValidAlpha(alpha2);
			isValidAlpha(alpha3);
			isValidAlpha(alpha4);
			isValidAlpha(alpha5);
		}
		else {
			printf("Error: invalid markov order at getTransProb()\n");
			exit(1);
		}
	}

	if(order == 0) {
		return transprob->zeroth[alpha0];
	}
	else if(order == 1) {
		return transprob->first[alpha1][alpha0];
	}
	else if(order == 2) {
		return transprob->second[alpha2][alpha1][alpha0];
	}
	else if(order == 3) {
		return transprob->third[alpha3][alpha2][alpha1][alpha0];
	}
	else if(order == 4) {
		return transprob->fourth[alpha4][alpha3][alpha2][alpha1][alpha0];
	}	
	else if(order == 5) {
		return transprob->fifth[alpha5][alpha4][alpha3][alpha2][alpha1][alpha0];
	}	
	else {
		printf("Error: invalid markov order at getTransProb()\n");
		exit(1);
	}
}



void copyTransProb(TransProb *src, TransProb *dest) {
	for(int a = 0; a < NUMALPHAS; a++) {
		dest->zeroth[a] = src->zeroth[a];
		for(int b = 0; b < NUMALPHAS; b++) {
			dest->first[a][b] = src->first[a][b];
			for(int c = 0; c < NUMALPHAS; c++) {
				dest->second[a][b][c] = src->second[a][b][c];
				for(int d = 0; d < NUMALPHAS; d++) {
					dest->third[a][b][c][d] = src->third[a][b][c][d];
					for(int e = 0; e < NUMALPHAS; e++) {
						dest->fourth[a][b][c][d][e] = src->fourth[a][b][c][d][e];
						for(int f = 0; f < NUMALPHAS; f++) {
							dest->fifth[a][b][c][d][e][f] = src->fifth[a][b][c][d][e][f];
						}
					}
				}
			}
		}
	}
}


//--------------------------------------------------------------------------------------------
// TransProb estimation
//
//--------------------------------------------------------------------------------------------
void resetTransCount(TransCount *transcount, const int bgPseudocount) {
	for(int a = 0; a < NUMCHARS; a++) {
		transcount->count0[a] = bgPseudocount;
		for(int b = 0; b < NUMCHARS; b++) {
			transcount->count1[a][b] = bgPseudocount;
			for(int c = 0; c < NUMCHARS; c++) {
				transcount->count2[a][b][c] = bgPseudocount;
				for(int d = 0; d < NUMCHARS; d++) {
					transcount->count3[a][b][c][d] = bgPseudocount;
					for(int e = 0; e < NUMCHARS; e++) {
						transcount->count4[a][b][c][d][e] = bgPseudocount;
						for(int f = 0; f < NUMCHARS; f++) {
							transcount->count5[a][b][c][d][e][f] = bgPseudocount;
						}
					}
				}
			}
		}
	}		
}

//*seq - array of sequence 
//len - sequence length
void tallyToTransCount(TransCount *transcount, int *seq, int len) {
	if(DEBUG0) {
		assert(transcount->maxorder <= MARKOV_ORDER_BOUND);
		assert(transcount->maxorder >= 0);
	}

	int maxorder = transcount->maxorder;

	//0th order
	if(maxorder >=0) {
		for(int j = 0; j < len; j++) {
			transcount->count0[ seq[j] ]++;
		}
	}

	//first-order
	if(maxorder >= 1) {
		//NOTE the difference between NUMCHARS and NUMALPHAS.
		//we count all the GAP_CHARS, but disregard them when computing the transprob
		for(int j = 1; j < len; j++) {
			transcount->count1[ seq[j-1] ][ seq[j] ]++;
		}
	}

	//second-order
	if(maxorder >= 2) {
		for(int j = 2; j < len; j++) {
			transcount->count2[ seq[j-2] ][ seq[j-1] ][ seq[j] ]++;
		}
	}

	if(maxorder >= 3) {
		for(int j = 3; j < len; j++) {
			transcount->count3[ seq[j-3] ][ seq[j-2] ][ seq[j-1] ][ seq[j] ]++;
		}
	}

	if(maxorder >= 4) {
		for(int j = 4; j < len; j++) {
			transcount->count4[ seq[j-4] ][ seq[j-3] ][ seq[j-2] ][ seq[j-1] ][ seq[j] ]++;
		}
	}

	if(maxorder >= 5) {
		for(int j = 5; j < len; j++) {
			transcount->count5[ seq[j-5] ][ seq[j-4] ][ seq[j-3] ][ seq[j-2] ][ seq[j-1] ][ seq[j] ]++;
		}
	}
}


TransProb* normalizeTransCount(TransCount *transcount) {
	TransProb *transprob = (TransProb*) malloc(sizeof(TransProb));
	transprob->maxorder = transcount->maxorder;
	int maxorder = transprob->maxorder;

	//0th order
	if(maxorder >=0) {
		int sum = 0;
		for(int a = 0; a < NUMALPHAS; a++) {
			sum += transcount->count0[a];
		}
		for(int a = 0; a < NUMALPHAS; a++) {
			transprob->zeroth[a] = ((double)transcount->count0[a]) / sum;
		}
	}

	//first-order
	if(maxorder >= 1) {
		//NOTE the difference between NUMCHARS and NUMALPHAS.
		//we count all the GAP_CHARS, but disregard them when computing the transprob

		for(int a = 0; a < NUMALPHAS; a++) {
			int sum = 0;
			for(int b = 0; b < NUMALPHAS; b++) {
				sum += transcount->count1[a][b];
			}
			for(int b = 0; b < NUMALPHAS; b++) {
				transprob->first[a][b] = ((double)transcount->count1[a][b]) / sum;
			}
		}
	}

	//second-order
	if(maxorder >= 2) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				int sum = 0;
				for(int c = 0; c < NUMALPHAS; c++) {
					sum += transcount->count2[a][b][c];
				}
				for(int c = 0; c < NUMALPHAS; c++) {
					transprob->second[a][b][c] = ((double)transcount->count2[a][b][c]) / sum;
				}
			}
		}

	}

	//third-order
	if(maxorder >= 3) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				for(int c = 0; c < NUMALPHAS; c++) {
					int sum = 0;
					for(int d = 0; d < NUMALPHAS; d++) {
						sum += transcount->count3[a][b][c][d];
					}
					for(int d = 0; d < NUMALPHAS; d++) {
						transprob->third[a][b][c][d] = ((double)transcount->count3[a][b][c][d]) / sum;
					}
				}
			}
		}
	}

	//fourth-order
	if(maxorder >= 4) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				for(int c = 0; c < NUMALPHAS; c++) {
					for(int d = 0; d < NUMALPHAS; d++) {
						int sum = 0;
						for(int e = 0; e < NUMALPHAS; e++) {
							sum += transcount->count4[a][b][c][d][e];
						}
						for(int e = 0; e < NUMALPHAS; e++) {
							transprob->fourth[a][b][c][d][e] = ((double)transcount->count4[a][b][c][d][e]) / sum;
						}
					}
				}
			}
		}
	}

	//fifth-order
	if(maxorder >= 5) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				for(int c = 0; c < NUMALPHAS; c++) {
					for(int d = 0; d < NUMALPHAS; d++) {
						for(int e = 0; e < NUMALPHAS; e++) {
							int sum = 0;
							for(int f = 0; f < NUMALPHAS; f++) {
								sum += transcount->count5[a][b][c][d][e][f];
							}
							for(int f = 0; f < NUMALPHAS; f++) {
								transprob->fifth[a][b][c][d][e][f] = ((double)transcount->count5[a][b][c][d][e][f]) / sum;
							}
						}
					}
				}
			}
		}
	}

	return transprob;
}


void marginalizeTransCount(TransCount *transcount) {
	int maxorder = transcount->maxorder;

	if(maxorder >= 5) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				for(int c = 0; c < NUMALPHAS; c++) {
					for(int d = 0; d < NUMALPHAS; d++) {
						for(int e = 0; e < NUMALPHAS; e++) {
							int sum = 0;
							for(int f = 0; f < NUMALPHAS; f++) {
								sum += transcount->count5[a][b][c][d][e][f];
							}
							transcount->count4[a][b][c][d][e] = sum;
						}
					}
				}
			}
		}
	}

	if(maxorder >= 4) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				for(int c = 0; c < NUMALPHAS; c++) {
					for(int d = 0; d < NUMALPHAS; d++) {
						int sum = 0;
						for(int e = 0; e < NUMALPHAS; e++) {
							sum += transcount->count4[a][b][c][d][e];
						}
						transcount->count3[a][b][c][d] = sum;
					}
				}
			}
		}
	}

	if(maxorder >= 3) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				for(int c = 0; c < NUMALPHAS; c++) {
					int sum = 0;
					for(int d = 0; d < NUMALPHAS; d++) {
						sum += transcount->count3[a][b][c][d];
					}
					transcount->count2[a][b][c] = sum;
				}
			}
		}
	}

	if(maxorder >= 2) {
		for(int a = 0; a < NUMALPHAS; a++) {
			for(int b = 0; b < NUMALPHAS; b++) {
				int sum = 0;
				for(int c = 0; c < NUMALPHAS; c++) {
					sum += transcount->count2[a][b][c];
				}
				transcount->count1[a][b] = sum;
			}
		}
	}

	if(maxorder >= 1) {
		for(int a = 0; a < NUMALPHAS; a++) {
			int sum = 0;
			for(int b = 0; b < NUMALPHAS; b++) {
				sum += transcount->count1[a][b];
			}
			transcount->count0[a] = sum;
		}
	}
}


//--------------------------------------------------------------------------------------------
// Print/Output functions
//
//--------------------------------------------------------------------------------------------

static
char numToChar(int num) {
	switch (num) {
		case 0: return 'A'; 
		case 1: return 'C'; 
		case 2: return 'G'; 
		case 3: return 'T'; 

		case GAP_CHAR: return 'X'; 
	}
	fprintf(stderr, "ERROR: numToChar has incorrect input\n");
	exit(1);

	return (char)NULL;
}

void printTransCount(FILE *fptr, TransCount *transcount, bool displayAll){
	for(int a = 0; a < NUMALPHAS; a++) {
		fprintf(fptr, "%c %7d ",  numToChar(a), transcount->count0[a]);

		if(displayAll) {
			for(int b = 0; b < NUMALPHAS && transcount->maxorder >= 1; b++) {
				fprintf(fptr, "%c%c %7d ", numToChar(b), numToChar(a), transcount->count1[b][a]);
				for(int c = 0; c < NUMALPHAS && transcount->maxorder >= 2; c++) {
					fprintf(fptr, "%c%c%c %7d ", numToChar(c), numToChar(b), numToChar(a), transcount->count2[c][b][a]);
					for(int d = 0; d < NUMALPHAS && transcount->maxorder >= 3; d++) {
						fprintf(fptr, "%c%c%c%c %7d ", numToChar(d), numToChar(c), numToChar(b), numToChar(a), transcount->count3[d][c][b][a]);
						for(int e = 0; e < NUMALPHAS && transcount->maxorder >= 4; e++) {
							fprintf(fptr, "%c%c%c%c%c %7d ", numToChar(e), numToChar(d), numToChar(c), numToChar(b), numToChar(a), transcount->count4[e][d][c][b][a]);
							for(int f = 0; f < NUMALPHAS && transcount->maxorder >= 5; f++) {
								fprintf(fptr, "%c%c%c%c%c%c %7d ", numToChar(f), numToChar(e), numToChar(d), numToChar(c), numToChar(b), numToChar(a), transcount->count5[f][e][d][c][b][a]);
							}
						}
					}
				}
			}
			fprintf(fptr, "\n");
		}
	}
	fprintf(fptr, "\n");
}


void printTransProb(FILE *fptr, TransProb *transprob, bool displayAll){
	for(int a = 0; a < NUMALPHAS; a++) {
		fprintf(fptr, "%c %7.4lf ",  numToChar(a), transprob->zeroth[a]);

		if(displayAll) {
			for(int b = 0; b < NUMALPHAS && transprob->maxorder >= 1; b++) {
				fprintf(fptr, "%c%c %7.4lf ", numToChar(b), numToChar(a), transprob->first[b][a]);
				for(int c = 0; c < NUMALPHAS && transprob->maxorder >= 2; c++) {
					fprintf(fptr, "%c%c%c %7.4lf ", numToChar(c), numToChar(b), numToChar(a), transprob->second[c][b][a]);
					for(int d = 0; d < NUMALPHAS && transprob->maxorder >= 3; d++) {
						fprintf(fptr, "%c%c%c%c %7.4lf ", numToChar(d), numToChar(c), numToChar(b), numToChar(a), transprob->third[d][c][b][a]);
						for(int e = 0; e < NUMALPHAS && transprob->maxorder >= 4; e++) {
							fprintf(fptr, "%c%c%c%c%c %7.4lf ", numToChar(e), numToChar(d), numToChar(c), numToChar(b), numToChar(a), transprob->fourth[e][d][c][b][a]);
							for(int f = 0; f < NUMALPHAS && transprob->maxorder >= 5; f++) {
								fprintf(fptr, "%c%c%c%c%c%c %7.4lf ", numToChar(f), numToChar(e), numToChar(d), numToChar(c), numToChar(b), numToChar(a), transprob->fifth[f][e][d][c][b][a]);
							}
						}
					}
				}
			}
			fprintf(fptr, "\n");
		}
	}
	fprintf(fptr, "\n");
}
