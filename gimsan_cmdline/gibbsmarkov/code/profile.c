#include "profile.h"
#include "dataset.h"

Profile *initProfile(int initspan, int maxspan) {
	Profile *p;
	if(initspan > maxspan) {
		fprintf(stderr, "Error: initspan > maxspan in initProfile()\n");
	}
	p = (Profile*) malloc(sizeof(Profile));
	p->cols = initspan;
	p->span = initspan;
	p->maxspan = maxspan;
	p->maxspan = p->maxspan;

	p->isgap = (bool*) calloc(p->maxspan, sizeof(bool)); 
	p->mat = (double**) malloc(sizeof(double*) * p->maxspan);
	for(int i = 0; i < p->maxspan; i++) {
		p->mat[i] = (double*) malloc(sizeof(double) * NUMALPHAS);
	}
	for(int i = 0; i < p->span; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			p->mat[i][j] = 0.0;
		}
	}
	return p;
}

void nilProfile(Profile *profile) {
	int i;
	free(profile->isgap);
	for( i = 0; i < profile->maxspan; i++) {
		free(profile->mat[i]);
	}
	free(profile->mat);
	free(profile);
}

void copyProfile(Profile *dest, Profile *src) {
	if(DEBUG0) {
		if(dest->maxspan != src->maxspan) {
			fprintf(stderr, "Error: dest->maxspan != src->maxspan in copyProfile()\n");
			exit(1);
		}
		for(int i = 0; i < src->span; i++) {
			for(int j = 0; j < NUMALPHAS; j++) {
				if(src->isgap[i] && fabs(src->mat[i][j] - 1.0) > 0.000001){
					fprintf(stderr, "Error: values on gaps are not set to 1.0 in copyProfile()\n");
					exit(1);
				}
			}
		}
	}
	dest->cols = src->cols;
	dest->span = src->span;
	for(int i = 0; i < src->span; i++) {
		dest->isgap[i] = src->isgap[i];
	}
	for(int i = 0; i < src->span; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			dest->mat[i][j] = src->mat[i][j];
		}
	}
}

void resetProfile(Profile *profile, int initspan) {
	if(DEBUG0) {
		if(initspan > profile->maxspan) {
			fprintf(stderr, "Error: initspan > maxspan in resetProfile()\n");
			exit(1);
		}
	}
	profile->cols = initspan;
	profile->span = initspan;
	for(int i = 0; i < profile->span; i++) {
		profile->isgap[i] = 0;
	}
	for(int i = 0; i < profile->span; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			profile->mat[i][j] = 0.0;
		}
	}
}
void setProfile(Profile *profile, double val) {
	for(int i = 0; i < profile->span; i++) {
		for(int j = 0; j < NUMALPHAS; j++) {
			if(!profile->isgap[i]) {
				profile->mat[i][j] = val;
			}
		}
	}
	if(DEBUG0) {
		if(profile->isgap[0]) {
			fprintf(stderr, "Error: profile column 0 is a gap in setProfile()\n");
			exit(1);
		}
		if(profile->span > profile->maxspan) {
			fprintf(stderr, "Error: span out of bound in setProfile().\n");
			exit(1);
		}
		for(int i = 0; i < profile->span; i++) {
			for(int j = 0; j < NUMALPHAS; j++) {
				if(profile->isgap[i] && fabs(profile->mat[i][j] - 1.0) > 0.000001){
					fprintf(stderr, "Error: values on gaps are not set to 1.0 in setProfile()\n");
					exit(1);
				}
			}
		}
	}
}

//In-place reverse-complement of a PWM
void revcomplProfile(Profile *profile) {
	//reverse isgap and profile
	int left = 0; 
	int right = profile->span - 1;
	while(left < right) {
		//reverse profile
		double *tempcol = profile->mat[left];
		profile->mat[left] = profile->mat[right];
		profile->mat[right] = tempcol;

		//reverse isgap
		bool temp = profile->isgap[left];
		profile->isgap[left] = profile->isgap[right];
		profile->isgap[right] = temp;

		left++;
		right--;
	}

	//complement
	for(int i = 0; i < profile->span; i++) {
		double temp;

		//swap A/T (0 and 3)
		temp = profile->mat[i][0];
		profile->mat[i][0] = profile->mat[i][3];
		profile->mat[i][3] = temp;

		//swap C/G (1 and 2)
		temp = profile->mat[i][1];
		profile->mat[i][1] = profile->mat[i][2];
		profile->mat[i][2] = temp;
	}
}

//----------------------------------------------------------------------
// Column Shift Profile
//----------------------------------------------------------------------
//Shift profile on a set of sequences to the left
void shiftProfileLeft(Profile *profile, double *initval) {
	rmProfileCol(profile, profile->span-1);
	addFrontCol(profile, initval);
}
void shiftProfileRight(Profile *profile, double *initval) {
	rmProfileCol(profile, 0);
	addBackCol(profile, initval);
}

//----------------------------------------------------------------------
// Dimension change
//----------------------------------------------------------------------
void addFrontCol(Profile *profile, double *initval) {
	if(DEBUG0) {
		if(profile->maxspan == profile->span) { //cannot expand span
			fprintf(stderr, "Error: maxspan == span in addFrontCol");
			exit(1);
		}
	}
	for(int  i = profile->span; i >= 1; i--) {
		for(int j = 0; j < NUMALPHAS; j++) {
			profile->mat[i][j] = profile->mat[i-1][j];
		}
		profile->isgap[i] = profile->isgap[i-1];	
	}

	//set values for the new column
	for(int i = 0; i < NUMALPHAS; i++) {
		profile->mat[0][i] = initval[i];
	}
	profile->isgap[0] = 0;

	profile->span++;
	profile->cols++;

}
void addBackCol(Profile *profile, double *initval) {
	if(DEBUG0) {
		if(profile->maxspan == profile->span) { //cannot expand span
			fprintf(stderr, "Error: maxspan == span in addFrontCol");
			exit(1);
		}
	}

	profile->isgap[profile->span] = 0;
	for(int i = 0; i < NUMALPHAS; i++) {
		profile->mat[profile->span][i] = initval[i];
	}

	profile->span++;
	profile->cols++;

}
void addProfileCol(Profile *profile, int pos, double *initval) {
	if(DEBUG0) {
		if(profile->cols == profile->span) { //no gaps to fill
			fprintf(stderr, "Error: no gaps to fill in addProfileCol");
			exit(1);
		}
		if(pos <= 0 || pos >= profile->span - 1) { //addition case of boundaries
			fprintf(stderr, "Error: pos out of range in addProfileCol() %d\n", pos);
			exit(1);
		}
		if(!profile->isgap[pos]) {
			fprintf(stderr, "Error: pos is not a gap in addProfileCol()\n");
			exit(1);
		}
	}
	profile->isgap[pos] = 0;
	profile->cols++;
	for(int i = 0; i < NUMALPHAS; i++) {
		profile->mat[pos][i] = initval[i];
	}
}
void rmProfileCol(Profile *profile, int pos) {
	if(DEBUG0) {
		if(profile->cols <= 1) {
			fprintf(stderr, "Error: profile->cols <= 1 in rmProfilecol");
			exit(1);
		}
		if(pos < 0 || pos >= profile->span) {
			fprintf(stderr, "Error: pos out of range in rmProfileCol() %d\n", pos);
			exit(1);
		}
		if(profile->isgap[pos]) {
			fprintf(stderr, "Error: pos is a gap in rmProfileCol()\n");
			exit(1);
		}
	}
	if(pos == 0) {
		int startpos;
		for(startpos = 1; startpos < profile->span; startpos++) {
			if(!profile->isgap[startpos]) {
				break;
			}
		}
		for(int i = startpos; i < profile->span; i++) {
			for(int a = 0; a < NUMALPHAS; a++) {
				profile->mat[i - startpos][a] = profile->mat[i][a];
			}
			profile->isgap[i - startpos] = profile->isgap[i];

		}
		profile->span -= startpos;
	}
	else if(pos == profile->span - 1) {
		profile->span--;
		while(profile->isgap[profile->span -1]) {
			profile->span--;
		}
	}
	else { //pos is in the middle
		for(int a = 0; a < NUMALPHAS; a++) {
			profile->mat[pos][a] = 1.0;
		}
		profile->isgap[pos] = 1; 
	}
	profile->cols--;

	if(DEBUG0) {
		if(profile->cols > profile->span) {
			fprintf(stderr, "Error: profile->cols > profile->span in rmProfileCol()\n");
			fprintf(stderr, "span %d, cols %d\n", profile->span, profile->cols);
			exit(1);
		}
	}
}
//----------------------------------------------------------------------
// Print
//----------------------------------------------------------------------
void printProfile(FILE *fptr, Profile *profile) {
	fprintf(fptr, "Pos ");
	for(int x=0; x< profile->span; x++) {
		fprintf(fptr, " [%2d]", x);
		fprintf(fptr, " ");
	}
	fprintf(fptr, "\n");

	for(int y=0; y< NUMALPHAS;  y++)
	{
		fprintf(fptr, "%c   ", numToChar(y));
		for(int x=0; x< profile->span; x++) {
			if(profile->isgap[x]) {
				fprintf(fptr, "  ** ");
				fprintf(fptr, " ");
			}
			else {
				fprintf(fptr, "%5.2lf",profile->mat[x][y]);
				fprintf(fptr, " ");
			}
		}
		fprintf(fptr, "\n");
	}
	fprintf(fptr, "\n");
}

void printCountmat(FILE *fptr, Profile *countmat) {
	fprintf(fptr, "Pos ");
	for(int x=0; x< countmat->span; x++) {
		fprintf(fptr, " [%2d]", x);
		fprintf(fptr, " ");
	}
	fprintf(fptr, "\n\n");

	for(int y=0; y< NUMALPHAS;  y++)
	{
		fprintf(fptr, "%c   ", numToChar(y));
		for(int x=0; x< countmat->span; x++) {
			if(countmat->isgap[x]) {
				fprintf(fptr, "  ** ");
				fprintf(fptr, " ");
			}
			else {
				//no decimal points
				fprintf(fptr, "%5.0lf",countmat->mat[x][y]);
				fprintf(fptr, " ");
			}
		}
		fprintf(fptr, "\n\n");
	}
}

