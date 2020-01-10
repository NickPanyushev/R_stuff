#include "siteloc.h"

Siteloc* constructDegenerateSiteloc(int numsites, int default_dspos) {
	Siteloc *siteloc = (Siteloc*) malloc(sizeof(Siteloc));
	siteloc->numsites = numsites;

	siteloc->seqInd = (int*) malloc(sizeof(int) * numsites);
	siteloc->dsPos = (int*) malloc(sizeof(int) * numsites);

	if(DEBUG0) {
		if(default_dspos == 0) {
			//range is [1, m] union [-m, -1]
			abort();
		}
	}

	for(int i = 0; i < numsites; i++) {
		//assume one site per sequence and they all start at default_pos
		siteloc->seqInd[i] = i;
		siteloc->dsPos[i] = default_dspos;
	}
	return siteloc;
}

static 
int determineNumSites(FILE *fptr) {
	int numsites = 0;

	rewind(fptr);
	char c;
	//if the line contains a digit
	bool hasDigit = false;
	while((c = fgetc(fptr)) != EOF) {
		if(c >= '0' && c <= '9') {
			hasDigit = true;
		}
		else if( c == '\n' || c== '\r') {
			if(hasDigit) {
				numsites++;
			}
			hasDigit = false;
		}
		else if( ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z')) {
			//if contain alphabet
			fprintf(stderr, "Error: site location file contains alphabet.\n");
			abort();
		}
	}
	return numsites;

}

Siteloc* openSiteloc(char *filename) {
	FILE *fptr;
	if((fptr = fopen(filename,"r")) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n",filename);
		fprintf(stderr,"File does not exist!\n");
		exit(1);
	}
	fclose(fptr);

	int numsites = determineNumSites(fptr);

	//construct the siteloc
	Siteloc *siteloc = (Siteloc*) malloc(sizeof(Siteloc));
	siteloc->numsites = numsites;
	siteloc->seqInd = (int*) malloc(sizeof(int) * numsites);
	siteloc->dsPos = (int*) malloc(sizeof(int) * numsites);

	//read the sites
	rewind(fptr);
	int siteCount = 0;
	int intCount = 0;
	while(1) {
		int num = 0;
		int ret = fscanf(fptr, "%d", &num);
		if(ret == EOF || ret == 0) {
			break;
		}
		else {
			if(siteCount >= siteloc->numsites) {
				fprintf(stderr, "Error: site location file is not in valid format.\n");
				fprintf(stderr, "Number of sites found and newline count are inconsistent.\n");
				exit(1);
			}

			if(intCount % 2 == 0) {
				//seq index
				siteloc->seqInd[siteCount] = num;
			}
			else {
				//dspos index
				siteloc->dsPos[siteCount] = num;
				siteCount++;
			}
			intCount++;
		}
	}
	return siteloc;
}

void destructSiteloc(Siteloc *siteloc) {
	free(siteloc->dsPos);
	free(siteloc->seqInd);
	free(siteloc);
}

