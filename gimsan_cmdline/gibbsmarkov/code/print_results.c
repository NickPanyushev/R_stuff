#include "print_results.h"
#include "gibbs_util.h"
//#include "heap.h"
//#include "multi_motifs_perm.h"

#define FLANKING_LEN 5

/* Prints out a given motif profile along with its occurences in the sample  */

static
void printGapPos(FILE *fptr, Profile *profile ) {
	for(int i = 0; i < profile->span; i++) {
		if(profile->isgap[i]) {
			fprintf(fptr, "Gap position: %d\n", i);
		}
	}
}

static
void printAlignment(FILE *fptr, int *sites, int span, double *scorePerSite, Dataset *data ) {
	int extra = 5; //extra to view on left and right
	char upToLow = 'a' - 'A'; //offset to convert uppercase to lowercase
	for(int i = 0; i < data->numseqs; i++) {
		if(sites[i] >= 0) {
			fprintf(fptr, "%3d   ", i);
			for(int j = -extra; j < span + extra; j++) {
				//print site index
				if(j == 0) {
					fprintf(fptr, " %6d ", abs(getConcat2DoubleStrandPos(data, i, sites[i])) );
				}
				if( j == span) { 
					fprintf(fptr, " %-6d ", abs(getConcat2DoubleStrandPos(data, i, sites[i] + span - 1)) );
				}


				//print alignment chars
				int pos = j + sites[i];
				if(pos < 0 || pos >= data->seqlen[i] 
						|| (isForwardStrand(data, i, sites[i]) && !isForwardStrand(data,i, pos)) 
						|| (!isForwardStrand(data, i, sites[i]) && isForwardStrand(data,i, pos)) 
				  )
				{
					fprintf(fptr, " "); //out of sequence range
				}
				else if(sites[i] <= pos && pos < sites[i] + span) {
					fprintf(fptr, "%c", numToChar(data->seqs[i][pos])); //upper case
				}
				else {
					//lower case
					fprintf(fptr, "%c", numToChar(data->seqs[i][pos]) + upToLow); 
				}
			}
			fprintf(fptr, "  (%c)", isForwardStrand(data, i, sites[i]) ? '+' : '-' );

			fprintf(fptr, " Score=%8.3lf", scorePerSite[i]);

			//not safe - does not prevent buffer overflow
			//but the header has capacity from dataset.h
			fprintf(fptr, "  %s", data->headers[i]);

			fprintf(fptr, "\n");
		}
	}
}

static
void printZeroOccurrencesSeqHeaders(FILE *fptr, int *sites, Dataset *data ) {
	for(int i = 0; i < data->numseqs; i++) {
		if(sites[i] < 0) {
			fprintf(fptr, "%3d ", i);
			//not safe - does not prevent buffer overflow
			fprintf(fptr, " %s", data->headers[i]);

			fprintf(fptr, "\n");
		}
	}
}


static
void printFlankingRegionProfile(FILE *fptr, int *sites, int span, Dataset *data) {
	double left[FLANKING_LEN][NUMCHARS];
	double right[FLANKING_LEN][NUMCHARS];

	for(int m = 0; m < FLANKING_LEN; m++) {
		for(int a = 0; a < NUMCHARS; a++) {
			left[m][a] = 0.0;
			right[m][a] = 0.0;
		}
	}

	for(int i = 0; i < data->numseqs; i++) {
		if(sites[i] >= 0) { //non-empty site under ZOOPS
			for(int m = 0; m < FLANKING_LEN; m++) {
				int lpos = sites[i] - m - 1;
				int rpos = sites[i] + span + m;

				if(lpos < 0 //out-of-bound on the left
					|| (isForwardStrand(data, i, sites[i]) && !isForwardStrand(data,i, lpos)) 
					|| (!isForwardStrand(data, i, sites[i]) && isForwardStrand(data,i, lpos))
				  )
				{
					left[m][GAP_CHAR] += 1.0; //out of range are counted as 'X'
				}
				else {
					left[m][data->seqs[i][lpos]] += 1.0;
				}

				if(rpos >= data->seqlen[i] //not out-of-bound on the right
					|| (isForwardStrand(data, i, sites[i]) && !isForwardStrand(data,i, rpos)) 
					|| (!isForwardStrand(data, i, sites[i]) && isForwardStrand(data,i, rpos))
				  )
				{
					right[m][GAP_CHAR] += 1.0; //out of range are counted as 'X'
				}
				else {
					right[m][data->seqs[i][rpos]] += 1.0;
				}
			}
		}
	}

	//normalize (including gap-chars)
	for(int m = 0; m < FLANKING_LEN; m++) {
		double lsum = 0.0; 
		double rsum = 0.0; 
		for(int a = 0; a < NUMCHARS; a++) {
			lsum += left[m][a];
			rsum += right[m][a];
		}
		for(int a = 0; a < NUMCHARS; a++) {
			left[m][a] /= lsum;
			right[m][a] /= rsum;
		}
	}

	fprintf(fptr, "Pos ");
	for(int x = -FLANKING_LEN; x <= -1; x++) {
		fprintf(fptr, " [%2d] ", x);
	}
	fprintf(fptr, "  ||  ");
	for(int x=0; x< FLANKING_LEN; x++) {
		fprintf(fptr, " [%2d] ", x + span);
	}
	fprintf(fptr, "\n");

	for(int y=0; y< NUMCHARS;  y++)
	{
		fprintf(fptr, "%c   ", numToChar(y));

		//print left flank
		for(int x = FLANKING_LEN - 1; x >= 0; x--) {
			fprintf(fptr, "%5.2lf ", left[x][y]);
		}

		fprintf(fptr, "  ||  ");

		//print right flank
		for(int x=0; x< FLANKING_LEN; x++) {
			fprintf(fptr, "%5.2lf ", right[x][y]);
		}
		fprintf(fptr, "\n");
	}
	fprintf(fptr, "\n");
}

static
void printCountmatAndSites(Profile *countmat, int *sites, double *scorePerSite, Dataset *data)
{
	updateCountmatFromSites(countmat, sites, data);
	Profile *profile = initProfile(countmat->span, countmat->maxspan);
	copyProfile(profile, countmat);

	//normalize "profile"
	for(int x = 0; x < countmat->span; x++) {
		if( !profile->isgap[x]) {
			double sum = 0.0;
			for(int y = 0; y < NUMALPHAS; y++) {
				sum += profile->mat[x][y];
			}
			for(int y = 0; y < NUMALPHAS; y++) {
				profile->mat[x][y] /= sum;
			}
		}
	}

	/* print profile */
	printf("Motif profile generated from sites:\n"); 
	printProfile(stdout, profile);

	printf("Flanking regions:\n"); 
	printFlankingRegionProfile(stdout, sites, countmat->span, data);


	/* NOW, PRINT CONSENSUS MOTIF */
	printf("Motif span: %d\n", profile->span);
	printf("Motif weight/width: %d\n", profile->cols);
	printf("Consensus motif: ");
	for(int x=0; x< profile->span; x++) {
		if(!profile->isgap[x]) {
			if(profile->mat[x][0] > 0.5) 
				printf("A");
			else if(profile->mat[x][1] > 0.5) 
				printf("C");
			else if(profile->mat[x][2] > 0.5) 
				printf("G");
			else if(profile->mat[x][3] > 0.5) 
				printf("T");
			else 
				printf("N");
		}
		else {
			printf("_");
		}
	}
	printf("\n");

	int numsites = getNumSites(sites, data->numseqs);
	printf("Number of predicted sites: %d (%.2lf%%)\n", numsites, (numsites * 100.0) / data->numseqs);
	printf("\n");

	if(profile->span != profile->cols) {
		printGapPos(stdout, profile);
		printf("\n");
	}

	printf("Sequence range is from 0 to (number of sequences - 1).\n");
	printf("Motif sites are in range 1 to length or -1 to -length.\n");
	for(int i=0; i< data->numseqs; i++) {
		if(sites[i] >= 0) {
			int pos = getConcat2DoubleStrandPos(data,i,sites[i]);
			int revcomplPos = getConcat2DoubleStrandPos(data, i, getConcatPosOfOppStrand(data, i, sites[i], profile->span));
			printf("Sequence %3d: motif site %6d  (%6d)\n",i, pos, revcomplPos);
		}
		else {
			printf("Sequence %3d: motif site %6s  (%6s)\n",i, "none", "none");
		}
	}
	

	printf("\n");
	printf("Alignments and FASTA-headers:\n");
	printAlignment(stdout, sites, profile->span, scorePerSite, data);
	printf("\n");
	if(numsites < data->numseqs) {
		printf("FASTA-header of sequences with zero occurrences:\n");
		printZeroOccurrencesSeqHeaders(stdout, sites, data);
		printf("\n");
	}
	nilProfile(profile);
}

void printRunNode(RunNode *node, enum ScoreMetric metric, Markov *markov, Dataset *data, Zoops *zoops) {
	if(DEBUG0) {
		assert(node->countmat->span == node->pswm->span);
	}

	printf("All of the following scores are computed after EM (if applicable).\n");
	if(metric == ILR) {
		//ILR from PWM
		double ilr_pwm = computeIlr(markov, data, node->pswm->mat, node->pswm->span);
		printf("Log Markovian-ILR of PWM: %.2lf\n", ilr_pwm);


		if(DEBUG0) {
			assert(fabs(ilr_pwm - node->score) < 0.00000001);
		}
	}

	//ILR from sites
	double ilr_sites = computeIlr(markov, data, zoops, node->sites, node->countmat->span);
	printf("Log Markovian-ILR generated from sites without pseudocount: %.2lf\n", ilr_sites);

	//CLR from sites
	double entscore= computeEntropy(markov, data, zoops, node->sites, node->countmat->span);
	printf("Log Markovian-CLR generated from sites without pseudocount: %.2lf\n", entscore);

	printf("Log(CLR) normalized by number of sequences: %.4lf\n", entscore / data->numseqs);

	//Score that was used as the best run
	printf("Score for ranking runs: %.2lf\n", node->score);
	printf("\n");

	//printf("PWM after EM:\n");
	//printProfile(stdout, node->pswm); //hopefully, no one screwed around with this matrix
	//printf("\n");

	//compute scores for each site
	double *scorePerSite = (double*) malloc(sizeof(double) * data->numseqs);
	for(int i = 0; i < data->numseqs; i++) {
		if(node->sites[i] >= 0) {
			scorePerSite[i] = computeEntropyPerSite(markov, data, zoops, node->sites, i, node->countmat->span);
		}
		else {
			scorePerSite[i] = NAN;
		}
	}

	printCountmatAndSites(node->countmat, node->sites, scorePerSite, data);

	free(scorePerSite);
}

////heap_rank is the rank if all motifs were sorted disregarding overlap/non-overlap
//static
//void recur_printTopRankRunNodes(int rank, int heap_rank, int remain, Heap *heap, MultiMotifs *multiMotifs,
//								enum ScoreMetric metric, Markov *markov, Dataset *data, Zoops *zoops) 
//{
//	const double threshold = multiMotifs->distCutoff; //increase to allow more similarity for motifs
//
//	//stopping condition
//	if(remain == 0 || heap->size == 0) {
//		return;
//	}
//
//	if(DEBUG0) {
//		if(multiMotifs->size != rank - 1) {
//			fprintf(stderr, "Inconsistent multiMotifs size with ranking");
//			exit(1);
//		}
//	}
//
//	//remove next best
//	RunNode *node = deleteRootNode(heap);
//	
//	//reject motif for print if > threshold percentage of sites are marked
//	bool rejectMotif = false;
//
//	//Decide to accept/reject motif
//	double current_min_dist = DBL_MAX;
//	for(MotifNode *m = multiMotifs->head; m!= NULL; m= m->next) {
//		double dist = computeMotifDist(node->countmat, m->countmat, multiMotifs);
//		if(dist < current_min_dist) {
//			current_min_dist = dist;
//		}
//		if(dist < threshold) {
//			rejectMotif = true;
//			break;
//		}
//	}
//
//	if(!rejectMotif) {
//		//print ranking
//		printf("-------------------------------------------------------------------\n");
//		printf("  RANK %d   (rank %d out of all %d motifs)\n", rank, heap_rank, heap->capacity);
//		if(rank > 1) {
//			printf("  Min distance over motifs with better ranks: %.5lg (cutoff: %.5lg)\n", current_min_dist, threshold);
//		}
//		else {
//			printf("  Min distance over motifs with better ranks: NaN (cutoff: %.5lg)\n", threshold);
//		}
//		printf("-------------------------------------------------------------------\n");
//		printRunNode(node, metric, markov, data, zoops);
//
//		//Add motif to the list of printed motif
//		prependToMultiMotifsList(node->countmat, multiMotifs);
//	}
//
//	if(DEBUG0) {
//		assert(heap_rank == heap->capacity - heap->size);
//	}
//
//	//recursion
//	if(!rejectMotif) {
//		recur_printTopRankRunNodes(rank+1, heap_rank+1, remain-1, heap, multiMotifs, metric, markov, data, zoops);
//	}
//	else {
//		recur_printTopRankRunNodes(rank, heap_rank+1, remain, heap, multiMotifs, metric, markov, data, zoops);
//	}
//
//}

//void printTopRankRunNodes(int numTop, RunSet *runset, enum ScoreMetric metric, 
//						  Markov *markov, Dataset *data, Zoops *zoops, double pvalCutoff)
//{
//	//const double pvalCutoff = 0.30;
//	const int numPermPairs = 1000; //number of permutation pairs
//
//	MultiMotifs *multiMotifs = initMultiMotifs(numTop, data->numseqs, data->maxspan);
//	setDistCutoffForMultiMotifs(multiMotifs, runset, numPermPairs, pvalCutoff);
//
//	//sort and print
//	Heap *heap = constructHeap(runset->len);
//	for(RunNode *node = runset->head; node!= NULL; node = node->next) {
//		processNode(heap, node);
//	}
//	recur_printTopRankRunNodes(1, 1, numTop, heap, multiMotifs, metric, markov, data, zoops);
//	destructHeap(heap);
//
//	nilMultiMotifs(multiMotifs);
//}

