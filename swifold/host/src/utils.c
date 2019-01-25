#include "utils.h"




/*
void merge_scores(int * scores, char ** titles, unsigned long int size) {
	unsigned long int i1 = 0;
	unsigned long int i2 = size / 2;
	unsigned long int it = 0;
	// allocate memory for temporary buffers
	char ** tmp2 = (char **) malloc(size*sizeof(char *));
	int * tmp3 = (int *) malloc (size*sizeof(int));

	while(i1 < size/2 && i2 < size) {
		if (scores[i1] > scores[i2]) {
			tmp2[it] = titles[i1];
			tmp3[it] = scores[i1];
			i1++;
		}
		else {
			tmp2[it] = titles[i2];
			tmp3[it] = scores[i2];
			i2 ++;
		}
		it ++;
	}

	while (i1 < size/2) {
		tmp2[it] = titles[i1];
		tmp3[it] = scores[i1];
	    i1++;
	    it++;
	}
	while (i2 < size) {
		tmp2[it] = titles[i2];
		tmp3[it] = scores[i2];
	    i2++;
	    it++;
	}

	memcpy(titles, tmp2, size*sizeof(char *));
	memcpy(scores, tmp3, size*sizeof(int));

	free(tmp2);
	free(tmp3);

}


void mergesort_scores_serial(int * scores, char ** titles, unsigned long int size) {
	int tmp_score;
	char * tmp_seq;

	if (size == 2) { 
		if (scores[0] <= scores[1]) {
			// swap scores
			tmp_score = scores[0];
			scores[0] = scores[1];
			scores[1] = tmp_score;
			// swap titles
			tmp_seq = titles[0];
			titles[0] = titles[1];
			titles[1] = tmp_seq;
		}
	} else {
		if (size > 2){
			mergesort_scores_serial(scores, titles, size/2);
			mergesort_scores_serial(scores + size/2, titles + size/2, size - size/2);
			merge_scores(scores, titles, size);
		}
	}
}

void sort_scores (int * scores, char ** titles, unsigned long int size, int threads) {
    if ( threads == 1) {
	      mergesort_scores_serial(scores, titles, size);
    }
    else if (threads > 1) {
        #pragma omp parallel sections num_threads(threads)
        {
            #pragma omp section
            sort_scores(scores, titles, size/2, threads/2);
            #pragma omp section
            sort_scores(scores + size/2, titles  + size/2, size-size/2, threads-threads/2);
        }

        merge_scores(scores, titles, size);
    } // threads > 1
}
*/
// Wall time
double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}
