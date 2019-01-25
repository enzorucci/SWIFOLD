/*
############################################
SW kernel (int_bw512 version). Single inner-loop computation.
############################################
*/

#define BLOCK_WIDTH 512 // update this value in host/src/arguments.h in case of modification

// AOC kernel demonstrating device-side printf call
__attribute__((reqd_work_group_size(1,1,1)))
__attribute__((num_compute_units(1)))
__attribute__((task))
__kernel void sw (	__global const char * restrict a, const int m, 
					__global const char * restrict b, const int nbb,
					const int open_extend_gap, const int extend_gap, const int match, const int mismatch,
					__global const int * restrict prev_lastCol,  __global int * restrict next_lastCol,
					__global const int * restrict prev_maxRow, __global int * restrict next_maxRow,
					__global int * restrict maxScore, const int jj){

	// auxiliar buffers
	int row1[BLOCK_WIDTH]={0}, maxCol1[BLOCK_WIDTH]={0};
	int score=0, auxLastCol=0;

	// private buffer for sequence b
	char private_b1[BLOCK_WIDTH];

	// pointer to sequence b
	__global const char * ptr_b = b + jj*BLOCK_WIDTH;

	// copy sequence b to private memory
	for(int i = 0; i < BLOCK_WIDTH; i++)
		private_b1[i] = ptr_b[i];

	for(int i = 0; i < m; i++){

		// copy resiude to local memory
		char a_i = a[i];
		int previous = 0;
		int maxRow_i = 0;
		int score_i = 0;
		if (jj != 0){
			// update row[0] with lastCols[i]
			row1[0] = prev_lastCol[i];
			maxRow_i = prev_maxRow[i];
		}

		#pragma unroll
		for (int j=0; j < BLOCK_WIDTH ; j++){
			//calcuate the diagonal value
			int current = row1[j] + (a_i==private_b1[j] ? match : -mismatch);
			// calculate current max value
			current = max(current, maxRow_i);
			current = max(current, maxCol1[j]);
			current = max(current, 0);
			// update max score
			score_i = max(score_i, current);
			// update maxRow and maxCol
			int aux1 = maxRow_i - extend_gap;
			int aux2 = maxCol1[j] - extend_gap;
			int aux3 = current - open_extend_gap;
			maxRow_i = max(aux1, aux3);
			maxCol1[j] =  max(aux2, aux3);	
			// update row buffer
			row1[j] = previous;
			previous = current;
		}

		if (jj != nbb-1) {
			// update lastCol
			next_lastCol[i] = auxLastCol;
			auxLastCol = previous;
			// update maxRow
			next_maxRow[i] = maxRow_i;
		}

		// update max score
		score = max(score, score_i);
	}

	// save max score of current block in global memory
	maxScore[0] = max(maxScore[0], score);

}
