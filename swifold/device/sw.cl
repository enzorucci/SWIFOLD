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

/*
############################################
SW kernel (int_bw1152 version). Multiple inner-loop computation.
############################################

As mentioned in the paper,  larger BW values produce better performance results but at the expense of higher resource consumption. 
So, you have to seek the largest BW value supported on your device in order to obtain the maximum performance.
It is important to remark that the AOC reports the appearance of non-real read-write dependences in private memory associated to matrices H (row) and F (maxCol) 
after a certain BW value, which aborts binary kernel generation. In order to solve this issue, the innermost loop is split into two or more loops to 
carry out the execution of wider blocks. The kernel above splits the inner loop into 3 parts to get a wider block.



#define BLOCK_WIDTH 1152  // update this value in host/src/arguments.h in case of modification

#define BLOCK_WIDTH1 512 
#define BLOCK_WIDTH2 512
#define BLOCK_WIDTH3 128

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
	int row1[BLOCK_WIDTH1]={0}, maxCol1[BLOCK_WIDTH1]={0}, row2[BLOCK_WIDTH2]={0}, maxCol2[BLOCK_WIDTH2]={0}, row3[BLOCK_WIDTH3]={0}, maxCol3[BLOCK_WIDTH3]={0};
	int score=0, auxLastCol=0;

	// private buffer for sequence b
	char private_b1[BLOCK_WIDTH1], private_b2[BLOCK_WIDTH2], private_b3[BLOCK_WIDTH3];

	// pointer to sequence b
	__global const char * ptr_b = b + jj*BLOCK_WIDTH;

	// copy sequence b to private memory
	for(int i = 0; i < BLOCK_WIDTH1; i++)
		private_b1[i] = ptr_b[i];

	ptr_b += BLOCK_WIDTH1;

	// copy sequence b to private memory
	for(int i = 0; i < BLOCK_WIDTH2; i++)
		private_b2[i] = ptr_b[i];

	ptr_b += BLOCK_WIDTH2;

	// copy sequence b to private memory
	for(int i = 0; i < BLOCK_WIDTH3; i++)
		private_b3[i] = ptr_b[i];

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
		for (int j=0; j < BLOCK_WIDTH1 ; j++){
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

		#pragma unroll
		for (int j=0; j < BLOCK_WIDTH2 ; j++){
			//calcuate the diagonal value
			int current = row2[j] + (a_i==private_b2[j] ? match : -mismatch);
			// calculate current max value
			current = max(current, maxRow_i);
			current = max(current, maxCol2[j]);
			current = max(current, 0);
			// update max score
			score_i = max(score_i, current);
			// update maxRow and maxCol
			int aux1 = maxRow_i - extend_gap;
			int aux2 = maxCol2[j] - extend_gap;
			int aux3 = current - open_extend_gap;
			maxRow_i = max(aux1, aux3);
			maxCol2[j] =  max(aux2, aux3);	
			// update row buffer
			row2[j] = previous;
			previous = current;
		}

		#pragma unroll
		for (int j=0; j < BLOCK_WIDTH3 ; j++){
			//calcuate the diagonal value
			int current = row3[j] + (a_i==private_b3[j] ? match : -mismatch);
			// calculate current max value
			current = max(current, maxRow_i);
			current = max(current, maxCol3[j]);
			current = max(current, 0);
			// update max score
			score_i = max(score_i, current);
			// update maxRow and maxCol
			int aux1 = maxRow_i - extend_gap;
			int aux2 = maxCol3[j] - extend_gap;
			int aux3 = current - open_extend_gap;
			maxRow_i = max(aux1, aux3);
			maxCol3[j] =  max(aux2, aux3);	
			// update row buffer
			row3[j] = previous;
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
*/