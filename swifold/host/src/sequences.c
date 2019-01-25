#include "sequences.h"

// Load sequence from file 
void load_sequence(char * sequence_filename, char ** ptr_sequence, char ** ptr_header, int * ptr_length, int * ptr_extended_length, int extend) {

	int i, j, extended_length, title_length, length, tmp_length;
	char * sequence=NULL, *title, buffer[BUFFER_SIZE], filename[BUFFER_SIZE], new_line='\n', *res;
	FILE * sequence_file;

	// open query sequence filename 
	sequence_file = fopen(sequence_filename,"r");

	if (sequence_file == NULL)	{
		printf("SWIFOLD: An error occurred while opening input sequence file.\n");
		exit(2);
	}

	// read title and residues	
	res = fgets(buffer,BUFFER_SIZE,sequence_file);
	title_length = 0;
	// read title
	while (strrchr(buffer,new_line) == NULL) {
		title_length += strlen(buffer);
		res = fgets(buffer,BUFFER_SIZE,sequence_file);
	}
	title_length += strlen(buffer) + 1;
	// read sequence
	length = 0;
	res = fgets(buffer,BUFFER_SIZE,sequence_file);
	while ((res != NULL)) {
		length += strlen(buffer)-1;
		res = fgets(buffer,BUFFER_SIZE,sequence_file);
	}

	if (extend) 
		extended_length = ceil((double) length / BLOCK_WIDTH) *BLOCK_WIDTH;
	else
		extended_length = length;

	// Rewind sequences database file 
	rewind(sequence_file);

	// Allocate memory for title buffer
	posix_memalign((void**)&title, AOCL_ALIGNMENT, title_length*sizeof(char));
	if (title == NULL) { printf("SWIFOLD: An error occurred while allocating memory.\n"); exit(1); }

	res = fgets(buffer,BUFFER_SIZE,sequence_file);
	// read header
	tmp_length = 1;
	do{
		strncpy(title+(tmp_length-1),buffer,strlen(buffer)-1);
		tmp_length += strlen(buffer)-1;
		res = fgets(buffer,BUFFER_SIZE,sequence_file);
	} while (strrchr(buffer,new_line) == NULL);
	title[tmp_length] = '\0';

	// Rewind sequences database file 
	rewind(sequence_file);

	// Allocate memory for sequence
	posix_memalign((void**)&sequence, AOCL_ALIGNMENT, extended_length*sizeof(char));
	if (sequence == NULL) { printf("SWIFOLD: An error occurred while allocating memory.\n"); exit(1); }

	// Read sequence from the file and load them in sequence buffer 
	res = fgets(buffer,BUFFER_SIZE,sequence_file);
	// read title
	while (strrchr(buffer,new_line) == NULL)
		res = fgets(buffer,BUFFER_SIZE,sequence_file);
	// read sequence
	tmp_length = 1;
	res = fgets(buffer,BUFFER_SIZE,sequence_file);
	while ((res != NULL)) {
		//printf("%s %d\n",buffer,strlen(buffer));
		strncpy(sequence+(tmp_length-1),buffer,strlen(buffer)-1);
		tmp_length += strlen(buffer)-1;
		res = fgets(buffer,BUFFER_SIZE,sequence_file);
	}

	// Close sequences database file 
	fclose(sequence_file);

	for (i=length; i < extended_length ; i++)
		sequence[i] = DUMMY_ELEMENT;

	*ptr_sequence = sequence;
	*ptr_length = length;
	*ptr_extended_length = extended_length;
	*ptr_header = title;

}

