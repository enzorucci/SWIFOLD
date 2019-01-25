#ifndef DB_H_INCLUDED
#define DB_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>
#include "arguments.h"
#include "utils.h"

#define BUFFER_SIZE 1000
#define ALLOCATION_CHUNK 1000
#define AOCL_ALIGNMENT 64
#define DUMMY_ELEMENT 'Z'

void load_sequence(char * sequence_filename, char ** ptr_sequence, char ** ptr_header, int * length, int * extended_length, int extend);


#endif
