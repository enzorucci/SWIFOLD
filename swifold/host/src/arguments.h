#ifndef ARGUMENTS_H_INCLUDED
#define ARGUMENTS_H_INCLUDED

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <argp.h>

#define VERSION "1.0"

#define CPU_THREADS 8 
#define MATCH 1
#define MISMATCH 3
#define OPEN_GAP 5
#define EXTEND_GAP 2
#define BLOCK_WIDTH 512
#define DEBUG 1
#define NUM_DEVICES 1
#define MAX_NUM_DEVICES 16

// Arguments parsing
void program_arguments_processing (int argc, char * argv[]);
static int parse_opt (int key, char *arg, struct argp_state *state);

// Global options
extern char * target_filename, * query_filename;
extern int cpu_threads, open_gap, extend_gap, match, mismatch;
extern unsigned int max_num_devices, num_devices;

#endif
