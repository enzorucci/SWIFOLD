#include "arguments.h"

const char *argp_program_bug_address =
  "<erucci@lidi.info.unlp.edu.ar>";

/* Program documentation. */
static char doc[] =
  "\nSWIFOLD";

void program_arguments_processing (int argc, char * argv[]) {

	int i, result, arg_count=4;

	struct argp_option options[] =	{
		{ 0, 0, 0, 0, "SWIFOLD execution", 1},
		{ "input", 'i', "<string>", 0,  "Query sequence filename (FASTA format). [REQUIRED]", 1},
		{ "target", 'j', "<string>", 0, "Target sequence filename (FASTA format). [REQUIRED]", 1},
		{ "match", 'M', "<integer>", 0, "Match score (default: 1).", 2},
		{ "mismatch", 'X', "<integer>", 0, "Mismatch penalty (default: 3).", 2},
		{ "gap_open", 'g', "<integer>", 0, "Gap open penalty (default: 5).", 2},
		{ "gap_extend", 'e', "<integer>", 0, "Gap extend penalty (default: 2).", 2},
		{ "cpu_threads", 'c', "<integer>", 0, "Number of CPU threads.", 2},
		{ 0 }
	};

	struct argp argp = { options, parse_opt, 0, doc};
	result = argp_parse (&argp, argc, argv, 0, 0, &arg_count); 
}

static int parse_opt (int key, char *arg, struct argp_state *state) {

	int *arg_count = (int *) state->input;

	switch(key) {
		case 'i':
			query_filename = arg;
			break;
		case 'j':
			target_filename = arg;
			break;
		case 'M':
			match = atoi(arg);
			break;
		case 'X':
			mismatch = atoi(arg);
			break;
		case 'g':
			open_gap = atoi(arg);
			break;
		case 'e':
			extend_gap = atoi(arg);
			break;
		case 'c':
			cpu_threads = atoi(arg);
			if (cpu_threads < 0)
				argp_failure (state, 1, 0, "The number of CPU threads must be greater than 0.");
			break;
		case ARGP_KEY_END:
			if (*arg_count == 1)
				argp_failure (state, 1, 0, "Missing options");
			if (query_filename == NULL)
				argp_failure (state, 1, 0, "Target sequence filename is required");
			if (target_filename == NULL)
				argp_failure (state, 1, 0, "Query sequence filename is required");
/*	    default:
			return ARGP_ERR_UNKNOWN;
*/	}

	return 0;
} 

