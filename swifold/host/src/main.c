#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include "CL/opencl.h"
#include "AOCL_Utils.h"
#include "arguments.h"
#include "sequences.h"
#include "utils.h"

using namespace aocl_utils;

// OpenCL runtime configuration
#define STRING_BUFFER_LEN 1024
#define PRECOMPILED_BINARY "sw"

// OpenCL runtime configuration
cl_platform_id platform = NULL;
scoped_array<cl_device_id> devices; // num_devices elements
cl_context context = NULL;
scoped_array<cl_command_queue> queues; // num_devices elements
cl_program program = NULL;
scoped_array<cl_kernel> kernels; // num_devices elements
scoped_array<cl_event> kernel_events;

// Function prototypes
bool init();
void cleanup();
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name);
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name);
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name);
static void device_info_string( cl_device_id device, cl_device_info param, const char* name);
static void display_device_info( cl_device_id device );
void checkErr(cl_int err, const char * name);

// Program options
char *target_filename=NULL, * query_filename=NULL;
int match=MATCH, mismatch=MISMATCH, open_gap=OPEN_GAP, extend_gap=EXTEND_GAP, cpu_threads=CPU_THREADS;
unsigned int num_devices=NUM_DEVICES, max_num_devices=MAX_NUM_DEVICES;


// Entry point.
int main(int argc, char * argv[]) {

	char * a, * a_header, *b, * b_header;
	int m, n, ext_m, ext_n, nbb, open_extend_gap, num_active_devices, * scores[MAX_NUM_DEVICES], score=0;
	int i, ii, j, jj, flag=0, d;
    time_t current_time = time(NULL);
	double workTime, tick;
	// CL vars
	cl_int status;
	cl_mem cl_a[MAX_NUM_DEVICES], cl_b[MAX_NUM_DEVICES], cl_scores[MAX_NUM_DEVICES], cl_prev_maxRow[MAX_NUM_DEVICES], cl_next_maxRow[MAX_NUM_DEVICES];
	cl_mem cl_prev_lastCol[MAX_NUM_DEVICES], cl_next_lastCol[MAX_NUM_DEVICES];
	 // Configure work set over which the kernel will execute
	size_t wgSize[3] = {1, 1, 1};
	size_t gSize[3] = {1, 1, 1};

    /* Process program arguments */
    program_arguments_processing(argc,argv);

	// init device and create kernel
	if(!init()) {
		return -1;
	}

	// Load query sequence
	load_sequence(query_filename,&a,&a_header,&m,&ext_m,0);

	// Load target sequence
	load_sequence(target_filename,&b,&b_header,&n,&ext_n,1);

	// Calculate nba and nbb
	nbb = ext_n / BLOCK_WIDTH;

	// Print database search information
	printf("\nSWIFOLD v%s\n \n",VERSION);
	printf("Input sequence:\t\t\t%s (%ld residues) \n",a_header,m);
	printf("Target sequence:\t\t%s (%ld residues) \n",b_header,n);
	printf("Match score:\t\t\t%d\n",match);
	printf("Mismatch penalty:\t\t%d\n",mismatch);
	printf("Gap open penalty:\t\t%d\n",open_gap);
	printf("Gap extend penalty:\t\t%d\n",extend_gap);
/*
	// allocate scores buffer
	for (d=0; d<num_devices ; d++)
		posix_memalign((void**)&scores[d], AOCL_ALIGNMENT, sizeof(int));
*/
		
	tick = dwalltime();

	// Create buffers in device 
	for (d=0; d<num_devices ; d++) {

		// input sequence
		cl_a[d] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, m* sizeof(char), a, &status);
		checkErr(status,"clCreateBuffer cl_a");
		// target sequence
		cl_b[d] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, ext_n* sizeof(char), b, &status);
		checkErr(status,"clCreateBuffer cl_b");
		// buffer that stores the last column of H from the previous block
		cl_prev_lastCol[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
		checkErr(status,"clCreateBuffer cl_prev_lastCol");
		// buffer that stores the last column of H from current block
		cl_next_lastCol[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
		checkErr(status,"clCreateBuffer cl_next_lastCol");
		// buffer that stores the last column of E from the previous block
		cl_prev_maxRow[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
		checkErr(status,"clCreateBuffer cl_prev_maxRow");
		// buffer that stores the last column of E from current block
		cl_next_maxRow[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
		checkErr(status,"clCreateBuffer cl_next_maxRow");
		// max score
		cl_scores[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), &score, &status);
		checkErr(status,"clCreateBuffer cl_scores");

		// Set the kernel arguments
		status = clSetKernelArg(kernels[d], 0, sizeof(cl_mem), &cl_a[d]);
		checkError(status, "Failed to set kernels[d] arg 0");

		status = clSetKernelArg(kernels[d], 1, sizeof(int), &m);
		checkError(status, "Failed to set kernels[d] arg 1");

		status = clSetKernelArg(kernels[d], 2, sizeof(cl_mem), &cl_b[d]);
		checkError(status, "Failed to set kernels[d] arg 2");

		// number of vertical blocks
		status = clSetKernelArg(kernels[d], 3, sizeof(int), &nbb);
		checkError(status, "Failed to set kernels[d] arg 3");

		open_extend_gap = open_gap + extend_gap;
		status = clSetKernelArg(kernels[d], 4, sizeof(int), &open_extend_gap);
		checkError(status, "Failed to set kernels[d] arg 4");

		status = clSetKernelArg(kernels[d], 5, sizeof(int), &extend_gap);
		checkError(status, "Failed to set kernels[d] arg 5");

		status = clSetKernelArg(kernels[d], 6, sizeof(int), &match);
		checkError(status, "Failed to set kernels[d] arg 6");

		status = clSetKernelArg(kernels[d], 7, sizeof(int), &mismatch);
		checkError(status, "Failed to set kernels[d] arg 7");

		status = clSetKernelArg(kernels[d], 12, sizeof(cl_mem), &cl_scores[d]);
		checkError(status, "Failed to set kernels[d] arg 12");

	}


	for(d = 0; d < num_devices; d++) {

		for (jj=0; jj<nbb ; jj++) {

			// swap buffers lastCol and maxRow 
			if (jj % 2 == 0){

				status = clSetKernelArg(kernels[d], 8, sizeof(cl_mem), &cl_prev_lastCol[d]);
				checkError(status, "Failed to set kernels[d] arg 8");

				status = clSetKernelArg(kernels[d], 9, sizeof(cl_mem), &cl_next_lastCol[d]);
				checkError(status, "Failed to set kernels[d] arg 9");

				status = clSetKernelArg(kernels[d], 10, sizeof(cl_mem), &cl_prev_maxRow[d]);
				checkError(status, "Failed to set kernels[d] arg 10");

				status = clSetKernelArg(kernels[d], 11, sizeof(cl_mem), &cl_next_maxRow[d]);
				checkError(status, "Failed to set kernels[d] arg 11");
			} else {
				status = clSetKernelArg(kernels[d], 9, sizeof(cl_mem), &cl_prev_lastCol[d]);
				checkError(status, "Failed to set kernels[d] arg 9");

				status = clSetKernelArg(kernels[d], 8, sizeof(cl_mem), &cl_next_lastCol[d]);
				checkError(status, "Failed to set kernels[d] arg 8");

				status = clSetKernelArg(kernels[d], 11, sizeof(cl_mem), &cl_prev_maxRow[d]);
				checkError(status, "Failed to set kernels[d] arg 11");

				status = clSetKernelArg(kernels[d], 10, sizeof(cl_mem), &cl_next_maxRow[d]);
				checkError(status, "Failed to set kernels[d] arg 10");
			}

			// vertical block number
			status = clSetKernelArg(kernels[d], 13, sizeof(int), &jj);
			checkError(status, "Failed to set kernels[d] arg 13");

			// Launch the kernel
			status = clEnqueueNDRangeKernel(queues[d], kernels[d], 1, NULL, gSize, wgSize, 0, NULL, NULL);
			checkError(status, "Failed to launch kernel");

			status = clFinish(queues[d]);
			checkError(status, "Failed to finish");

		}

	}

	for(d = 0; d < num_devices; d++) {
		// Copy alignment score to host array 
		status = clEnqueueReadBuffer(queues[d], cl_scores[d], CL_TRUE, 0, sizeof(int), &score, 0, NULL, NULL);
		checkErr(status,"clEnqueueReadBuffer: Couldn't read cl_scores buffer");
	}

	workTime = dwalltime() - tick;

	// Wait for command queue to complete pending events
	for (d=0; d<num_devices ; d++)
		status = clFinish(queues[d]);
	checkError(status, "Failed to finish");

	// Free allocated memory
	free(a);
	free(b);

	printf("\nScore:\t\t\t\t%d\n",score);
	printf("\nSearch date:\t\t\t%s",ctime(&current_time));
	printf("Search time:\t\t\t%lf seconds\n",workTime);
	printf("Search speed:\t\t\t%.2lf GCUPS\n",((long int)n*(long int)m) / (workTime*1000000000));

	// free FPGA resources
	for (d=0; d<num_devices ; d++){
		clReleaseMemObject(cl_a[d]);
		clReleaseMemObject(cl_b[d]);
		clReleaseMemObject(cl_prev_maxRow[d]);
		clReleaseMemObject(cl_next_maxRow[d]);
		clReleaseMemObject(cl_prev_lastCol[d]);
		clReleaseMemObject(cl_next_lastCol[d]);
		clReleaseMemObject(cl_scores[d]);
	}

	// Free the resources allocated
	cleanup();

	return 0;
}

/////// HELPER FUNCTIONS ///////

bool init() {
  cl_int status;

  if(!setCwdToExeDir()) {
    return false;
  }

  // Get the OpenCL platform.
  platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
  if(platform == NULL) {
    printf("ERROR: Unable to find Intel(R) FPGA OpenCL platform.\n");
    return false;
  }

  // User-visible output - Platform information
/*  {
    char char_buffer[STRING_BUFFER_LEN]; 
    printf("Querying platform for info:\n");
    printf("==========================\n");
    clGetPlatformInfo(platform, CL_PLATFORM_NAME, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n", "CL_PLATFORM_NAME", char_buffer);
    clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n", "CL_PLATFORM_VENDOR ", char_buffer);
    clGetPlatformInfo(platform, CL_PLATFORM_VERSION, STRING_BUFFER_LEN, char_buffer, NULL);
    printf("%-40s = %s\n\n", "CL_PLATFORM_VERSION ", char_buffer);
  }
*/
  // Query the available OpenCL devices.
  devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &max_num_devices));

  // Display some device information.
//  display_device_info(devices[0]);

  // Create the context.
  context = clCreateContext(NULL, num_devices, devices, NULL, NULL, &status);
  checkError(status, "Failed to create context");

  queues.reset(num_devices);
  kernels.reset(num_devices);
  kernel_events.reset(num_devices*20);

  // Create the command queue.
  for (int i = 0; i<num_devices ; i++){
	  queues[i] = clCreateCommandQueue(context, devices[i], CL_QUEUE_PROFILING_ENABLE, &status);
	  checkError(status, "Failed to create command queue");
  }

  // Create the program.
  std::string binary_file = getBoardBinaryFile(PRECOMPILED_BINARY, devices[0]);
//  printf("Using AOCX: %s\n", binary_file.c_str());
  program = createProgramFromBinary(context, binary_file.c_str(), devices, num_devices);

  // Build the program that was just created.
  status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
  checkError(status, "Failed to build program");

  for (int i = 0; i<num_devices ; i++){
	  // Create the kernel - name passed in here must match kernel name in the
	  // original CL file, that was compiled into an AOCX file using the AOC tool
	  const char * kernel_name = "sw";  // Kernel name, as defined in the CL file
	  kernels[i] = clCreateKernel(program, kernel_name, &status);
	  checkError(status, "Failed to create kernel");
  }

  return true;
}

// Free the resources allocated during initialization
void cleanup() {
  for(unsigned i = 0; i < num_devices; ++i) {
    if(kernels && kernels[i]) {
      clReleaseKernel(kernels[i]);
    }
    if(queues && queues[i]) {
      clReleaseCommandQueue(queues[i]);
    }
  }
  if(program) {
    clReleaseProgram(program);
  }
  if(context) {
    clReleaseContext(context);
  }
}

// Helper functions to display parameters returned by OpenCL queries
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name) {
   cl_ulong a;
   clGetDeviceInfo(device, param, sizeof(cl_ulong), &a, NULL);
   printf("%-40s = %lu\n", name, a);
}
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name) {
   cl_uint a;
   clGetDeviceInfo(device, param, sizeof(cl_uint), &a, NULL);
   printf("%-40s = %u\n", name, a);
}
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name) {
   cl_bool a;
   clGetDeviceInfo(device, param, sizeof(cl_bool), &a, NULL);
   printf("%-40s = %s\n", name, (a?"true":"false"));
}
static void device_info_string( cl_device_id device, cl_device_info param, const char* name) {
   char a[STRING_BUFFER_LEN]; 
   clGetDeviceInfo(device, param, STRING_BUFFER_LEN, &a, NULL);
   printf("%-40s = %s\n", name, a);
}

// Query and display OpenCL information on device and runtime environment
static void display_device_info( cl_device_id device ) {

   printf("Querying device for info:\n");
   printf("========================\n");
   device_info_string(device, CL_DEVICE_NAME, "CL_DEVICE_NAME");
   device_info_string(device, CL_DEVICE_VENDOR, "CL_DEVICE_VENDOR");
   device_info_uint(device, CL_DEVICE_VENDOR_ID, "CL_DEVICE_VENDOR_ID");
   device_info_string(device, CL_DEVICE_VERSION, "CL_DEVICE_VERSION");
   device_info_string(device, CL_DRIVER_VERSION, "CL_DRIVER_VERSION");
   device_info_uint(device, CL_DEVICE_ADDRESS_BITS, "CL_DEVICE_ADDRESS_BITS");
   device_info_bool(device, CL_DEVICE_AVAILABLE, "CL_DEVICE_AVAILABLE");
   device_info_bool(device, CL_DEVICE_ENDIAN_LITTLE, "CL_DEVICE_ENDIAN_LITTLE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_SIZE, "CL_DEVICE_GLOBAL_MEM_SIZE");
   device_info_bool(device, CL_DEVICE_IMAGE_SUPPORT, "CL_DEVICE_IMAGE_SUPPORT");
   device_info_ulong(device, CL_DEVICE_LOCAL_MEM_SIZE, "CL_DEVICE_LOCAL_MEM_SIZE");
   device_info_ulong(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, "CL_DEVICE_MAX_CLOCK_FREQUENCY");
   device_info_ulong(device, CL_DEVICE_MAX_COMPUTE_UNITS, "CL_DEVICE_MAX_COMPUTE_UNITS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_ARGS, "CL_DEVICE_MAX_CONSTANT_ARGS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE");
   device_info_uint(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, "CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE");

   {
      cl_command_queue_properties ccp;
      clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(cl_command_queue_properties), &ccp, NULL);
      printf("%-40s = %s\n", "Command queue out of order? ", ((ccp & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)?"true":"false"));
      printf("%-40s = %s\n", "Command queue profiling enabled? ", ((ccp & CL_QUEUE_PROFILING_ENABLE)?"true":"false"));
   }
}

/* Error checking */
void checkErr(cl_int err, const char * name)
{
	if (err != CL_SUCCESS) {
		printf("\n ERROR (%d): %s\n",err,name);
		exit(EXIT_FAILURE);
	}
}
