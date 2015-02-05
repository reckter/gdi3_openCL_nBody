#include <CL/cl.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>

class OpenCLManager {
public:
	OpenCLManager();

	int printAllDevicesInfo();
	int initializeOpenCLDevices();
	int prepareKernel(const char* file, const char*  kernelname);
	int checkNDRangeConfig(size_t* global, size_t* local, int NDRangeDimension);
	double printProfilingInfo(cl_event* events);

	cl_device_id* devices;
	cl_context context;
	// This program uses only one kernel and this serves as a handle to it
	cl_kernel  kernel;
	cl_command_queue    commandQueue;

private:
	std::string convertToString(const char *filename);
	cl_int status;
	cl_platform_id          platform;
	cl_platform_id*			platforms;
	cl_program program;
	cl_uint blockSize;
};