#include "OpenCLManager.h"


/**
	Name: OpenCLManager
	Parameters: None
	Return values: None
	Description: Simple constructorm that initializes to default values
**/
OpenCLManager::OpenCLManager()
{
	this->platform = 0;
	status = 0;
	blockSize = 8; // Default size of 8. If one wants a bigger size for graphic card this mus be set manually
}

/**
	Name: printAllDevicesInfo
	Parameters: None
	Return values: 0 if everything was allright
	Description:	Function from OpenCLCookbook that lists all platforms and devices in the system
					It uses internal variables for platforms and devices which are freed after execution
**/
int OpenCLManager::printAllDevicesInfo()
{
	unsigned int i, j;
    char* value;
    size_t valueSize;
    cl_uint platformCount;
    cl_platform_id* platforms;
    cl_uint deviceCount;
    cl_device_id* devices;
    cl_uint maxComputeUnits;
 
    // get all platforms
    clGetPlatformIDs(0, NULL, &platformCount); // find out how many platforms are in the systems
    platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount); // allocate a buffer for all these platforms
    clGetPlatformIDs(platformCount, platforms, NULL); // read the informations about the platforms and write them to plaforms-fields.
 
	// run over all platforms and list the devices in them
    for (i = 0; i < platformCount; i++) {
        // get all devices for that platform
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount); // find out how many devices are associated with this platform
        devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount); // allocate a buffer for a apropriate count of devices
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL); // read the information about each device

		// for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {
            // print device name
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
            printf("%d. Device: %s\n", j+1, value);
            free(value);

			// print hardware device version
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
            printf(" %d.%d Hardware version: %s\n", j+1, 1, value);
            free(value);

			// print software driver version
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
            printf(" %d.%d Software version: %s\n", j+1, 2, value);
            free(value);
 
            // print c version supported by compiler for device
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            printf(" %d.%d OpenCL C version: %s\n", j+1, 3, value);
            free(value);

            // print parallel compute units
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(maxComputeUnits), &maxComputeUnits, NULL);
            printf(" %d.%d Parallel compute units: %d\n", j+1, 4, maxComputeUnits);

			// print some other values about the device-memory system:
			cl_ulong tmp=0;
			clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(cl_ulong), (void*)&tmp,NULL);
			printf(" %d.%d Global mem cache size: %lu bytes.\n", j+1, 4,  tmp);

			clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, sizeof(cl_ulong), (void*)&tmp,NULL);
			printf(" %d.%d Global mem cacheline size: %lu bytes.\n", j+1, 4,  tmp);

			clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), (void*)&tmp,NULL);
			printf(" %d.%d Global mem size: %lu bytes.\n", j+1, 4,  tmp);

			clGetDeviceInfo(devices[j],CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), (void*)&tmp,NULL);
			printf(" %d.%d Local mem size: %lu bytes.\n", j+1, 4,  tmp);

			int stop=0;



        }
        free(devices);
    }
    free(platforms);
    return 0;
}

int OpenCLManager::initializeOpenCLDevices() {
	cl_int status = 0;
	size_t deviceListSize;
	bool untestedPlattform = true;
	std::string vendor;

	// Search for all Platforms (AMD, Intel, NVidia, ...) and detect how many there are
	cl_uint numPlatforms;
	platform = NULL;
	status = clGetPlatformIDs(0,NULL, &numPlatforms);

	if(status != CL_SUCCESS)
    {
            std::cout<<
                    "clGetPlatformIDs failed -> count of platforms";
            return 1;
    }

	// If there are any platform find detailed infos for them and save them in the platforms variable
	if (0 < numPlatforms)
	{
		platforms = new cl_platform_id[numPlatforms];
		status = clGetPlatformIDs(numPlatforms, platforms, NULL);
        if(status != CL_SUCCESS)
        {
                std::cout<<
                        "clGetPlatformIDs failed -> plattform IDs";
                return 1;
        }

        for (unsigned i = 0; i < numPlatforms; ++i)
        {
            char pbuf[100];
            status = clGetPlatformInfo(platforms[i],
                                                                CL_PLATFORM_NAME,
                                                                sizeof(pbuf),
                                                                pbuf,
                                                                NULL);
            if(status != CL_SUCCESS)
            {
                    std::cout<<
                            "clGetPlatformInfo failed";
                    return 1;
            }

            platform = platforms[i];

			// If there is a platform, that contains a graphics card, then take it!
			cl_uint tmpdevcnt = 0;
			clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &tmpdevcnt);
			if(tmpdevcnt > 0) {
				vendor = pbuf;
				break;
			}
            vendor = pbuf;
        }
		delete[] platforms;
	}
    
    // If we could find our platform, use it. Otherwise pass a NULL and get whatever the
    // implementation thinks we should be using.
    cl_context_properties cps[3] =
    {
            CL_CONTEXT_PLATFORM,
            (cl_context_properties)platform,
            0
    };
   //  Use NULL for backward compatibility
    cl_context_properties* cprops = (NULL == platform) ? NULL : cps;

	// Test for GPU
    context = clCreateContextFromType(cprops,
            CL_DEVICE_TYPE_GPU,
            NULL,
            NULL,
            &status);

	// If there is no GPU test for Xeon Phi
    if(status != CL_SUCCESS )
    {
            std::cout << "No Graphic Card found, searching for Xeon Phi ..." << std::endl;
            context = clCreateContextFromType(
                    cprops,
                    CL_DEVICE_TYPE_ACCELERATOR,
                    NULL,
                    NULL,
                    &status);
    }

	// If GPU is not available, then use CPU
    if(status != CL_SUCCESS )
    {
            std::cout << "Unsupported/no GPU and no Xeon Phi device; falling back to CPU ..." << std::endl;
            context = clCreateContextFromType(
                    cprops,
                    CL_DEVICE_TYPE_CPU,
                    NULL,
                    NULL,
                    &status);
            if(status != CL_SUCCESS)
            {
                    std::cout<<"Error: Creating Context. (clCreateContextFromType)\n";
                    return 1;
            }
    }

    // First, get the size of device lista data
    status = clGetContextInfo(context,
                            CL_CONTEXT_DEVICES,
                            0,
                            NULL,
                            &deviceListSize);
    if(status != CL_SUCCESS)
    {
            std::cout<<
                    "Error: Getting Context Info \
                    (device lista size, clGetContextInfo)\n";
            return 1;
    }

    /////////////////////////////////////////////////////////////////
    // Detect OpenCL devices
    /////////////////////////////////////////////////////////////////
    devices = (cl_device_id *)malloc(deviceListSize);
    if(devices == 0)
    {
            std::cout<<"Error: No devices found.\n";
            return 1;
    }

    /* Now, get the device lista data */
    status = clGetContextInfo(
                                context,
                                CL_CONTEXT_DEVICES,
                                deviceListSize,
                                devices,
                                NULL);
    if(status != CL_SUCCESS)
    {
            std::cout<<
                    "Error: Getting Context Info \
                    (device lista, clGetContextInfo)\n";
            return 1;
    }

    /////////////////////////////////////////////////////////////////
    // Create an OpenCL command queue
    /////////////////////////////////////////////////////////////////
    commandQueue = clCreateCommandQueue(
                                        context,
                                        devices[0],
                                        CL_QUEUE_PROFILING_ENABLE,
                                        &status);
    if(status != CL_SUCCESS)
    {
            std::cout<<"Creating Command Queue. (clCreateCommandQueue)\n";
            return 1;
    }

	return 0;
}

int OpenCLManager::prepareKernel(const char* file, const char*  kernelname) {

	std::string sourceStr;
    sourceStr.append(convertToString(file));

    const char * source    = sourceStr.c_str();
    size_t sourceSize[]    = { strlen(source) };

    program = clCreateProgramWithSource(
                                context,
                                1,
                                &source,
                                sourceSize,
                                &status);
    if(status != CL_SUCCESS)
    {
        std::cout<<
                        "Error: Loading Binary into cl_program \
                        (clCreateProgramWithBinary)\n";
        return 1;
    }

    status = clBuildProgram(program, 1, devices, NULL, NULL, NULL); // evt. 2. para auf 0
    if(status != CL_SUCCESS)
    {
            std::cout<<"Error: Building Program (clBuildProgram)\n";
            //print kernel compilation error
            char programLog[8192];

            status = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 8192, programLog, 0);

            std::cout<<programLog<<std::endl;
            return 1;
    }

    char device_name[128];
    status = clGetDeviceInfo(devices[0], CL_DEVICE_NAME, 128, device_name, NULL);
    if (status!=CL_SUCCESS)
    {
        printf("ERROR: Failed to get device information (device name)...\n");
        return -1;
    }
    printf("Using device %s...\n", device_name);

	kernel = clCreateKernel(program, kernelname, &status);

    if(status != CL_SUCCESS)
    {
            std::cout<<"Error: Creating Kernel from program. (clCreateKernel)\n";
            return 1;
    }
	return 0;
}

int OpenCLManager::checkNDRangeConfig(size_t* global, size_t* localThreads, int NDRangeDimension) {
    cl_uint maxDims;
	cl_uint maxComputeUnits;
	size_t maxWorkGroupSize;                                                                                                                                            
    size_t maxWorkItemSizes[3];                                                                                                                                         
                                                                                                                                                                            
    cl_ulong usedLocalMemory, totalLocalMemory, availableLocalMemory, neededLocalMemory;     
                                                                                                                                                          
    // Query device capabilities. Maximum                                                                                                                                
    // work item dimensions and the maximmum                                                                                                                             
    // work item sizes                                                                                                                                                   
                                                                                                                                                                
    clGetDeviceInfo(                                                                                                                                                    
            devices[0],                                                                                                                                                 
            CL_DEVICE_MAX_WORK_GROUP_SIZE,                                                                                                                              
            sizeof(size_t),                                                                                                                                             
            (void*)&maxWorkGroupSize,                                                                                                                                   
            NULL);                                                                                                                                                      
                                                                                                                                                                            
    clGetDeviceInfo(                                                                                                                                                    
            devices[0],                                                                                                                                                 
            CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,                                                                                                                         
            sizeof(cl_uint),                                                                                                                                            
            (void*)&maxDims,                                                                                                                                            
            NULL);

    clGetDeviceInfo(                                                                                                                                                    
            devices[0],                                                                                                                                                 
            CL_DEVICE_MAX_WORK_ITEM_SIZES,                                                                                                                              
            sizeof(size_t)*maxDims,                                                                                                                                     
            (void*)maxWorkItemSizes,                                                                                                                                    
            NULL);                                                                                                                                                      
                                                                                                                                                                            
    clGetDeviceInfo(                                                                                                                                                   
            devices[0],                                                                                                                                         
            CL_DEVICE_LOCAL_MEM_SIZE,                                                                                                                           
            sizeof(cl_ulong),                                                                                                                                   
            (void *)&totalLocalMemory,                                                                                                                          
            NULL);  

    clGetDeviceInfo(                                                                                                                                                   
            devices[0],                                                                                                                                         
            CL_DEVICE_MAX_COMPUTE_UNITS,                                                                                                                           
            sizeof(cl_uint),                                                                                                                                   
            (void *)&maxComputeUnits,                                                                                                                          
            NULL); 

	printf("Device has %lu bytes of local memory.\n", totalLocalMemory);
	printf("Device has %i compute units.\n", maxComputeUnits);

	if(NDRangeDimension==3)
	{
		if(	(localThreads[0] > maxWorkItemSizes[0]) || 
			(localThreads[1] > maxWorkItemSizes[1]) || 
			(localThreads[2] > maxWorkItemSizes[2]) ||
			(localThreads[2] * localThreads[1] * localThreads[0] > maxWorkGroupSize))
		{
			std::cout<<"Unsupported: Device does not support requested number of work items.";                                                                              
			return 1;    	
		}
	}

	if(NDRangeDimension==2) {
		if(	(localThreads[0] > maxWorkItemSizes[0]) || 
			(localThreads[1] > maxWorkItemSizes[1]) || 
			(localThreads[1] * localThreads[0] > maxWorkGroupSize))                                                                                                        
		{                                                                                                                                                                   
 			std::cout<<"Unsupported: Device does not support requested number of work items.";                                                                              
			return 1;    	                                                                                                                                            
		}   
	}

	if(NDRangeDimension==1) {
		if(	localThreads[0] > maxWorkItemSizes[0] ||                                                                                 
			localThreads[0] > maxWorkGroupSize)                                                                                                         
		{                                                                                                                                                                   
				std::cout<<"Unsupported: Device does not support requested number of work items.";                                                                              
				return 1;                                                                                                                                                 
		} 
	}

    clGetKernelWorkGroupInfo(                                                                                                                                           
                                    kernel,                                                                                                                             
                                    devices[0],                                                                                                                         
                                    CL_KERNEL_LOCAL_MEM_SIZE,                                                                                                           
                                    sizeof(cl_ulong),                                                                                                                   
                                    &usedLocalMemory,                                                                                                                   
                                    NULL);                                                                                                                              
                                                                                                                                                                            
    availableLocalMemory = totalLocalMemory - usedLocalMemory;                                                                                                          
                                                                                                                                                                            
    neededLocalMemory    = 16*blockSize*blockSize*sizeof(cl_float);                                                                                                     
                                                                                                                                                                            
    if(neededLocalMemory > availableLocalMemory)                                                                                                                        
    {                                                                                                                                                                   
            std::cout << "Unsupported: Insufficient local memory on device." << std::endl;     
			return 1;
    }

	return 0;
}

double OpenCLManager::printProfilingInfo(cl_event* events) {
		cl_ulong kernelsStartTime = 0;                                                                                                                                      
		cl_ulong kernelsEndTime = 0;                                                                                                                                        
		status = clGetEventProfilingInfo(                                                                                                                           
									events[0],                                                                                                                         
									CL_PROFILING_COMMAND_START,                                                                                                        
									sizeof(cl_ulong),                                                                                                                  
									&kernelsStartTime,                                                                                                                 
									NULL);                                                                                                                             
		if(status != CL_SUCCESS) {
			printf("clGetEventProfilingInfo failed.");     
			return 1;
		}
                                                                                                                                                                            
						                                                                                             
                                                                         
                                                                                                                                                                            
		status = clGetEventProfilingInfo(                                                                                                                           
									events[0],                                                                                                                         
									CL_PROFILING_COMMAND_END,                                                                                                          
									sizeof(cl_ulong),                                                                                                                  
									&kernelsEndTime,                                                                                                                   
									NULL);                                                                                                                             
		if(status != CL_SUCCESS) {
			printf("clGetEventProfilingInfo failed.");
			return 1;
		}                                                      
                                                                                                                                                                            
		// Compute total time (also convert from nanoseconds to seconds)                                                                                            
		double totalKernelTime = (double)(kernelsEndTime - kernelsStartTime)/1e9;
		//printf("OpenCL timing for kernel GPU Tracer %f sec \n",totalKernelTime);  
	return totalKernelTime;
}

std::string OpenCLManager::convertToString(const char *filename)
{
    size_t size;
    char*  str;
    std::string s;

    std::fstream f(filename, (std::fstream::in | std::fstream::binary));

    if(f.is_open())
    {
            size_t fileSize;
            f.seekg(0, std::fstream::end);
            size = fileSize = size_t(f.tellg());
            f.seekg(0, std::fstream::beg);

            str = new char[size+1];
            if(!str)
            {
                    f.close();
                    return NULL;
            }

            f.read(str, fileSize);
            f.close();
            str[size] = '\0';

            s = str;

            return s;
    }
    return "";
}