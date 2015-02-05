#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cmath>
#include <vector>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "OpenCLManager.h"

using namespace std;

OpenCLManager* oclm;

// global vars
int SIM_TYPE = 0; // termines which kernel will be loaded
int nthread = 64; // Number of local threads
int NDRangeDimension = 1; // 	
int nstep; // How many simulation steps will be done?
int nparticle; // How many particles
float dt; // Integration step size
float eps; // mininum distance between two particles

float* pos1; // Position Array of particles. Also contains mass information.
float* pos1g; cl_mem pos1g_buf; // Buffers for GPU solution

float* pos2; // Array for temporary position updates
float* pos2g; cl_mem pos2g_buf; // Buffers for GPU solution

float* vel; // Array of velocities in all three dimensions
float* velg; cl_mem velg_buf; // Buffers for GPU solution

float* ref_sol; // GPU reference solution, that is read from file

int is_even; cl_mem is_even_buf; // Inidcates whether pos1 or pos2 is the active buffer

// OpenCL related vars
cl_uint platformCount;

cl_uint deviceCount;
cl_device_id* devices;
cl_uint maxComputeUnits;

// Function Prototypes
int cleanCL();
int setupCL();
int runCL();

int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && !(x & (x - 1)));
}

void readSolutionFromFile() {
	std::ifstream myfile ("NBodySolutionDat.txt");
	std::string act_line;
	int i = 0;

	if (myfile.is_open()) {
		std::getline(myfile, act_line);

		std::istringstream iss(act_line);
		int particles;
		iss >> particles;

		if(particles != nparticle) {
			cout << "ERROR: Particel cound of input and reference-solution file does not match!";
			return;
		}

		ref_sol = new float[particles*4];

		while (std::getline(myfile, act_line))
		{
			std::istringstream iss(act_line);
			float x,y,z;
			if (!(iss >> x >> y >> z)) { break; } // error

			ref_sol[4*i + 0] = x;
			ref_sol[4*i + 1] = y;
			ref_sol[4*i + 2] = z;
			ref_sol[4*i + 3] = 0.0f;
			i++;
		}
		myfile.close();
	}
}

void readExampleBodyValsFile() {
	std::ifstream myfile ("NBodyInputDat.txt");
	std::string act_line;
	int i = 0;

	if (myfile.is_open()) {
		std::getline(myfile, act_line);
		std::istringstream iss(act_line);
		iss >> nparticle >> nstep;

		pos1 = new float[nparticle*4];
		pos1g = new float[nparticle*4];

		pos2 = new float[nparticle*4];
		pos2g = new float[nparticle*4];

		vel = new float[nparticle * 4];
		velg = new float[nparticle * 4];

		while (std::getline(myfile, act_line))
		{
			std::istringstream iss(act_line);
			float x,y,z,m;
			if (!(iss >> x >> y >> z >> m)) { break; } // error

			pos1[4*i + 0] = x; pos2[4*i + 0] = 0.0f;
			pos1[4*i + 1] = y; pos2[4*i + 1] = 0.0f;
			pos1[4*i + 2] = z; pos2[4*i + 2] = 0.0f;
			pos1[4*i + 3] = m; pos2[4*i + 3] = m;

			pos1g[4*i + 0] = x; pos2g[4*i + 0] = 0.0f;
			pos1g[4*i + 1] = y; pos2g[4*i + 1] = 0.0f;
			pos1g[4*i + 2] = z; pos2g[4*i + 2] = 0.0f;
			pos1g[4*i + 3] = m; pos2g[4*i + 3] = m;

			vel[4*i + 0] = 0.0f; velg[4*i + 0] = 0.0f;
			vel[4*i + 1] = 0.0f; velg[4*i + 1] = 0.0f;
			vel[4*i + 2] = 0.0f; velg[4*i + 2] = 0.0f;
			vel[4*i + 3] = 0.0f; velg[4*i + 3] = 0.0f;
			i++;
		}
		myfile.close();
	}
}

bool evalAgainstFile() {

	bool is_ok = true;
	vector<int> err_row;
	vector<float> err_vals;
	int err_cnt = 0;
	int printcnt = 0;
	int printcnt2 = 0;

	cout << "Checking Solution vector: \n";

	for(int i=0; i<nparticle*3; i++) {
		if(i%4==3)
			continue;

		float oval = ref_sol[i]; // original value from reference solution
		float cval = pos1g[i]; // calculated value from programm
		float gval = pos1g[i];
		float tmp;

		if(abs(oval) > abs(cval)) {
			tmp = oval;
			oval = cval;
			cval = tmp;
		}
		float diff = (cval - oval) / cval * 100.0f;

		if(printcnt2 < 10) {
			cout << i << ": " << oval << "   /   " << cval << ";" << diff << "\n";
			printcnt2++;
		}

		if(abs(diff) > 0.50) {
			// If both are very small then also don't print
			if(!(	abs(oval) < 0.1 && abs(cval) < 0.1))
			{
				is_ok = false;
				err_row.push_back(i-3);
				err_vals.push_back(abs(diff));
				err_vals.push_back(cval);
				err_vals.push_back(oval);
				err_cnt++;
			}
		}
	}

	if(is_ok) {
		std::cout << "Residuum Evalution successfull. No errors found.\n";
	}

	else {
		for(int i=0; i<err_cnt;i++) {
			if(printcnt < 25)
			{
				std::cout	<<  "Error of " << err_vals[i*3+0] << " @R:" << err_row[i] << "   "
							<<  err_vals[i*3+1] << " " << err_vals[i*3+2] << "\n";
				printcnt++;
			}
		}
	}

	std::cout << "\n";
	return is_ok;

}

void executeReference(int nstep, int nparticle, float dt, float eps, float* pos1, float* pos2, float* vel) {
	for(int n=0; n<nstep; n++)
	{
		for(int i=0; i<nparticle; i++) {
			float ax=0.0f;
			float ay=0.0f;
			float az=0.0f;

			for(int j=0; j<nparticle; j++) {
				float dx = pos1[j*4 + 0] - pos1[i*4 + 0];
				float dy = pos1[j*4 + 1] - pos1[i*4 + 1];
				float dz = pos1[j*4 + 2] - pos1[i*4 + 2];

				float invr = 1.0f / sqrt(dx*dx + dy*dy + dz*dz + eps);

				float invr3 = invr*invr*invr;

				float f = pos1[j*4 + 3] * invr3;

				ax += f*dx; /* accumulate the acceleration from gravitational attraction */
				ay += f*dy;
				az += f*dz;
			}

			pos2[i*4 + 0] = pos1[i*4 + 0] + dt*vel[i*4 + 0] + 0.5f*dt*dt*ax;
			pos2[i*4 + 1] = pos1[i*4 + 1] + dt*vel[i*4 + 1] + 0.5f*dt*dt*ay;
			pos2[i*4 + 2] = pos1[i*4 + 2] + dt*vel[i*4 + 2] + 0.5f*dt*dt*az;

			vel[i*4 + 0] +=dt*ax;
			vel[i*4 + 1] +=dt*ay;
			vel[i*4 + 2] +=dt*az;
		}

		for(int i=0; i<nparticle*4; i++) {
			pos1[i] = pos2[i];
		}
	}

	return;
}

int setupCL() {
	cl_int status = 0;	
	cl_context context;
	devices = NULL;

	// Context, CommandQueue und Devices
	oclm = new OpenCLManager();
	oclm->printAllDevicesInfo();
	oclm->initializeOpenCLDevices();

	switch(SIM_TYPE) {
	case 0:
		oclm->prepareKernel("nbody_std.cl", "nbody");
		std::cout<<"Using standard kernel.\n";  
		break;
	case 1:
		oclm->prepareKernel("nbody_vec.cl", "nbody");
		std::cout<<"Using vectorized kernel\n";  
		break;
	case 2:
		oclm->prepareKernel("nbody_mem.cl", "nbody");
		std::cout<<"Using local memory kernel\n";  
		break;
	default:
		oclm->prepareKernel("nbody_std.cl", "nbody");
		std::cout<<"Using standard kernel.\n"; 
		break;
	}	

	context = oclm->context;

    // Initialize buffers 
	// Adapt for every project
	pos2g_buf = clCreateBuffer(                                                                                                                             
                context,                                                                                                                         
                CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,                                                                                                               
				sizeof(cl_float) * (nparticle*4),                                                                                              
				pos2g,                                                                                                                            
                            &status);                                                                                                                        
    if(status != CL_SUCCESS)                                                                                                                                    
    {                                                                                                                                                           
            std::cout<<"Error: clCreateBuffer (pos2g_buf)\n";      
			return 0;
    }

	pos1g_buf = clCreateBuffer(                                                                                                                             
                context,                                                                                                                         
                CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,                                                                                                               
				sizeof(cl_float) * (nparticle*4),                                                                                              
				pos1g,                                                                                                                            
                            &status);                                                                                                                        
    if(status != CL_SUCCESS)                                                                                                                                    
    {                                                                                                                                                           
            std::cout<<"Error: clCreateBuffer (pos1g_buf)\n";      
			return 0;
    }

	velg_buf = clCreateBuffer(                                                                                                                             
                context,                                                                                                                         
                CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,                                                                                                               
				sizeof(cl_float) * (nparticle*4),                                                                                              
				velg,                                                                                                                            
                            &status);                                                                                                                        
    if(status != CL_SUCCESS)                                                                                                                                    
    {                                                                                                                                                           
            std::cout<<"Error: clCreateBuffer (velg_buf)\n";      
			return 0;
    }

	return 0;
}

int runCL() {
    cl_int status;  
	size_t globalThreads[1];                                                                                                                                            
    size_t localThreads[1];
	cl_event events[2]; 

	globalThreads[0] = nparticle;
	localThreads[0] = nthread;

	oclm->checkNDRangeConfig(globalThreads, localThreads, NDRangeDimension);

	cl_kernel kernel = oclm->kernel;

	double elapsed_time = 0.0;
	is_even = 1;


	/***!TODO!: Complete the argument-list of the clSetKernelArg call below.
	This Argument should be the velocity vector!***/
	// velocity input
	status = clSetKernelArg(
            kernel,
            2,
            sizeof(cl_mem),
            (void *)velg_buf);
	if(status != CL_SUCCESS)                                                                                                                                            
	{                                                                                                                                                                   
			std::cout<<"Error: Setting kernel argument. (velg_buf): " << status << "\n";
			return 1;                                                                                                                                                   
	}

	// number of particles
	status = clSetKernelArg(                                                                                                                                           
									kernel,                                                                                                                             
									4,                                                                                                                                  
									sizeof(cl_int),                                                                                                                     
									(void *)&nparticle);                                                                                                            
	if(status != CL_SUCCESS)                                                                                                                                            
	{                                                                                                                                                                   
			std::cout<<"Error: Setting kernel argument. (nparticle)\n";                                                                                                    
			return 1;                                                                                                                                                   
	}

	// step-width during integration
	status = clSetKernelArg(                                                                                                                                           
									kernel,                                                                                                                             
									5,                                                                                                                                  
									sizeof(cl_int),                                                                                                                     
									(void *)&dt);                                                                                                            
	if(status != CL_SUCCESS)                                                                                                                                            
	{                                                                                                                                                                   
			std::cout<<"Error: Setting kernel argument. (dt)\n";                                                                                                    
			return 1;                                                                                                                                                   
	}

	// step-width during integration
	status = clSetKernelArg(                                                                                                                                           
									kernel,                                                                                                                             
									6,                                                                                                                                  
									sizeof(cl_int),                                                                                                                     
									(void *)&eps);                                                                                                            
	if(status != CL_SUCCESS)                                                                                                                                            
	{                                                                                                                                                                   
			std::cout<<"Error: Setting kernel argument. (eps)\n";                                                                                                    
			return 1;                                                                                                                                                   
	}

	// is-even trigger
	status = clSetKernelArg(                                                                                                                                           
									kernel,                                                                                                                             
									3,                                                                                                                                  
									sizeof(cl_int),                                                                                                                     
									(void *)&is_even);                                                                                                            
	if(status != CL_SUCCESS)                                                                                                                                            
	{                                                                                                                                                                   
			std::cout<<"Error: Setting kernel argument. (is_even)\n";                                                                                                    
			return 1;                                                                                                                                                   
	}

	// only needed in the local-memory version
	if(SIM_TYPE==2)
	{
		// step-width during integration
		status = clSetKernelArg(                                                                                                                                           
										kernel,                                                                                                                             
										7,                                                                                                                                  
										sizeof(cl_float) * localThreads[0] * 4,                                                                                                                     
										NULL);                                                                                                            
		if(status != CL_SUCCESS)                                                                                                                                            
		{                                                                                                                                                                   
				std::cout<<"Error: Setting kernel argument. (blockeddata)\n";                                                                                                    
				return 1;                                                                                                                                                   
		}
	}



	for(int i=0; i<nstep; i++)
	{
		/*** Set appropriate arguments to the kernel ***/	
		// pos11 input-output
		status = clSetKernelArg(                                                                                                                                           
										kernel,                                                                                                                             
										0,                                                                                                                                  
										sizeof(cl_mem),                                                                                                                     
										(void *)&pos1g_buf);                                                                                                            
		if(status != CL_SUCCESS)                                                                                                                                            
		{                                                                                                                                                                   
				std::cout<<"Error: Setting kernel argument. (pos1g_buf)\n";                                                                                                    
				return 1;                                                                                                                                                   
		}


		// pos2 input-output
		status = clSetKernelArg(                                                                                                                                           
										kernel,                                                                                                                             
										1,                                                                                                                                  
										sizeof(cl_mem),                                                                                                                     
										(void *)&pos2g_buf);                                                                                                            
		if(status != CL_SUCCESS)                                                                                                                                            
		{                                                                                                                                                                   
				std::cout<<"Error: Setting kernel argument. (f2g_buf)\n";                                                                                                    
				return 1;                                                                                                                                                   
		}

		// is-even trigger
		status = clSetKernelArg(                                                                                                                                           
										kernel,                                                                                                                             
										3,                                                                                                                                  
										sizeof(cl_int),                                                                                                                     
										(void *)&is_even);                                                                                                            
		if(status != CL_SUCCESS)                                                                                                                                            
		{                                                                                                                                                                   
				std::cout<<"Error: Setting kernel argument. (is_even)\n";                                                                                                    
				return 1;                                                                                                                                                   
		}

		/***!TODO!: COMPLETE clEnqueueNDRangeKernel argument list ***/
		status = clEnqueueNDRangeKernel(                                                                                                                                    
			oclm->commandQueue,                                                                                                                              
									kernel,
									1,
									NULL,                                                                                                                                      
									globalThreads,
									localThreads,
									0,                                                                                                                                         
									NULL,                                                                                                                                      
									&events[0]); //
		if(status != CL_SUCCESS)                                                                                                                                            
		{                                                                                                                                                                   
				std::cout<<                                                                                                                                                 
						"Error: Enqueueing kernel onto command queue. (clEnqueueNDRangeKernel):" << status << "\n";
		} 

		/* wait for the kernel call to finish execution */                                                                                                                  
		status = clWaitForEvents(1, &events[0]);                                                                                                                            
		if(status != CL_SUCCESS)                                                                                                                                            
		{                                                                                                                                                                   
				std::cout<<                                                                                                                                                 
						"Error: Waiting for kernel run to finish.(clWaitForEvents)\n";                                                                                                                               
		}

		// Timing stuff to meassure the time spent for kernel
		elapsed_time += oclm->printProfilingInfo(events);
		clReleaseEvent(events[0]);

		// set output/input for next round
		if(is_even==1)
			is_even=0;
		else
			is_even=1;
	}
	///////////////// End of iterative for

	/* Enqueue readBuffer*/                                                                                                                                             
	status = clEnqueueReadBuffer(                                                                                                                                       
							oclm->commandQueue,                                                                                                                               
							pos1g_buf,                                                                                                                              
							CL_TRUE,                                                                                                                                    
							0,                                                                                                                                          
							nparticle * 4 * sizeof(cl_float),                                                                                                         
							pos1g,                                                                                                                                     
							0,                                                                                                                                          
							NULL,                                                                                                                                       
							&events[1]);                                                                                                                                
                                                                                                                                                                            
	if(status != CL_SUCCESS)                                                                                                                                            
	{                                                                                                                                                                   
			std::cout <<                                                                                                                                                
					"Error: clEnqueueReadBuffer failed. (clEnqueueReadBuffer)\n";                                                                                                                          
	}  

    /* Wait for the read buffer to finish execution */                                                                                                                  
    status = clWaitForEvents(1, &events[1]);                                                                                                                            
    if(status != CL_SUCCESS)                                                                                                                                            
    {                                                                                                                                                                   
            std::cout<<                                                                                                                                                 
                    "Error: Waiting for read buffer call to finish. \n";                                                                                                
    }                                                                                                                                                                   
                                                                                                                                                                            
    clReleaseEvent(events[1]);
	printf("GPU Kernels needed %fl seconds\n", elapsed_time);

	return 0;
}

int cleanCL() {
	clReleaseMemObject(pos1g_buf);
	clReleaseMemObject(pos2g_buf);
	clReleaseMemObject(velg_buf);

	delete [] pos1g;
	delete [] pos2g;
	delete [] velg;

	return 0;
}

int main(int argc, char *argv[]) {
	SIM_TYPE = 0;

    cout << "starting!" << endl;
	for (int i = 1; i < argc; ++i) {

		std::string arg = std::string(argv[i]);

		// Search for type of kernel
		if ((arg == "-k") || (arg == "--kernel")) {
			if (i + 1 < argc) {
				SIM_TYPE = atoi(argv[i+1]);
				std::cout << "Setting kernel type to " << SIM_TYPE << ".\n";
			}
		}
	}

	dt = 0.0001f;
	eps = 0.0001f;

	cout << "Reading Input Values from File.\n";
	readExampleBodyValsFile();
	cout << "Read Solution from File.\n";
	readSolutionFromFile();	

	cout << "setting up OpenCL.\n";
	setupCL();

    cout << "running OpenCL.\n";
	runCL();

	cout << "Starting Evaluation.\n";
	evalAgainstFile();

	cleanCL();

	delete[] pos1;
	delete[] pos2;
	delete[] vel;

	return 0;
}
