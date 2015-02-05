__kernel void nbody (	__global float* pos1 , __global float* pos2, __global float* vel, int is_even,
							int nparticle, float dt, float eps
							)
{

	// Temp Pointers to make code clearer
	__global float* pos1_t;
	__global float* pos2_t;

	int i = get_global_id(0);
	int n = get_global_size(0);

	return;
}