__kernel void nbody (	__global float4* pos1 , __global float4* pos2, __global float4* vel, int is_even,
							int nparticle, float dt, float eps, __local float4* blockdata
							)
{

	// Temp Pointers to make code clearer
	__global float4* pos1_t;
	__global float4* pos2_t;

	int i = get_global_id(0);
	int ti = get_local_id(0);

	int n = get_global_size(0);
	int nt = get_local_size(0);

	int nb = n/nt; // number of Work-groups

	return;

}