__kernel void nbody (	__global float4* pos1 , __global float4* pos2, __global float4* vel, int is_even,
							int nparticle, float dt, float eps
							)
{

	int i = get_global_id(0);
	int n = get_global_size(0);


	__global float4* source;
	__global float4* destination;
	if(is_even) {
	    source = pos1;
	    destination = pos2;
	} else {
	    source = pos2;
	    destination = pos1;
	}

	float4 accel = (float4)(0.0, 0.0, 0.0, 0.0);
	for(int j = 0; j < n; j++) {
	    // the 4th diff is the mass, must be ignored later on!
	    float4 diff = source[j] - source[i];

	    float invr = 1.0 / sqrt( diff.x * diff.x + diff.y * diff.y * diff.z * diff.z + eps);
	    
	    float invr3 = invr * invr * invr;
	    
	    float f = source[j].w * invr3;
	    
	    accel += f * diff;

	}
	
	destination[i].xyz = source[i].xyz + dt * vel[i].xyz + 0.5 * dt * dt * accel.xyz;

	vel[i].xyz + dt * accel.xyz;

	return;
}