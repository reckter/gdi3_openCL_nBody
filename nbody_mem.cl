__kernel void nbody (	__global float4* pos1 , __global float4* pos2, __global float4* vel, int is_even,
							int nparticle, float dt, float eps, __local float4* blockdata
							)
{

	// Temp Pointers to make code clearer
	__global float4* source;
	__global float4* destination;

	if(is_even) {
	    source = pos1;
	    destination = pos2;
	} else {
	    source = pos2;
	    destination = pos1;
	}

	int i = get_global_id(0);
	int ti = get_local_id(0);

	int n = get_global_size(0);
	int nt = get_local_size(0); //  blockdata.length

	int nb = n/nt; // number of Work-groups


    float3 dt_f3 = (float3)(dt, dt, dt);
    float3 pointFive = (float3)(0.5f, 0.5f, 0.5f);

    float4 accel = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

	for(int j = 0; j < nb; j++) {
	    event_t event = async_work_group_copy(
	                    blockdata,
	                    source + j * nt,
	                    nt,
	                    0);

        barrier(CLK_LOCAL_MEM_FENCE);
        wait_group_events(1, &event);


        for(int k = 0; k < nt; k++) {
            // the 4th diff is the mass, must be ignored later on!
            float4 diff = blockdata[k] - source[i];

            float invr = 1.0 / sqrt( diff.x * diff.x + diff.y * diff.y + diff.z * diff.z + eps);

            float invr3 = invr * invr * invr;

            float f = blockdata[k].w * invr3;

            float4 f_f4 = (float4)(f, f, f, f);

            accel += f_f4 * diff;

        }

	}

	destination[i].xyz = source[i].xyz + dt_f3 * vel[i].xyz + pointFive * dt_f3 * dt_f3 * accel.xyz;

	vel[i].xyz += dt_f3 * accel.xyz;


	return;

}