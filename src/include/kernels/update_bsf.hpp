#ifndef UPDATE_BSF_HPP
#define UPDATE_BSF_HPP

__device__ __forceinline__
float  update_bsf_system(float * threshold, float value) {
    
    // exploiting that positive ints and floats map monotonously

    if (value < 0)
	return INFINITY;

    unsigned int * ptr =   (unsigned int *) threshold;
    unsigned int   val = *((unsigned int *) (&value));    

    const auto old = atomicMin_system(ptr, val);
    const auto result = *((float*) (&old));

    return result;
}

#endif
