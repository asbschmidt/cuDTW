#ifndef ED_HPP
#define ED_HPP

namespace ED {

template < 
    typename value_t,
    typename index_t> __global__
void dist_kernel(
 // value_t * query stored in constant memory cquery
    value_t * data, 
    value_t * dist, 
    index_t   num_entries,
    index_t   num_features) {

    for (index_t entry = blockIdx.x; entry < num_entries; entry += gridDim.x) {
        
        value_t accum = 0.0;
        for (index_t index = threadIdx.x; 
                     index < num_features; index += blockDim.x)
            accum += (cQuery[index]-data[entry*num_features+index])
                  *  (cQuery[index]-data[entry*num_features+index]);

        for (index_t offset = warpSize >> 1; offset > 0; offset >>= 1)
            accum += __shfl_down_sync(0xFFFFFFFF, accum, offset, warpSize);

        if (threadIdx.x % warpSize == 0)
            atomicAdd(dist+entry, accum);  
    }
}


}
#endif
