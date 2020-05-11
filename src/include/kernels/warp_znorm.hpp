#ifndef WARP_ZNORM_HPP
#define WARP_ZNORM_HPP

#include <cmath>

template <
    typename value_t,
    int warpsize =32> __device__
value_t sum_warp(value_t value) {
    
    value_t accum = value;

    for (int delta = warpsize>>1; delta > 0; delta >>= 1)
        accum += __shfl_xor_sync(0xFFFFFFFF, accum, delta, warpsize);

    return accum;
}

template <
    typename value_t,
    int warpsize =32> __device__
value_t sum_square_warp(value_t value) {
    
    value_t accum = value*value;

    for (int delta = warpsize>>1; delta > 0; delta >>= 1)
        accum += __shfl_xor_sync(0xFFFFFFFF, accum, delta, warpsize);

    return accum;
}

__device__
float cuda_rsqrt(float x) {
    return x > 0 ? rsqrtf(x) : 1;
}

__device__
double cuda_rsqrt(double x) {
    return x > 0 ? rsqrt(x)  : 1;
}

template <
    typename value_t,
    typename index_t>
void cpu_znorm(value_t * series, index_t length) {

    value_t X = 0, Y = 0;

    for (index_t index = 0; index < length; index++) {
        const value_t value = series[index];
        X += value;
        Y += value*value;
    }

    X /= length;
    Y  = Y/length-X*X;
    Y  = Y > 0 ? 1.0/std::sqrt(Y) : 1;

    std::cout << "STATUS: query stats mu = " << X << " sigma = " << 1.0/Y << std::endl;

    for (index_t index = 0; index < length; index++) {
        series[index]=(series[index]-X)*Y;
    }
}

#endif