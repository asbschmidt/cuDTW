#ifndef DTW_HPP
#define DTW_HPP

namespace FullDTW {

#include "./kernels/SHFL_FULLDTW_127.cuh"
#include "./kernels/SHFL_FULLDTW_255.cuh"
#include "./kernels/SHFL_FULLDTW_511.cuh"
#include "./kernels/SHFL_FULLDTW_1023.cuh"
#include "./kernels/SHFL_FULLDTW_2047.cuh"
#include "./kernels/SUB_WARP.cuh"
#include "./kernels/SUB_WARP_MULTI_QUERY.cuh"
#include "./kernels/SUB_WARP_FULLDTW_ECG.cuh"
#include "./kernels/DTW_ECG_2047.cuh"
#include "./kernels/ECG_1023_EARLY_EXIT_DTW.cuh"
#include "./kernels/DTW_ECG_2047_EARLY_EXIT.cuh"
#include "./kernels/ECG_511_EARLY_EXIT_DTW.cuh"
#include "./kernels/ECG_255_EARLY_EXIT_DTW.cuh"
#include "./kernels/ECG_127_EARLY_EXIT_DTW.cuh"

template <
    typename value_t,
    typename index_t> __host__
void dist (
    value_t * Subject,
    value_t * Dist,
    index_t num_entries,
    index_t num_features,
    index_t num_queries,
    bool subwarp,
    uint8_t query_type,
    bool lower_bound_stream,
    value_t * bsf,
    cudaStream_t stream=0) {

    if (query_type == 0) {
        if (num_features == 127) {
             if (subwarp) {
                 const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
                 const dim3 block(   32, 1, 1);
                 if (num_queries == 1) sub_warp_DTW<4><<<grid, block, 0, stream>>>
                    (Subject, Dist, num_entries, num_features);
                 else if (num_queries == 2) sub_warp_DTW_multi_query<4, 2><<<grid, block, 0, stream>>>
                    (Subject, Dist, num_entries, num_features);
                 else if (num_queries == 4) sub_warp_DTW_multi_query<4, 4><<<grid, block, 0, stream>>>
                    (Subject, Dist, num_entries, num_features);
                 else if (num_queries == 8) sub_warp_DTW_multi_query<4, 8><<<grid, block, 0, stream>>>
                    (Subject, Dist, num_entries, num_features);
                 else if (num_queries == 16) sub_warp_DTW_multi_query<4, 16><<<grid, block, 0, stream>>>
                    (Subject, Dist, num_entries, num_features);
                 return;
            } else {
                const dim3 grid (num_entries, 1, 1);
                const dim3 block(         32, 1, 1);
                shfl_FullDTW_127 <<<grid, block, 0, stream>>>
                    (Subject, Dist, num_entries, num_features);
                return;
            }
        }

        if (num_features == 255) {
        if (subwarp) {
             const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
             const dim3 block(   32, 1, 1);
             if (num_queries == 1) sub_warp_DTW<8><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
             else if (num_queries == 2) sub_warp_DTW_multi_query<8, 2><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
             else if (num_queries == 4) sub_warp_DTW_multi_query<8, 4><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
             else if (num_queries == 8) sub_warp_DTW_multi_query<8, 8><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
             else if (num_queries == 16) sub_warp_DTW_multi_query<8, 16><<<grid, block, 0, stream>>>
                  (Subject, Dist, num_entries, num_features);
             return;
        } else {
            const dim3 grid (num_entries, 1, 1);
            const dim3 block(         32, 1, 1);
            shfl_FullDTW_255 <<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
            return;
        }
        }

        if (num_features == 511) {
        if (subwarp) {
             const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
             const dim3 block(   32, 1, 1);
             if (num_queries == 1) sub_warp_DTW<16><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
            else if (num_queries == 2) sub_warp_DTW_multi_query<16, 2><<<grid, block, 0, stream>>>
                 (Subject, Dist, num_entries, num_features);
            else if (num_queries == 4) sub_warp_DTW_multi_query<16, 4><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
            else if (num_queries == 8) sub_warp_DTW_multi_query<16, 8><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
            else if (num_queries == 16) sub_warp_DTW_multi_query<16, 16><<<grid, block, 0, stream>>>
                 (Subject, Dist, num_entries, num_features);

            return;
        } else {
            const dim3 grid (num_entries, 1, 1);
            const dim3 block(         32, 1, 1);
            shfl_FullDTW_511 <<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
            return;
        }
        }

        if (num_features == 1023) {
        if (subwarp) {
             const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
             const dim3 block(   32, 1, 1);
             if (num_queries == 1) sub_warp_DTW<32><<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
             else if (num_queries == 2) sub_warp_DTW_multi_query<32, 2><<<grid, block, 0, stream>>>
                  (Subject, Dist, num_entries, num_features);
             else if (num_queries == 4) sub_warp_DTW_multi_query<32, 4><<<grid, block, 0, stream>>>
                  (Subject, Dist, num_entries, num_features);
             else if (num_queries == 8) sub_warp_DTW_multi_query<32, 8><<<grid, block, 0, stream>>>
                 (Subject, Dist, num_entries, num_features);
             else if (num_queries == 16) sub_warp_DTW_multi_query<32, 16><<<grid, block, 0, stream>>>
                   (Subject, Dist, num_entries, num_features);

             return;
        } else {
            const dim3 grid (num_entries, 1, 1);
            const dim3 block(         32, 1, 1);
            shfl_FullDTW_1023<<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
            return;
        }
        }

        if (num_features == 2047) {
            const dim3 grid (num_entries, 1, 1);
            const dim3 block(         32, 1, 1);
            shfl_FullDTW_2047<<<grid, block, 0, stream>>>
                (Subject, Dist, num_entries, num_features);
            return;
        }

        std::cout << "ERROR: length not supported, exiting." << std::endl;
        exit(1);
    }
    else if (query_type == 1) {  // ECG

        if (num_features == 1023) {
             if (subwarp) {
                  const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
                  const dim3 block(   32, 1, 1);
                  sub_warp_FullDTW_ECG<32><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features);
            }
            else {
                //value_t threshold = 0.4;
                const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
                const dim3 block(   32, 1, 1);
                //ECG_1023_early_exit_DTW<32><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, threshold);
                ECG_1023_early_exit_DTW<32><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, bsf);
            }
        }

        if (num_features == 511) {
             if (subwarp) {
                  const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
                  const dim3 block(   32, 1, 1);
                  sub_warp_FullDTW_ECG<16><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features);
            }
            else {
                //value_t threshold = 0.4;
                const dim3 grid (num_entries, 1, 1);
                const dim3 block(   32, 1, 1);
                //ECG_1023_early_exit_DTW<32><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, threshold);
                ECG_511_early_exit_DTW<<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, bsf);
            }

        }

        if (num_features == 255) {
             if (subwarp) {
                  const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
                  const dim3 block(   32, 1, 1);
                  sub_warp_FullDTW_ECG<8><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features);
            }
            else {
                //value_t threshold = 0.4;
                const dim3 grid (num_entries, 1, 1);
                const dim3 block(   32, 1, 1);
                //ECG_1023_early_exit_DTW<32><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, threshold);
                ECG_255_early_exit_DTW<<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, bsf);
            }

        }

        if (num_features == 127) {
             if (subwarp) {
                  const dim3 grid (num_entries/(1024/(num_features+1)), 1, 1);
                  const dim3 block(   32, 1, 1);
                  sub_warp_FullDTW_ECG<4><<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features);
            }
            else {
                const dim3 grid (num_entries, 1, 1);
                const dim3 block(   32, 1, 1);
                ECG_127_early_exit_DTW<<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, bsf);
            }

        }

        if (num_features == 2047) {
            if (subwarp) {
                  const dim3 grid (num_entries, 1, 1);
                  const dim3 block(   32, 1, 1);
                  DTW_ecg_2047<<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features);
             }
             else {
                 //value_t threshold = 1.0;
                 const dim3 grid (num_entries, 1, 1);
                 const dim3 block(   32, 1, 1);
                 DTW_ecg_2047_early_exit<<<grid, block, 0, stream>>>(Subject, Dist, num_entries, num_features, bsf);

             }
        }

    }
}

}



#endif
