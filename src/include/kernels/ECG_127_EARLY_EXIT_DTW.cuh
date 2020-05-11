#ifndef ECG_127_EARLY_EXIT_DTW
#define ECG_127_EARLY_EXIT_DTW
// Using 4 values per thread no shared memory
template <
    typename index_t,
    typename value_t> __global__
void ECG_127_early_exit_DTW(  // was void DTW_fast_128_shuffle_kernel_no_shared_memory
    value_t * Subject,
    value_t * Dist,
    index_t num_entries,
    index_t num_features,
    value_t * threshold) {

    const index_t blid = blockIdx.x;
    const index_t thid = threadIdx.x;
    const index_t lane = num_features+1;
    const index_t base = blid;
    const index_t WARP_SIZE = 32;
    const index_t l = thid;

    //extern __shared__ value_t Subject_cache[];

    value_t penalty_left = INFINITY;
    value_t penalty_diag = 0;  // INFINITY;
    value_t penalty_here0 = INFINITY; // 0;
    value_t penalty_here1 = INFINITY; // 0;
    value_t penalty_here2 = INFINITY; // 0;
    value_t penalty_here3 = INFINITY; // 0;
    value_t penalty_temp0;
    value_t penalty_temp1;


        if (thid == 0) {
            //if (iter > 0) penalty_left = Subject_cache[1]; else penalty_left = INFINITY;
            penalty_left = INFINITY;
            penalty_diag = INFINITY;
            penalty_here0 = INFINITY; // 0;
            penalty_here1 = INFINITY; // 0;
            penalty_here2 = INFINITY; // 0;
            penalty_here3 = INFINITY; // 0;

        }

        value_t subject_value0 = l == 0 ? 0 : Subject[base+4*l-1];
        value_t subject_value1 = Subject[base+4*l-0];
        value_t subject_value2 = Subject[base+4*l+1];
        value_t subject_value3 = Subject[base+4*l+2];


                #ifdef LOCAL_ZNORM_STREAM

                value_t X  = sum_warp(subject_value0
                           +          subject_value1
                           +          subject_value2
                           +          subject_value3);

                // mean value
                X /= num_features;

                value_t Y  = sum_warp(subject_value0  * subject_value0
                           +          subject_value1  * subject_value1
                           +          subject_value2  * subject_value2
                           +          subject_value3  * subject_value3);

                // inverse standard deviation
                Y = cuda_rsqrt(Y/num_features - X*X);

                subject_value0  = (subject_value0 -X)*Y;
                subject_value1  = (subject_value1 -X)*Y;
                subject_value2  = (subject_value2 -X)*Y;
                subject_value3  = (subject_value3 -X)*Y;

                #endif

                unsigned mask = 1;

                value_t S_first = subject_value1; // Subject[base];
                S_first = __shfl_sync(0xffffffff, S_first, 0);
                value_t S_second = subject_value2; // Subject[base+1];
                S_second = __shfl_sync(0xffffffff, S_second, 0);
                value_t S_third = subject_value3; // Subject[base+2];
                S_third = __shfl_sync(0xffffffff, S_third, 0);

                value_t S_last = subject_value3; // Subject[base+num_features-1];
                S_last = __shfl_sync(0xffffffff, S_last, 31);
                value_t S_second_last = subject_value2; // // Subject[base+num_features-2];
                S_second_last = __shfl_sync(0xffffffff, S_second_last, 31);
                value_t S_third_last = subject_value1; // Subject[base+num_features-3];
                S_third_last = __shfl_sync(0xffffffff, S_third_last, 31);


                const value_t Q_first = cQuery[0];
                const value_t Q_second = cQuery[1];
                const value_t Q_third = cQuery[2];
                const value_t Q_last = cQuery[num_features-1];
                const value_t Q_second_last = cQuery[num_features-2];
                const value_t Q_third_last = cQuery[num_features-3];

                value_t upper_left = (S_first-Q_first)*(S_first-Q_first);
                value_t right_upper_left = upper_left + (S_first-Q_second)*(S_first-Q_second);
                value_t lower_upper_left = upper_left + (S_second-Q_first)*(S_second-Q_first);
                upper_left += (S_second-Q_second)*(S_second-Q_second);
                value_t lower_middle_upper_left = (S_third-Q_second)*(S_third-Q_second) + min(lower_upper_left,upper_left);
                value_t right_middle_upper_left = (S_second-Q_third)*(S_second-Q_third) + min(right_upper_left,upper_left);
                upper_left = min(upper_left + (S_third-Q_third)*(S_third-Q_third), min(right_upper_left + (S_first-Q_third)*(S_first-Q_third), lower_upper_left + (S_third-Q_first)*(S_third-Q_first)));
                upper_left = min(upper_left, min(lower_middle_upper_left,right_middle_upper_left));

                value_t lower_right = (S_last-Q_last)*(S_last-Q_last);
                right_upper_left = lower_right + (S_last-Q_second_last)*(S_last-Q_second_last);
                lower_upper_left = lower_right + (S_second_last-Q_last)*(S_second_last-Q_last);
                lower_right += (S_second_last-Q_second_last)*(S_second_last-Q_second_last);
                lower_middle_upper_left = (S_third_last-Q_second_last)*(S_third_last-Q_second_last) + min(lower_upper_left,lower_right);
                right_middle_upper_left = (S_second_last-Q_third_last)*(S_second_last-Q_third_last) + min(right_upper_left,lower_right);
                lower_right = min(lower_right + (S_third_last-Q_third_last)*(S_third_last-Q_third_last), min(right_upper_left + (S_last-Q_third_last)*(S_last-Q_third_last), lower_upper_left + (S_third_last-Q_last)*(S_third_last-Q_last)));
                lower_right = min(lower_right, min(lower_middle_upper_left,right_middle_upper_left));

                if (upper_left + lower_right >= *threshold) {
                    Dist[blid] = 100000;
                } else {


        index_t counter = 1;
        value_t query_value = INFINITY;
        value_t new_query_value = cQuery[thid];
        if (thid == 0) query_value = new_query_value;
        if (thid == 0) penalty_here1 = 0; // (query_value - subject_value0)*(query_value - subject_value0);
        //penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here1, 1, 32);
        new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
        //const index_t j = l;
        //if (blid == 0 && thid == 31 && iter == 0) Dist[0] = penalty_here0;
        //if (blid == 0 && thid == 31 && iter == 0) Dist[1] = penalty_here1;
        //if (blid == 0 && iter == 0) Dist[2*thid] = penalty_here0;
        //if (blid == 0 && iter == 0) Dist[2*thid+1] = penalty_here1;

        penalty_temp0 = penalty_here0;
        penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
        //if (i==2) penalty_temp1 = INFINITY; else penalty_temp1 = penalty_here1;  // -> move before main loop!!!
        penalty_temp1 = INFINITY;
        penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
        penalty_temp0 = penalty_here2;
        penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
        penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));

        query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
        if (thid == 0) query_value = new_query_value;
        new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
        counter++;

        penalty_diag = penalty_left;
        penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here3, 1, 32);

        //if (iter && thid == 0) penalty_left = Subject_cache[2+1];
        if (thid == 0) penalty_left = INFINITY;


        for (index_t k = 3; k < lane+WARP_SIZE-1; k++) {
            const index_t i = k-l;
            //outside = k <= l || i >= lane;

            //const value_t residue = outside ? INFINITY : cQuery[i-1]-subject_value;
            //const value_t residue = outside ? INFINITY : query_value-subject_value;
            //if (thid == 0 && iter == 0 && k == 2) penalty_temp = INFINITY; else
            penalty_temp0 = penalty_here0;
            penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
            //if (i==2) penalty_temp1 = INFINITY; else penalty_temp1 = penalty_here1;  // -> move before main loop!!!
            penalty_temp1 = penalty_here1;
            penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
            penalty_temp0 = penalty_here2;
            penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
            penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));

            //if (counter%32 == 0) new_query_value = cQuery[i+2*thid-1];

            if (counter%32 == 0) {
                penalty_temp0 = min(min(min(min(penalty_left, penalty_here0), penalty_here1), penalty_here2), penalty_here3);
                mask = __any_sync(0xFFFFFFFF, penalty_temp0 < *threshold - lower_right);
                if (!mask) {Dist[blid] = 100*(counter/32); break; }
                new_query_value = cQuery[i+2*thid-1];
            }

            query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
            if (thid == 0) query_value = new_query_value;
            new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
            //if (thid == 0) if (!outside) Dist[counter] = query_value; else Dist[counter] = 0;
            counter++;

            // save the right column
            //if (!outside && thid == 31 && iter < lane/WARP_SIZE-1) Subject_cache[i] = penalty_here; // TO DO: replace this by shhffles
            //if (iter < lane/WARP_SIZE-1 && thid == 31 && k>l) Subject_cache[i] = penalty_here3;

            // shuffle the penalty
            penalty_diag = penalty_left;
            penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here3, 1, 32);
            //if (thid == 0 && !outside) penalty_left = Subject_cache[i+1]; // TO DO: replace by shuffles
        //    if (iter > 0 && thid == 0 && k>l) penalty_left = Subject_cache[i+1];
            //if (iter && thid == 0) penalty_left = Subject_cache[i+1];
            //if (!iter && thid == 0) penalty_left = INFINITY;
            if (thid == 0) penalty_left = INFINITY;
            //if (thid == 0 && k>l) if (iter > 0) penalty_left = Subject_cache[i+1]; else penalty_left = INFINITY;
            //if (thid == 0 && k>l && iter > 0) penalty_left = Subject_cache[i+1];
        }
        penalty_temp0 = penalty_here0;
        penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
        penalty_temp1 = penalty_here1;
        penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
        penalty_temp0 = penalty_here2;
        penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
        penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));

        //if (blid == 15) {
        //    Dist[128*lane + 4*l] = penalty_here0;
        //    Dist[128*lane + 4*l+1] = penalty_here1;
        //    Dist[128*lane + 4*l+2] = penalty_here2;
        //    Dist[128*lane + 4*l+3] = penalty_here3;
        //}

        //if (thid == 31 && iter < lane/WARP_SIZE-1) Subject_cache[lane+(iter+1)*WARP_SIZE-1-l] = penalty_here3;
        //iter++;
    if(thid == blockDim.x-1 && mask) {
        Dist[blid] = penalty_here3;
        if (penalty_here3 < *threshold)
            update_bsf_system(threshold, penalty_here3);
    }

    //if(thid == blockDim.x-1)  Dist[blid] = penalty_here3;
}
}

#endif
