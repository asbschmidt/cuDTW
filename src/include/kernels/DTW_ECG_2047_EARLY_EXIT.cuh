#ifndef DTW_ECG_2047_EARLY_EXIT
#define DTW_ECG_2047_EARLY_EXIT

#include "update_bsf.hpp"
#include "warp_znorm.hpp"
#include <stdio.h>

// 64 values per thread no shared memory
template <
    typename index_t,
    typename value_t> __global__
void DTW_ecg_2047_early_exit(
    value_t * Subject,
    value_t * Dist,
    index_t num_entries,
    index_t num_features,
    value_t * threshold) {

    const index_t blid = blockIdx.x;
    const index_t thid = threadIdx.x;
    const index_t lane = num_features+1;
    //const index_t base = blid*num_features;
    const index_t base = blid;

    const index_t l = thid;
    //extern __shared__ value_t Subject_cache[];

    value_t penalty_left = INFINITY;
    value_t penalty_diag = 0;  // INFINITY;
    value_t penalty_here0 = INFINITY; // 0;
    value_t penalty_here1 = INFINITY; // 0;
    value_t penalty_here2 = INFINITY; // 0;
    value_t penalty_here3 = INFINITY; // 0;
    value_t penalty_here4 = INFINITY; // 0;
    value_t penalty_here5 = INFINITY; // 0;
    value_t penalty_here6 = INFINITY; // 0;
    value_t penalty_here7 = INFINITY; // 0;
    value_t penalty_here8 = INFINITY; // 0;
    value_t penalty_here9 = INFINITY; // 0;
    value_t penalty_here10 = INFINITY; // 0;
    value_t penalty_here11 = INFINITY; // 0;
    value_t penalty_here12 = INFINITY; // 0;
    value_t penalty_here13 = INFINITY; // 0;
    value_t penalty_here14 = INFINITY; // 0;
    value_t penalty_here15 = INFINITY; // 0;
    value_t penalty_here16 = INFINITY; // 0;
    value_t penalty_here17 = INFINITY; // 0;
    value_t penalty_here18 = INFINITY; // 0;
    value_t penalty_here19 = INFINITY; // 0;
    value_t penalty_here20 = INFINITY; // 0;
    value_t penalty_here21 = INFINITY; // 0;
    value_t penalty_here22 = INFINITY; // 0;
    value_t penalty_here23 = INFINITY; // 0;
    value_t penalty_here24 = INFINITY; // 0;
    value_t penalty_here25 = INFINITY; // 0;
    value_t penalty_here26 = INFINITY; // 0;
    value_t penalty_here27 = INFINITY; // 0;
    value_t penalty_here28 = INFINITY; // 0;
    value_t penalty_here29 = INFINITY; // 0;
    value_t penalty_here30 = INFINITY; // 0;
    value_t penalty_here31 = INFINITY; // 0;
    value_t penalty_here32 = INFINITY; // 0;
    value_t penalty_here33 = INFINITY; // 0;
    value_t penalty_here34 = INFINITY; // 0;
    value_t penalty_here35 = INFINITY; // 0;
    value_t penalty_here36 = INFINITY; // 0;
    value_t penalty_here37 = INFINITY; // 0;
    value_t penalty_here38 = INFINITY; // 0;
    value_t penalty_here39 = INFINITY; // 0;
    value_t penalty_here40 = INFINITY; // 0;
    value_t penalty_here41 = INFINITY; // 0;
    value_t penalty_here42 = INFINITY; // 0;
    value_t penalty_here43 = INFINITY; // 0;
    value_t penalty_here44 = INFINITY; // 0;
    value_t penalty_here45 = INFINITY; // 0;
    value_t penalty_here46 = INFINITY; // 0;
    value_t penalty_here47 = INFINITY; // 0;
    value_t penalty_here48 = INFINITY; // 0;
    value_t penalty_here49 = INFINITY; // 0;
    value_t penalty_here50 = INFINITY; // 0;
    value_t penalty_here51 = INFINITY; // 0;
    value_t penalty_here52 = INFINITY; // 0;
    value_t penalty_here53 = INFINITY; // 0;
    value_t penalty_here54 = INFINITY; // 0;
    value_t penalty_here55 = INFINITY; // 0;
    value_t penalty_here56 = INFINITY; // 0;
    value_t penalty_here57 = INFINITY; // 0;
    value_t penalty_here58 = INFINITY; // 0;
    value_t penalty_here59 = INFINITY; // 0;
    value_t penalty_here60 = INFINITY; // 0;
    value_t penalty_here61 = INFINITY; // 0;
    value_t penalty_here62 = INFINITY; // 0;
    value_t penalty_here63 = INFINITY; // 0;
    value_t penalty_temp0;
    value_t penalty_temp1;

    //unsigned mask = 1;
    //const value_t S_first = Subject[base];
    //const value_t S_second = Subject[base+1];
    //const value_t S_third = Subject[base+2];
    //const value_t S_last = Subject[base+num_features-1];
    //const value_t S_second_last = Subject[base+num_features-2];
    //const value_t S_third_last = Subject[base+num_features-3];

    //const value_t Q_first = cQuery[0];
    //const value_t Q_second = cQuery[1];
    //const value_t Q_third = cQuery[2];
    //const value_t Q_last = cQuery[num_features-1];
    //const value_t Q_second_last = cQuery[num_features-2];
    //const value_t Q_third_last = cQuery[num_features-3];

    //value_t upper_left = (S_first-Q_first)*(S_first-Q_first);
    //value_t right_upper_left = upper_left + (S_first-Q_second)*(S_first-Q_second);
    //value_t lower_upper_left = upper_left + (S_second-Q_first)*(S_second-Q_first);
    //upper_left += (S_second-Q_second)*(S_second-Q_second);
    //value_t lower_middle_upper_left = (S_third-Q_second)*(S_third-Q_second) + min(lower_upper_left,upper_left);
    //value_t right_middle_upper_left = (S_second-Q_third)*(S_second-Q_third) + min(right_upper_left,upper_left);
    //upper_left = min(upper_left + (S_third-Q_third)*(S_third-Q_third), min(right_upper_left + (S_first-Q_third)*(S_first-Q_third), lower_upper_left + (S_third-Q_first)*(S_third-Q_first)));
    //upper_left = min(upper_left, min(lower_middle_upper_left,right_middle_upper_left));

    //value_t lower_right = (S_last-Q_last)*(S_last-Q_last);
    //right_upper_left = lower_right + (S_last-Q_second_last)*(S_last-Q_second_last);
    //lower_upper_left = lower_right + (S_second_last-Q_last)*(S_second_last-Q_last);
    //lower_right += (S_second_last-Q_second_last)*(S_second_last-Q_second_last);
    //lower_middle_upper_left = (S_third_last-Q_second_last)*(S_third_last-Q_second_last) + min(lower_upper_left,lower_right);
    //right_middle_upper_left = (S_second_last-Q_third_last)*(S_second_last-Q_third_last) + min(right_upper_left,lower_right);
    //lower_right = min(lower_right + (S_third_last-Q_third_last)*(S_third_last-Q_third_last), min(right_upper_left + (S_last-Q_third_last)*(S_last-Q_third_last), lower_upper_left + (S_third_last-Q_last)*(S_third_last-Q_last)));
    //lower_right = min(lower_right, min(lower_middle_upper_left,right_middle_upper_left));

    //upper_left += min((S_first-Q_second)*(S_first-Q_second)+(S_first-Q_third)*(S_first-Q_third), min((S_second-Q_first)*(S_second-Q_first)+(S_third-Q_first)*(S_third-Q_first),(S_second-Q_second)*(S_second-Q_second)));
    //lower_right += min((S_last-Q_second_last)*(S_last-Q_second_last)+(S_last-Q_third_last)*(S_last-Q_third_last), min((S_second_last-Q_last)*(S_second_last-Q_last)+(S_third_last-Q_last)*(S_third_last-Q_last),(S_second_last-Q_second_last)*(S_second_last-Q_second_last)));
    //value_t trivial_lower_bound = upper_left + lower_right;
    //value_t trivial_lower_bound = (S_first-Q_first)*(S_first-Q_first) + (S_last-Q_last)*(S_last-Q_last);

    //if (upper_left + lower_right >= *threshold) {
    //    Dist[blid] = 100000;
    //} else {



    if (thid == 0) {
            penalty_left = INFINITY;
            penalty_diag = INFINITY;
            penalty_here0 = INFINITY; // 0;
            penalty_here1 = INFINITY; // 0;
            penalty_here2 = INFINITY; // 0;
            penalty_here3 = INFINITY; // 0;
            penalty_here4 = INFINITY; // 0;
            penalty_here5 = INFINITY; // 0;
            penalty_here6 = INFINITY; // 0;
            penalty_here7 = INFINITY; // 0;
            penalty_here8 = INFINITY; // 0;
            penalty_here9 = INFINITY; // 0;
            penalty_here10 = INFINITY; // 0;
            penalty_here11 = INFINITY; // 0;
            penalty_here12 = INFINITY; // 0;
            penalty_here12 = INFINITY; // 0;
            penalty_here13 = INFINITY; // 0;
            penalty_here14 = INFINITY; // 0;
            penalty_here15 = INFINITY; // 0;
            penalty_here16 = INFINITY; // 0;
            penalty_here17 = INFINITY; // 0;
            penalty_here18 = INFINITY; // 0;
            penalty_here19 = INFINITY; // 0;
            penalty_here20 = INFINITY; // 0;
            penalty_here21 = INFINITY; // 0;
            penalty_here22 = INFINITY; // 0;
            penalty_here23 = INFINITY; // 0;
            penalty_here24 = INFINITY; // 0;
            penalty_here25 = INFINITY; // 0;
            penalty_here26 = INFINITY; // 0;
            penalty_here27 = INFINITY; // 0;
            penalty_here28 = INFINITY; // 0;
            penalty_here29 = INFINITY; // 0;
            penalty_here30 = INFINITY; // 0;
            penalty_here31 = INFINITY; // 0;
    }

        //const value_t subject_value = Subject[base+l-1];
        //const value_t subject_value0 = Subject[base+64*l-1];
        value_t subject_value0 = l == 0 ? 0 : Subject[base+64*l-1];
        value_t subject_value1 = Subject[base+64*l-0];
        value_t subject_value2 = Subject[base+64*l+1];
        value_t subject_value3 = Subject[base+64*l+2];
        value_t subject_value4 = Subject[base+64*l+3];
        value_t subject_value5 = Subject[base+64*l+4];
        value_t subject_value6 = Subject[base+64*l+5];
        value_t subject_value7 = Subject[base+64*l+6];
        value_t subject_value8 = Subject[base+64*l+7];
        value_t subject_value9 = Subject[base+64*l+8];
        value_t subject_value10 = Subject[base+64*l+9];
        value_t subject_value11 = Subject[base+64*l+10];
        value_t subject_value12 = Subject[base+64*l+11];
        value_t subject_value13 = Subject[base+64*l+12];
        value_t subject_value14 = Subject[base+64*l+13];
        value_t subject_value15 = Subject[base+64*l+14];

        value_t subject_value16 = Subject[base+64*l+15];
        value_t subject_value17 = Subject[base+64*l+16];
        value_t subject_value18 = Subject[base+64*l+17];
        value_t subject_value19 = Subject[base+64*l+18];
        value_t subject_value20 = Subject[base+64*l+19];
        value_t subject_value21 = Subject[base+64*l+20];
        value_t subject_value22 = Subject[base+64*l+21];
        value_t subject_value23 = Subject[base+64*l+22];
        value_t subject_value24 = Subject[base+64*l+23];
        value_t subject_value25 = Subject[base+64*l+24];
        value_t subject_value26 = Subject[base+64*l+25];
        value_t subject_value27 = Subject[base+64*l+26];
        value_t subject_value28 = Subject[base+64*l+27];
        value_t subject_value29 = Subject[base+64*l+28];
        value_t subject_value30 = Subject[base+64*l+29];
        value_t subject_value31 = Subject[base+64*l+30];

        value_t subject_value32 = Subject[base+64*l+31];
        value_t subject_value33 = Subject[base+64*l+32];
        value_t subject_value34 = Subject[base+64*l+33];
        value_t subject_value35 = Subject[base+64*l+34];
        value_t subject_value36 = Subject[base+64*l+35];
        value_t subject_value37 = Subject[base+64*l+36];
        value_t subject_value38 = Subject[base+64*l+37];
        value_t subject_value39 = Subject[base+64*l+38];
        value_t subject_value40 = Subject[base+64*l+39];
        value_t subject_value41 = Subject[base+64*l+40];
        value_t subject_value42 = Subject[base+64*l+41];
        value_t subject_value43 = Subject[base+64*l+42];
        value_t subject_value44 = Subject[base+64*l+43];
        value_t subject_value45 = Subject[base+64*l+44];
        value_t subject_value46 = Subject[base+64*l+45];
        value_t subject_value47 = Subject[base+64*l+46];

        value_t subject_value48 = Subject[base+64*l+47];
        value_t subject_value49 = Subject[base+64*l+48];
        value_t subject_value50 = Subject[base+64*l+49];
        value_t subject_value51 = Subject[base+64*l+50];
        value_t subject_value52 = Subject[base+64*l+51];
        value_t subject_value53 = Subject[base+64*l+52];
        value_t subject_value54 = Subject[base+64*l+53];
        value_t subject_value55 = Subject[base+64*l+54];
        value_t subject_value56 = Subject[base+64*l+55];
        value_t subject_value57 = Subject[base+64*l+56];
        value_t subject_value58 = Subject[base+64*l+57];
        value_t subject_value59 = Subject[base+64*l+58];
        value_t subject_value60 = Subject[base+64*l+59];
        value_t subject_value61 = Subject[base+64*l+60];
        value_t subject_value62 = Subject[base+64*l+61];
        value_t subject_value63 = Subject[base+64*l+62];
        //if (blid == 0) Dist[2*l] = subject_value0;
        //if (blid == 0) Dist[2*l+1] = subject_value1;


            #ifdef LOCAL_ZNORM_STREAM

            value_t X  = sum_warp(subject_value0
                       +          subject_value1
                       +          subject_value2
                       +          subject_value3
                       +          subject_value4
                       +          subject_value5
                       +          subject_value6
                       +          subject_value7
                       +          subject_value8
                       +          subject_value9
                       +          subject_value10
                       +          subject_value11
                       +          subject_value12
                       +          subject_value13
                       +          subject_value14
                       +          subject_value15
                       +          subject_value16
                       +          subject_value17
                       +          subject_value18
                       +          subject_value19
                       +          subject_value20
                       +          subject_value21
                       +          subject_value22
                       +          subject_value23
                       +          subject_value24
                       +          subject_value25
                       +          subject_value26
                       +          subject_value27
                       +          subject_value28
                       +          subject_value29
                       +          subject_value30
                       +          subject_value31
                       +                     subject_value32
                                  +          subject_value33
                                  +          subject_value34
                                  +          subject_value35
                                  +          subject_value36
                                  +          subject_value37
                                  +          subject_value38
                                  +          subject_value39
                                  +          subject_value40
                                  +          subject_value41
                                  +          subject_value42
                                  +          subject_value43
                                  +          subject_value44
                                  +          subject_value45
                                  +          subject_value46
                                  +          subject_value47
                                  +          subject_value48
                                  +          subject_value49
                                  +          subject_value50
                                  +          subject_value51
                                  +          subject_value52
                                  +          subject_value53
                                  +          subject_value54
                                  +          subject_value55
                                  +          subject_value56
                                  +          subject_value57
                                  +          subject_value58
                                  +          subject_value59
                                  +          subject_value60
                                  +          subject_value61
                                  +          subject_value62
                                  +          subject_value63);

            // mean value
            X /= num_features;

            value_t Y  = sum_warp(subject_value0  * subject_value0
                       +          subject_value1  * subject_value1
                       +          subject_value2  * subject_value2
                       +          subject_value3  * subject_value3
                       +          subject_value4  * subject_value4
                       +          subject_value5  * subject_value5
                       +          subject_value6  * subject_value6
                       +          subject_value7  * subject_value7
                       +          subject_value8  * subject_value8
                       +          subject_value9  * subject_value9
                       +          subject_value10 * subject_value10
                       +          subject_value11 * subject_value11
                       +          subject_value12 * subject_value12
                       +          subject_value13 * subject_value13
                       +          subject_value14 * subject_value14
                       +          subject_value15 * subject_value15
                       +          subject_value16 * subject_value16
                       +          subject_value17 * subject_value17
                       +          subject_value18 * subject_value18
                       +          subject_value19 * subject_value19
                       +          subject_value20 * subject_value20
                       +          subject_value21 * subject_value21
                       +          subject_value22 * subject_value22
                       +          subject_value23 * subject_value23
                       +          subject_value24 * subject_value24
                       +          subject_value25 * subject_value25
                       +          subject_value26 * subject_value26
                       +          subject_value27 * subject_value27
                       +          subject_value28 * subject_value28
                       +          subject_value29 * subject_value29
                       +          subject_value30 * subject_value30
                       +          subject_value31 * subject_value31
                       +                     subject_value32  * subject_value32
                                  +          subject_value33  * subject_value33
                                  +          subject_value34  * subject_value34
                                  +          subject_value35  * subject_value35
                                  +          subject_value36  * subject_value36
                                  +          subject_value37  * subject_value37
                                  +          subject_value38  * subject_value38
                                  +          subject_value39  * subject_value39
                                  +          subject_value40  * subject_value40
                                  +          subject_value41  * subject_value41
                                  +          subject_value42 * subject_value42
                                  +          subject_value43 * subject_value43
                                  +          subject_value44 * subject_value44
                                  +          subject_value45 * subject_value45
                                  +          subject_value46 * subject_value46
                                  +          subject_value47 * subject_value47
                                  +          subject_value48 * subject_value48
                                  +          subject_value49 * subject_value49
                                  +          subject_value50 * subject_value50
                                  +          subject_value51 * subject_value51
                                  +          subject_value52 * subject_value52
                                  +          subject_value53 * subject_value53
                                  +          subject_value54 * subject_value54
                                  +          subject_value55 * subject_value55
                                  +          subject_value56 * subject_value56
                                  +          subject_value57 * subject_value57
                                  +          subject_value58 * subject_value58
                                  +          subject_value59 * subject_value59
                                  +          subject_value60 * subject_value60
                                  +          subject_value61 * subject_value61
                                  +          subject_value62 * subject_value62
                                  +          subject_value63 * subject_value63 );

            // inverse standard deviation
            Y = cuda_rsqrt(Y/num_features - X*X);

            subject_value0  = (subject_value0 -X)*Y;
            subject_value1  = (subject_value1 -X)*Y;
            subject_value2  = (subject_value2 -X)*Y;
            subject_value3  = (subject_value3 -X)*Y;
            subject_value4  = (subject_value4 -X)*Y;
            subject_value5  = (subject_value5 -X)*Y;
            subject_value6  = (subject_value6 -X)*Y;
            subject_value7  = (subject_value7 -X)*Y;
            subject_value8  = (subject_value8 -X)*Y;
            subject_value9  = (subject_value9 -X)*Y;
            subject_value10 = (subject_value10-X)*Y;
            subject_value11 = (subject_value11-X)*Y;
            subject_value12 = (subject_value12-X)*Y;
            subject_value13 = (subject_value13-X)*Y;
            subject_value14 = (subject_value14-X)*Y;
            subject_value15 = (subject_value15-X)*Y;
            subject_value16 = (subject_value16-X)*Y;
            subject_value17 = (subject_value17-X)*Y;
            subject_value18 = (subject_value18-X)*Y;
            subject_value19 = (subject_value19-X)*Y;
            subject_value20 = (subject_value20-X)*Y;
            subject_value21 = (subject_value21-X)*Y;
            subject_value22 = (subject_value22-X)*Y;
            subject_value23 = (subject_value23-X)*Y;
            subject_value24 = (subject_value24-X)*Y;
            subject_value25 = (subject_value25-X)*Y;
            subject_value26 = (subject_value26-X)*Y;
            subject_value27 = (subject_value27-X)*Y;
            subject_value28 = (subject_value28-X)*Y;
            subject_value29 = (subject_value29-X)*Y;
            subject_value30 = (subject_value30-X)*Y;
            subject_value31 = (subject_value31-X)*Y;
            subject_value32  = (subject_value32 -X)*Y;
            subject_value33  = (subject_value33 -X)*Y;
            subject_value34  = (subject_value34 -X)*Y;
            subject_value35  = (subject_value35 -X)*Y;
            subject_value36  = (subject_value36 -X)*Y;
            subject_value37  = (subject_value37 -X)*Y;
            subject_value38  = (subject_value38 -X)*Y;
            subject_value39  = (subject_value39 -X)*Y;
            subject_value40  = (subject_value40 -X)*Y;
            subject_value41  = (subject_value41 -X)*Y;
            subject_value42 = (subject_value42-X)*Y;
            subject_value43 = (subject_value43-X)*Y;
            subject_value44 = (subject_value44-X)*Y;
            subject_value45 = (subject_value45-X)*Y;
            subject_value46 = (subject_value46-X)*Y;
            subject_value47 = (subject_value47-X)*Y;
            subject_value48 = (subject_value48-X)*Y;
            subject_value49 = (subject_value49-X)*Y;
            subject_value50 = (subject_value50-X)*Y;
            subject_value51 = (subject_value51-X)*Y;
            subject_value52 = (subject_value52-X)*Y;
            subject_value53 = (subject_value53-X)*Y;
            subject_value54 = (subject_value54-X)*Y;
            subject_value55 = (subject_value55-X)*Y;
            subject_value56 = (subject_value56-X)*Y;
            subject_value57 = (subject_value57-X)*Y;
            subject_value58 = (subject_value58-X)*Y;
            subject_value59 = (subject_value59-X)*Y;
            subject_value60 = (subject_value60-X)*Y;
            subject_value61 = (subject_value61-X)*Y;
            subject_value62 = (subject_value62-X)*Y;
            subject_value63 = (subject_value63-X)*Y;

            #endif

            unsigned mask = 1;

            value_t S_first = subject_value1; // Subject[base];
            S_first = __shfl_sync(0xffffffff, S_first, 0);
            value_t S_second = subject_value2; // Subject[base+1];
            S_second = __shfl_sync(0xffffffff, S_second, 0);
            value_t S_third = subject_value3; // Subject[base+2];
            S_third = __shfl_sync(0xffffffff, S_third, 0);

            value_t S_last = subject_value63; // Subject[base+num_features-1];
            S_last = __shfl_sync(0xffffffff, S_last, 31);
            value_t S_second_last = subject_value62; // // Subject[base+num_features-2];
            S_second_last = __shfl_sync(0xffffffff, S_second_last, 31);
            value_t S_third_last = subject_value61; // Subject[base+num_features-3];
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
    new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
        //const index_t j = l;
        //if (blid == 0 && thid == 31 && iter == 0) Dist[0] = penalty_here0;
        //if (blid == 0 && thid == 31 && iter == 0) Dist[1] = penalty_here1;
        //if (blid == 0 && iter == 0) Dist[2*thid] = penalty_here0;
        //if (blid == 0 && iter == 0) Dist[2*thid+1] = penalty_here1;

    penalty_temp0 = penalty_here0;
        penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
        penalty_temp1 = INFINITY;
        penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
        penalty_temp0 = penalty_here2;
        penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
        penalty_temp1 = penalty_here3;
        penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));
        penalty_temp0 = penalty_here4;
        penalty_here4 = (query_value-subject_value4) * (query_value-subject_value4) + min(penalty_here3, min(penalty_here4, penalty_temp1));
        penalty_temp1 = penalty_here5;
        penalty_here5 = (query_value-subject_value5) * (query_value-subject_value5) + min(penalty_here4, min(penalty_here5, penalty_temp0));
        penalty_temp0 = penalty_here6;
        penalty_here6 = (query_value-subject_value6) * (query_value-subject_value6) + min(penalty_here5, min(penalty_here6, penalty_temp1));
        penalty_temp1 = penalty_here7;
        penalty_here7 = (query_value-subject_value7) * (query_value-subject_value7) + min(penalty_here6, min(penalty_here7, penalty_temp0));
        penalty_temp0 = penalty_here8;
        penalty_here8 = (query_value-subject_value8) * (query_value-subject_value8) + min(penalty_here7, min(penalty_here8, penalty_temp1));
        penalty_temp1 = penalty_here9;
        penalty_here9 = (query_value-subject_value9) * (query_value-subject_value9) + min(penalty_here8, min(penalty_here9, penalty_temp0));
        penalty_temp0 = penalty_here10;
        penalty_here10 = (query_value-subject_value10) * (query_value-subject_value10) + min(penalty_here9, min(penalty_here10, penalty_temp1));
        penalty_temp1 = penalty_here11;
        penalty_here11 = (query_value-subject_value11) * (query_value-subject_value11) + min(penalty_here10, min(penalty_here11, penalty_temp0));
        penalty_temp0 = penalty_here12;
        penalty_here12 = (query_value-subject_value12) * (query_value-subject_value12) + min(penalty_here11, min(penalty_here12, penalty_temp1));
        penalty_temp1 = penalty_here13;
        penalty_here13 = (query_value-subject_value13) * (query_value-subject_value13) + min(penalty_here12, min(penalty_here13, penalty_temp0));
        penalty_temp0 = penalty_here14;
        penalty_here14 = (query_value-subject_value14) * (query_value-subject_value14) + min(penalty_here13, min(penalty_here14, penalty_temp1));
        penalty_temp1 = penalty_here15;
        penalty_here15 = (query_value-subject_value15) * (query_value-subject_value15) + min(penalty_here14, min(penalty_here15, penalty_temp0));

        penalty_temp0 = penalty_here16;
        penalty_here16 = (query_value-subject_value16) * (query_value-subject_value16) + min(penalty_here15, min(penalty_here16, penalty_temp1));
        penalty_temp1 = penalty_here17;
        penalty_here17 = (query_value-subject_value17) * (query_value-subject_value17) + min(penalty_here16, min(penalty_here17, penalty_temp0));
        penalty_temp0 = penalty_here18;
        penalty_here18 = (query_value-subject_value18) * (query_value-subject_value18) + min(penalty_here17, min(penalty_here18, penalty_temp1));
        penalty_temp1 = penalty_here19;
        penalty_here19 = (query_value-subject_value19) * (query_value-subject_value19) + min(penalty_here18, min(penalty_here19, penalty_temp0));
        penalty_temp0 = penalty_here20;
        penalty_here20 = (query_value-subject_value20) * (query_value-subject_value20) + min(penalty_here19, min(penalty_here20, penalty_temp1));
        penalty_temp1 = penalty_here21;
        penalty_here21 = (query_value-subject_value21) * (query_value-subject_value21) + min(penalty_here20, min(penalty_here21, penalty_temp0));
        penalty_temp0 = penalty_here22;
        penalty_here22 = (query_value-subject_value22) * (query_value-subject_value22) + min(penalty_here21, min(penalty_here22, penalty_temp1));
        penalty_temp1 = penalty_here23;
        penalty_here23 = (query_value-subject_value23) * (query_value-subject_value23) + min(penalty_here22, min(penalty_here23, penalty_temp0));
        penalty_temp0 = penalty_here24;
        penalty_here24 = (query_value-subject_value24) * (query_value-subject_value24) + min(penalty_here23, min(penalty_here24, penalty_temp1));
        penalty_temp1 = penalty_here25;
        penalty_here25 = (query_value-subject_value25) * (query_value-subject_value25) + min(penalty_here24, min(penalty_here25, penalty_temp0));
        penalty_temp0 = penalty_here26;
        penalty_here26 = (query_value-subject_value26) * (query_value-subject_value26) + min(penalty_here25, min(penalty_here26, penalty_temp1));
        penalty_temp1 = penalty_here27;
        penalty_here27 = (query_value-subject_value27) * (query_value-subject_value27) + min(penalty_here26, min(penalty_here27, penalty_temp0));
        penalty_temp0 = penalty_here28;
        penalty_here28 = (query_value-subject_value28) * (query_value-subject_value28) + min(penalty_here27, min(penalty_here28, penalty_temp1));
        penalty_temp1 = penalty_here29;
        penalty_here29 = (query_value-subject_value29) * (query_value-subject_value29) + min(penalty_here28, min(penalty_here29, penalty_temp0));
        penalty_temp0 = penalty_here30;
        penalty_here30 = (query_value-subject_value30) * (query_value-subject_value30) + min(penalty_here29, min(penalty_here30, penalty_temp1));
        penalty_temp1 = penalty_here31;
        penalty_here31 = (query_value-subject_value31) * (query_value-subject_value31) + min(penalty_here30, min(penalty_here31, penalty_temp0));

        penalty_temp0 = penalty_here32;
        penalty_here32 = (query_value-subject_value32) * (query_value-subject_value32) + min(penalty_here31, min(penalty_here32, penalty_temp1));
        penalty_temp1 = penalty_here33;
        penalty_here33 = (query_value-subject_value33) * (query_value-subject_value33) + min(penalty_here32, min(penalty_here33, penalty_temp0));
        penalty_temp0 = penalty_here34;
        penalty_here34 = (query_value-subject_value34) * (query_value-subject_value34) + min(penalty_here33, min(penalty_here34, penalty_temp1));
        penalty_temp1 = penalty_here35;
        penalty_here35 = (query_value-subject_value35) * (query_value-subject_value35) + min(penalty_here34, min(penalty_here35, penalty_temp0));
        penalty_temp0 = penalty_here36;
        penalty_here36 = (query_value-subject_value36) * (query_value-subject_value36) + min(penalty_here35, min(penalty_here36, penalty_temp1));
        penalty_temp1 = penalty_here37;
        penalty_here37 = (query_value-subject_value37) * (query_value-subject_value37) + min(penalty_here36, min(penalty_here37, penalty_temp0));
        penalty_temp0 = penalty_here38;
        penalty_here38 = (query_value-subject_value38) * (query_value-subject_value38) + min(penalty_here37, min(penalty_here38, penalty_temp1));
        penalty_temp1 = penalty_here39;
        penalty_here39 = (query_value-subject_value39) * (query_value-subject_value39) + min(penalty_here38, min(penalty_here39, penalty_temp0));
        penalty_temp0 = penalty_here40;
        penalty_here40 = (query_value-subject_value40) * (query_value-subject_value40) + min(penalty_here39, min(penalty_here40, penalty_temp1));
        penalty_temp1 = penalty_here41;
        penalty_here41 = (query_value-subject_value41) * (query_value-subject_value41) + min(penalty_here40, min(penalty_here41, penalty_temp0));
        penalty_temp0 = penalty_here42;
        penalty_here42 = (query_value-subject_value42) * (query_value-subject_value42) + min(penalty_here41, min(penalty_here42, penalty_temp1));
        penalty_temp1 = penalty_here43;
        penalty_here43 = (query_value-subject_value43) * (query_value-subject_value43) + min(penalty_here42, min(penalty_here43, penalty_temp0));
        penalty_temp0 = penalty_here44;
        penalty_here44 = (query_value-subject_value44) * (query_value-subject_value44) + min(penalty_here43, min(penalty_here44, penalty_temp1));
        penalty_temp1 = penalty_here45;
        penalty_here45 = (query_value-subject_value45) * (query_value-subject_value45) + min(penalty_here44, min(penalty_here45, penalty_temp0));
        penalty_temp0 = penalty_here46;
        penalty_here46 = (query_value-subject_value46) * (query_value-subject_value46) + min(penalty_here45, min(penalty_here46, penalty_temp1));
        penalty_temp1 = penalty_here47;
        penalty_here47 = (query_value-subject_value47) * (query_value-subject_value47) + min(penalty_here46, min(penalty_here47, penalty_temp0));

        penalty_temp0 = penalty_here48;
        penalty_here48 = (query_value-subject_value48) * (query_value-subject_value48) + min(penalty_here47, min(penalty_here48, penalty_temp1));
        penalty_temp1 = penalty_here49;
        penalty_here49 = (query_value-subject_value49) * (query_value-subject_value49) + min(penalty_here48, min(penalty_here49, penalty_temp0));
        penalty_temp0 = penalty_here50;
        penalty_here50 = (query_value-subject_value50) * (query_value-subject_value50) + min(penalty_here49, min(penalty_here50, penalty_temp1));
        penalty_temp1 = penalty_here51;
        penalty_here51 = (query_value-subject_value51) * (query_value-subject_value51) + min(penalty_here50, min(penalty_here51, penalty_temp0));
        penalty_temp0 = penalty_here52;
        penalty_here52 = (query_value-subject_value52) * (query_value-subject_value52) + min(penalty_here51, min(penalty_here52, penalty_temp1));
        penalty_temp1 = penalty_here53;
        penalty_here53 = (query_value-subject_value53) * (query_value-subject_value53) + min(penalty_here52, min(penalty_here53, penalty_temp0));
        penalty_temp0 = penalty_here54;
        penalty_here54 = (query_value-subject_value54) * (query_value-subject_value54) + min(penalty_here53, min(penalty_here54, penalty_temp1));
        penalty_temp1 = penalty_here55;
        penalty_here55 = (query_value-subject_value55) * (query_value-subject_value55) + min(penalty_here54, min(penalty_here55, penalty_temp0));
        penalty_temp0 = penalty_here56;
        penalty_here56 = (query_value-subject_value56) * (query_value-subject_value56) + min(penalty_here55, min(penalty_here56, penalty_temp1));
        penalty_temp1 = penalty_here57;
        penalty_here57 = (query_value-subject_value57) * (query_value-subject_value57) + min(penalty_here56, min(penalty_here57, penalty_temp0));
        penalty_temp0 = penalty_here58;
        penalty_here58 = (query_value-subject_value58) * (query_value-subject_value58) + min(penalty_here57, min(penalty_here58, penalty_temp1));
        penalty_temp1 = penalty_here59;
        penalty_here59 = (query_value-subject_value59) * (query_value-subject_value59) + min(penalty_here58, min(penalty_here59, penalty_temp0));
        penalty_temp0 = penalty_here60;
        penalty_here60 = (query_value-subject_value60) * (query_value-subject_value60) + min(penalty_here59, min(penalty_here60, penalty_temp1));
        penalty_temp1 = penalty_here61;
        penalty_here61 = (query_value-subject_value61) * (query_value-subject_value61) + min(penalty_here60, min(penalty_here61, penalty_temp0));
        penalty_temp0 = penalty_here62;
        penalty_here62 = (query_value-subject_value62) * (query_value-subject_value62) + min(penalty_here61, min(penalty_here62, penalty_temp1));
        penalty_here63 = (query_value-subject_value63) * (query_value-subject_value63) + min(penalty_here62, min(penalty_here63, penalty_temp0));


    query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
    if (thid == 0) query_value = new_query_value;
    new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
    counter++;

    penalty_diag = penalty_left;
    penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here31, 1, 32);

    if (thid == 0) penalty_left = INFINITY;

    for (index_t k = 3; k < lane+32-1; k++) {
        const index_t i = k-l;
            //outside = k <= l || i >= lane;

        penalty_temp0 = penalty_here0;
            penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
            penalty_temp1 = penalty_here1;
            penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
            penalty_temp0 = penalty_here2;
            penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
            penalty_temp1 = penalty_here3;
            penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));
            penalty_temp0 = penalty_here4;
            penalty_here4 = (query_value-subject_value4) * (query_value-subject_value4) + min(penalty_here3, min(penalty_here4, penalty_temp1));
            penalty_temp1 = penalty_here5;
            penalty_here5 = (query_value-subject_value5) * (query_value-subject_value5) + min(penalty_here4, min(penalty_here5, penalty_temp0));
            penalty_temp0 = penalty_here6;
            penalty_here6 = (query_value-subject_value6) * (query_value-subject_value6) + min(penalty_here5, min(penalty_here6, penalty_temp1));
            penalty_temp1 = penalty_here7;
            penalty_here7 = (query_value-subject_value7) * (query_value-subject_value7) + min(penalty_here6, min(penalty_here7, penalty_temp0));
            penalty_temp0 = penalty_here8;
            penalty_here8 = (query_value-subject_value8) * (query_value-subject_value8) + min(penalty_here7, min(penalty_here8, penalty_temp1));
            penalty_temp1 = penalty_here9;
            penalty_here9 = (query_value-subject_value9) * (query_value-subject_value9) + min(penalty_here8, min(penalty_here9, penalty_temp0));
            penalty_temp0 = penalty_here10;
            penalty_here10 = (query_value-subject_value10) * (query_value-subject_value10) + min(penalty_here9, min(penalty_here10, penalty_temp1));
            penalty_temp1 = penalty_here11;
            penalty_here11 = (query_value-subject_value11) * (query_value-subject_value11) + min(penalty_here10, min(penalty_here11, penalty_temp0));
            penalty_temp0 = penalty_here12;
            penalty_here12 = (query_value-subject_value12) * (query_value-subject_value12) + min(penalty_here11, min(penalty_here12, penalty_temp1));
            penalty_temp1 = penalty_here13;
            penalty_here13 = (query_value-subject_value13) * (query_value-subject_value13) + min(penalty_here12, min(penalty_here13, penalty_temp0));
            penalty_temp0 = penalty_here14;
            penalty_here14 = (query_value-subject_value14) * (query_value-subject_value14) + min(penalty_here13, min(penalty_here14, penalty_temp1));
            penalty_temp1 = penalty_here15;
            penalty_here15 = (query_value-subject_value15) * (query_value-subject_value15) + min(penalty_here14, min(penalty_here15, penalty_temp0));

            penalty_temp0 = penalty_here16;
            penalty_here16 = (query_value-subject_value16) * (query_value-subject_value16) + min(penalty_here15, min(penalty_here16, penalty_temp1));
            penalty_temp1 = penalty_here17;
            penalty_here17 = (query_value-subject_value17) * (query_value-subject_value17) + min(penalty_here16, min(penalty_here17, penalty_temp0));
            penalty_temp0 = penalty_here18;
            penalty_here18 = (query_value-subject_value18) * (query_value-subject_value18) + min(penalty_here17, min(penalty_here18, penalty_temp1));
            penalty_temp1 = penalty_here19;
            penalty_here19 = (query_value-subject_value19) * (query_value-subject_value19) + min(penalty_here18, min(penalty_here19, penalty_temp0));
            penalty_temp0 = penalty_here20;
            penalty_here20 = (query_value-subject_value20) * (query_value-subject_value20) + min(penalty_here19, min(penalty_here20, penalty_temp1));
            penalty_temp1 = penalty_here21;
            penalty_here21 = (query_value-subject_value21) * (query_value-subject_value21) + min(penalty_here20, min(penalty_here21, penalty_temp0));
            penalty_temp0 = penalty_here22;
            penalty_here22 = (query_value-subject_value22) * (query_value-subject_value22) + min(penalty_here21, min(penalty_here22, penalty_temp1));
            penalty_temp1 = penalty_here23;
            penalty_here23 = (query_value-subject_value23) * (query_value-subject_value23) + min(penalty_here22, min(penalty_here23, penalty_temp0));
            penalty_temp0 = penalty_here24;
            penalty_here24 = (query_value-subject_value24) * (query_value-subject_value24) + min(penalty_here23, min(penalty_here24, penalty_temp1));
            penalty_temp1 = penalty_here25;
            penalty_here25 = (query_value-subject_value25) * (query_value-subject_value25) + min(penalty_here24, min(penalty_here25, penalty_temp0));
            penalty_temp0 = penalty_here26;
            penalty_here26 = (query_value-subject_value26) * (query_value-subject_value26) + min(penalty_here25, min(penalty_here26, penalty_temp1));
            penalty_temp1 = penalty_here27;
            penalty_here27 = (query_value-subject_value27) * (query_value-subject_value27) + min(penalty_here26, min(penalty_here27, penalty_temp0));
            penalty_temp0 = penalty_here28;
            penalty_here28 = (query_value-subject_value28) * (query_value-subject_value28) + min(penalty_here27, min(penalty_here28, penalty_temp1));
            penalty_temp1 = penalty_here29;
            penalty_here29 = (query_value-subject_value29) * (query_value-subject_value29) + min(penalty_here28, min(penalty_here29, penalty_temp0));
            penalty_temp0 = penalty_here30;
            penalty_here30 = (query_value-subject_value30) * (query_value-subject_value30) + min(penalty_here29, min(penalty_here30, penalty_temp1));
            penalty_temp1 = penalty_here31;
            penalty_here31 = (query_value-subject_value31) * (query_value-subject_value31) + min(penalty_here30, min(penalty_here31, penalty_temp0));

            penalty_temp0 = penalty_here32;
            penalty_here32 = (query_value-subject_value32) * (query_value-subject_value32) + min(penalty_here31, min(penalty_here32, penalty_temp1));
            penalty_temp1 = penalty_here33;
            penalty_here33 = (query_value-subject_value33) * (query_value-subject_value33) + min(penalty_here32, min(penalty_here33, penalty_temp0));
            penalty_temp0 = penalty_here34;
            penalty_here34 = (query_value-subject_value34) * (query_value-subject_value34) + min(penalty_here33, min(penalty_here34, penalty_temp1));
            penalty_temp1 = penalty_here35;
            penalty_here35 = (query_value-subject_value35) * (query_value-subject_value35) + min(penalty_here34, min(penalty_here35, penalty_temp0));
            penalty_temp0 = penalty_here36;
            penalty_here36 = (query_value-subject_value36) * (query_value-subject_value36) + min(penalty_here35, min(penalty_here36, penalty_temp1));
            penalty_temp1 = penalty_here37;
            penalty_here37 = (query_value-subject_value37) * (query_value-subject_value37) + min(penalty_here36, min(penalty_here37, penalty_temp0));
            penalty_temp0 = penalty_here38;
            penalty_here38 = (query_value-subject_value38) * (query_value-subject_value38) + min(penalty_here37, min(penalty_here38, penalty_temp1));
            penalty_temp1 = penalty_here39;
            penalty_here39 = (query_value-subject_value39) * (query_value-subject_value39) + min(penalty_here38, min(penalty_here39, penalty_temp0));
            penalty_temp0 = penalty_here40;
            penalty_here40 = (query_value-subject_value40) * (query_value-subject_value40) + min(penalty_here39, min(penalty_here40, penalty_temp1));
            penalty_temp1 = penalty_here41;
            penalty_here41 = (query_value-subject_value41) * (query_value-subject_value41) + min(penalty_here40, min(penalty_here41, penalty_temp0));
            penalty_temp0 = penalty_here42;
            penalty_here42 = (query_value-subject_value42) * (query_value-subject_value42) + min(penalty_here41, min(penalty_here42, penalty_temp1));
            penalty_temp1 = penalty_here43;
            penalty_here43 = (query_value-subject_value43) * (query_value-subject_value43) + min(penalty_here42, min(penalty_here43, penalty_temp0));
            penalty_temp0 = penalty_here44;
            penalty_here44 = (query_value-subject_value44) * (query_value-subject_value44) + min(penalty_here43, min(penalty_here44, penalty_temp1));
            penalty_temp1 = penalty_here45;
            penalty_here45 = (query_value-subject_value45) * (query_value-subject_value45) + min(penalty_here44, min(penalty_here45, penalty_temp0));
            penalty_temp0 = penalty_here46;
            penalty_here46 = (query_value-subject_value46) * (query_value-subject_value46) + min(penalty_here45, min(penalty_here46, penalty_temp1));
            penalty_temp1 = penalty_here47;
            penalty_here47 = (query_value-subject_value47) * (query_value-subject_value47) + min(penalty_here46, min(penalty_here47, penalty_temp0));

            penalty_temp0 = penalty_here48;
            penalty_here48 = (query_value-subject_value48) * (query_value-subject_value48) + min(penalty_here47, min(penalty_here48, penalty_temp1));
            penalty_temp1 = penalty_here49;
            penalty_here49 = (query_value-subject_value49) * (query_value-subject_value49) + min(penalty_here48, min(penalty_here49, penalty_temp0));
            penalty_temp0 = penalty_here50;
            penalty_here50 = (query_value-subject_value50) * (query_value-subject_value50) + min(penalty_here49, min(penalty_here50, penalty_temp1));
            penalty_temp1 = penalty_here51;
            penalty_here51 = (query_value-subject_value51) * (query_value-subject_value51) + min(penalty_here50, min(penalty_here51, penalty_temp0));
            penalty_temp0 = penalty_here52;
            penalty_here52 = (query_value-subject_value52) * (query_value-subject_value52) + min(penalty_here51, min(penalty_here52, penalty_temp1));
            penalty_temp1 = penalty_here53;
            penalty_here53 = (query_value-subject_value53) * (query_value-subject_value53) + min(penalty_here52, min(penalty_here53, penalty_temp0));
            penalty_temp0 = penalty_here54;
            penalty_here54 = (query_value-subject_value54) * (query_value-subject_value54) + min(penalty_here53, min(penalty_here54, penalty_temp1));
            penalty_temp1 = penalty_here55;
            penalty_here55 = (query_value-subject_value55) * (query_value-subject_value55) + min(penalty_here54, min(penalty_here55, penalty_temp0));
            penalty_temp0 = penalty_here56;
            penalty_here56 = (query_value-subject_value56) * (query_value-subject_value56) + min(penalty_here55, min(penalty_here56, penalty_temp1));
            penalty_temp1 = penalty_here57;
            penalty_here57 = (query_value-subject_value57) * (query_value-subject_value57) + min(penalty_here56, min(penalty_here57, penalty_temp0));
            penalty_temp0 = penalty_here58;
            penalty_here58 = (query_value-subject_value58) * (query_value-subject_value58) + min(penalty_here57, min(penalty_here58, penalty_temp1));
            penalty_temp1 = penalty_here59;
            penalty_here59 = (query_value-subject_value59) * (query_value-subject_value59) + min(penalty_here58, min(penalty_here59, penalty_temp0));
            penalty_temp0 = penalty_here60;
            penalty_here60 = (query_value-subject_value60) * (query_value-subject_value60) + min(penalty_here59, min(penalty_here60, penalty_temp1));
            penalty_temp1 = penalty_here61;
            penalty_here61 = (query_value-subject_value61) * (query_value-subject_value61) + min(penalty_here60, min(penalty_here61, penalty_temp0));
            penalty_temp0 = penalty_here62;
            penalty_here62 = (query_value-subject_value62) * (query_value-subject_value62) + min(penalty_here61, min(penalty_here62, penalty_temp1));
            penalty_here63 = (query_value-subject_value63) * (query_value-subject_value63) + min(penalty_here62, min(penalty_here63, penalty_temp0));

            //if (blid == 0 && thid == 31 && iter == 0) Dist[2*(k-1)] = penalty_here0;
            //if (blid == 0 && thid == 31 && iter == 0) Dist[2*(k-1)+1] = penalty_here1;
            //if (blid == 0 && iter == 0) Dist[64*(k-1)+2*thid] = penalty_here0;
            //if (blid == 0 && iter == 0) Dist[64*(k-1)+2*thid+1] = penalty_here1;

        if (counter%32 == 0) {
            penalty_temp0 = min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(min(penalty_left, penalty_here0), penalty_here1), penalty_here2), penalty_here3), penalty_here4), penalty_here5), penalty_here6), penalty_here7),
            penalty_here8), penalty_here9), penalty_here10), penalty_here11), penalty_here12), penalty_here13), penalty_here14), penalty_here15),
            penalty_here16), penalty_here17), penalty_here18), penalty_here19), penalty_here20), penalty_here21), penalty_here22), penalty_here23),
            penalty_here24), penalty_here25), penalty_here26), penalty_here27), penalty_here28), penalty_here29), penalty_here30), penalty_here31);
            mask = __any_sync(0xFFFFFFFF, penalty_temp0 < *threshold - lower_right);
            if (!mask) {Dist[blid] = 1000*(counter/32); break; }

            new_query_value = cQuery[i+2*thid-1];
        }
        query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
        if (thid == 0) query_value = new_query_value;
        new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
            //if (thid == 0) if (!outside) Dist[counter] = query_value; else Dist[counter] = 0;
        counter++;

            // shuffle the penalty
        penalty_diag = penalty_left;
        penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here63, 1, 32);
            //if (thid == 0 && !outside) penalty_left = Subject_cache[i+1]; // TO DO: replace by shuffles
        //    if (iter > 0 && thid == 0 && k>l) penalty_left = Subject_cache[i+1];
            //if (iter && thid == 0) penalty_left = Subject_cache[i+1];
            //if (!iter && thid == 0) penalty_left = INFINITY;
        if (thid == 0) penalty_left = INFINITY;
    }
    penalty_temp0 = penalty_here0;
    penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
        penalty_temp1 = penalty_here1;
        penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
        penalty_temp0 = penalty_here2;
        penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
        penalty_temp1 = penalty_here3;
        penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));
        penalty_temp0 = penalty_here4;
        penalty_here4 = (query_value-subject_value4) * (query_value-subject_value4) + min(penalty_here3, min(penalty_here4, penalty_temp1));
        penalty_temp1 = penalty_here5;
        penalty_here5 = (query_value-subject_value5) * (query_value-subject_value5) + min(penalty_here4, min(penalty_here5, penalty_temp0));
        penalty_temp0 = penalty_here6;
        penalty_here6 = (query_value-subject_value6) * (query_value-subject_value6) + min(penalty_here5, min(penalty_here6, penalty_temp1));
        penalty_temp1 = penalty_here7;
        penalty_here7 = (query_value-subject_value7) * (query_value-subject_value7) + min(penalty_here6, min(penalty_here7, penalty_temp0));
        penalty_temp0 = penalty_here8;
        penalty_here8 = (query_value-subject_value8) * (query_value-subject_value8) + min(penalty_here7, min(penalty_here8, penalty_temp1));
        penalty_temp1 = penalty_here9;
        penalty_here9 = (query_value-subject_value9) * (query_value-subject_value9) + min(penalty_here8, min(penalty_here9, penalty_temp0));
        penalty_temp0 = penalty_here10;
        penalty_here10 = (query_value-subject_value10) * (query_value-subject_value10) + min(penalty_here9, min(penalty_here10, penalty_temp1));
        penalty_temp1 = penalty_here11;
        penalty_here11 = (query_value-subject_value11) * (query_value-subject_value11) + min(penalty_here10, min(penalty_here11, penalty_temp0));
        penalty_temp0 = penalty_here12;
        penalty_here12 = (query_value-subject_value12) * (query_value-subject_value12) + min(penalty_here11, min(penalty_here12, penalty_temp1));
        penalty_temp1 = penalty_here13;
        penalty_here13 = (query_value-subject_value13) * (query_value-subject_value13) + min(penalty_here12, min(penalty_here13, penalty_temp0));
        penalty_temp0 = penalty_here14;
        penalty_here14 = (query_value-subject_value14) * (query_value-subject_value14) + min(penalty_here13, min(penalty_here14, penalty_temp1));
        penalty_temp1 = penalty_here15;
        penalty_here15 = (query_value-subject_value15) * (query_value-subject_value15) + min(penalty_here14, min(penalty_here15, penalty_temp0));

        penalty_temp0 = penalty_here16;
        penalty_here16 = (query_value-subject_value16) * (query_value-subject_value16) + min(penalty_here15, min(penalty_here16, penalty_diag));
        penalty_temp1 = penalty_here17;
        penalty_here17 = (query_value-subject_value17) * (query_value-subject_value17) + min(penalty_here16, min(penalty_here17, penalty_temp0));
        penalty_temp0 = penalty_here18;
        penalty_here18 = (query_value-subject_value18) * (query_value-subject_value18) + min(penalty_here17, min(penalty_here18, penalty_temp1));
        penalty_temp1 = penalty_here19;
        penalty_here19 = (query_value-subject_value19) * (query_value-subject_value19) + min(penalty_here18, min(penalty_here19, penalty_temp0));
        penalty_temp0 = penalty_here20;
        penalty_here20 = (query_value-subject_value20) * (query_value-subject_value20) + min(penalty_here19, min(penalty_here20, penalty_temp1));
        penalty_temp1 = penalty_here21;
        penalty_here21 = (query_value-subject_value21) * (query_value-subject_value21) + min(penalty_here20, min(penalty_here21, penalty_temp0));
        penalty_temp0 = penalty_here22;
        penalty_here22 = (query_value-subject_value22) * (query_value-subject_value22) + min(penalty_here21, min(penalty_here22, penalty_temp1));
        penalty_temp1 = penalty_here23;
        penalty_here23 = (query_value-subject_value23) * (query_value-subject_value23) + min(penalty_here22, min(penalty_here23, penalty_temp0));
        penalty_temp0 = penalty_here24;
        penalty_here24 = (query_value-subject_value24) * (query_value-subject_value24) + min(penalty_here23, min(penalty_here24, penalty_temp1));
        penalty_temp1 = penalty_here25;
        penalty_here25 = (query_value-subject_value25) * (query_value-subject_value25) + min(penalty_here24, min(penalty_here25, penalty_temp0));
        penalty_temp0 = penalty_here26;
        penalty_here26 = (query_value-subject_value26) * (query_value-subject_value26) + min(penalty_here25, min(penalty_here26, penalty_temp1));
        penalty_temp1 = penalty_here27;
        penalty_here27 = (query_value-subject_value27) * (query_value-subject_value27) + min(penalty_here26, min(penalty_here27, penalty_temp0));
        penalty_temp0 = penalty_here28;
        penalty_here28 = (query_value-subject_value28) * (query_value-subject_value28) + min(penalty_here27, min(penalty_here28, penalty_temp1));
        penalty_temp1 = penalty_here29;
        penalty_here29 = (query_value-subject_value29) * (query_value-subject_value29) + min(penalty_here28, min(penalty_here29, penalty_temp0));
        penalty_temp0 = penalty_here30;
        penalty_here30 = (query_value-subject_value30) * (query_value-subject_value30) + min(penalty_here29, min(penalty_here30, penalty_temp1));
        penalty_temp1 = penalty_here31;
        penalty_here31 = (query_value-subject_value31) * (query_value-subject_value31) + min(penalty_here30, min(penalty_here31, penalty_temp0));

        penalty_temp0 = penalty_here32;
        penalty_here32 = (query_value-subject_value32) * (query_value-subject_value32) + min(penalty_here31, min(penalty_here32, penalty_temp1));
        penalty_temp1 = penalty_here33;
        penalty_here33 = (query_value-subject_value33) * (query_value-subject_value33) + min(penalty_here32, min(penalty_here33, penalty_temp0));
        penalty_temp0 = penalty_here34;
        penalty_here34 = (query_value-subject_value34) * (query_value-subject_value34) + min(penalty_here33, min(penalty_here34, penalty_temp1));
        penalty_temp1 = penalty_here35;
        penalty_here35 = (query_value-subject_value35) * (query_value-subject_value35) + min(penalty_here34, min(penalty_here35, penalty_temp0));
        penalty_temp0 = penalty_here36;
        penalty_here36 = (query_value-subject_value36) * (query_value-subject_value36) + min(penalty_here35, min(penalty_here36, penalty_temp1));
        penalty_temp1 = penalty_here37;
        penalty_here37 = (query_value-subject_value37) * (query_value-subject_value37) + min(penalty_here36, min(penalty_here37, penalty_temp0));
        penalty_temp0 = penalty_here38;
        penalty_here38 = (query_value-subject_value38) * (query_value-subject_value38) + min(penalty_here37, min(penalty_here38, penalty_temp1));
        penalty_temp1 = penalty_here39;
        penalty_here39 = (query_value-subject_value39) * (query_value-subject_value39) + min(penalty_here38, min(penalty_here39, penalty_temp0));
        penalty_temp0 = penalty_here40;
        penalty_here40 = (query_value-subject_value40) * (query_value-subject_value40) + min(penalty_here39, min(penalty_here40, penalty_temp1));
        penalty_temp1 = penalty_here41;
        penalty_here41 = (query_value-subject_value41) * (query_value-subject_value41) + min(penalty_here40, min(penalty_here41, penalty_temp0));
        penalty_temp0 = penalty_here42;
        penalty_here42 = (query_value-subject_value42) * (query_value-subject_value42) + min(penalty_here41, min(penalty_here42, penalty_temp1));
        penalty_temp1 = penalty_here43;
        penalty_here43 = (query_value-subject_value43) * (query_value-subject_value43) + min(penalty_here42, min(penalty_here43, penalty_temp0));
        penalty_temp0 = penalty_here44;
        penalty_here44 = (query_value-subject_value44) * (query_value-subject_value44) + min(penalty_here43, min(penalty_here44, penalty_temp1));
        penalty_temp1 = penalty_here45;
        penalty_here45 = (query_value-subject_value45) * (query_value-subject_value45) + min(penalty_here44, min(penalty_here45, penalty_temp0));
        penalty_temp0 = penalty_here46;
        penalty_here46 = (query_value-subject_value46) * (query_value-subject_value46) + min(penalty_here45, min(penalty_here46, penalty_temp1));
        penalty_temp1 = penalty_here47;
        penalty_here47 = (query_value-subject_value47) * (query_value-subject_value47) + min(penalty_here46, min(penalty_here47, penalty_temp0));

        penalty_temp0 = penalty_here48;
        penalty_here48 = (query_value-subject_value48) * (query_value-subject_value48) + min(penalty_here47, min(penalty_here48, penalty_temp1));
        penalty_temp1 = penalty_here49;
        penalty_here49 = (query_value-subject_value49) * (query_value-subject_value49) + min(penalty_here48, min(penalty_here49, penalty_temp0));
        penalty_temp0 = penalty_here50;
        penalty_here50 = (query_value-subject_value50) * (query_value-subject_value50) + min(penalty_here49, min(penalty_here50, penalty_temp1));
        penalty_temp1 = penalty_here51;
        penalty_here51 = (query_value-subject_value51) * (query_value-subject_value51) + min(penalty_here50, min(penalty_here51, penalty_temp0));
        penalty_temp0 = penalty_here52;
        penalty_here52 = (query_value-subject_value52) * (query_value-subject_value52) + min(penalty_here51, min(penalty_here52, penalty_temp1));
        penalty_temp1 = penalty_here53;
        penalty_here53 = (query_value-subject_value53) * (query_value-subject_value53) + min(penalty_here52, min(penalty_here53, penalty_temp0));
        penalty_temp0 = penalty_here54;
        penalty_here54 = (query_value-subject_value54) * (query_value-subject_value54) + min(penalty_here53, min(penalty_here54, penalty_temp1));
        penalty_temp1 = penalty_here55;
        penalty_here55 = (query_value-subject_value55) * (query_value-subject_value55) + min(penalty_here54, min(penalty_here55, penalty_temp0));
        penalty_temp0 = penalty_here56;
        penalty_here56 = (query_value-subject_value56) * (query_value-subject_value56) + min(penalty_here55, min(penalty_here56, penalty_temp1));
        penalty_temp1 = penalty_here57;
        penalty_here57 = (query_value-subject_value57) * (query_value-subject_value57) + min(penalty_here56, min(penalty_here57, penalty_temp0));
        penalty_temp0 = penalty_here58;
        penalty_here58 = (query_value-subject_value58) * (query_value-subject_value58) + min(penalty_here57, min(penalty_here58, penalty_temp1));
        penalty_temp1 = penalty_here59;
        penalty_here59 = (query_value-subject_value59) * (query_value-subject_value59) + min(penalty_here58, min(penalty_here59, penalty_temp0));
        penalty_temp0 = penalty_here60;
        penalty_here60 = (query_value-subject_value60) * (query_value-subject_value60) + min(penalty_here59, min(penalty_here60, penalty_temp1));
        penalty_temp1 = penalty_here61;
        penalty_here61 = (query_value-subject_value61) * (query_value-subject_value61) + min(penalty_here60, min(penalty_here61, penalty_temp0));
        penalty_temp0 = penalty_here62;
        penalty_here62 = (query_value-subject_value62) * (query_value-subject_value62) + min(penalty_here61, min(penalty_here62, penalty_temp1));
        penalty_here63 = (query_value-subject_value63) * (query_value-subject_value63) + min(penalty_here62, min(penalty_here63, penalty_temp0));



    if(thid == blockDim.x-1 && mask)  {
        Dist[blid] = penalty_here63;
        if (penalty_here63 < *threshold)
               update_bsf_system(threshold, penalty_here63);
    }
}
}


#endif
