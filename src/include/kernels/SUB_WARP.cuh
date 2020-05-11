#ifndef SUB_WARP_DTW
#define SUB_WARP_DTW


// 32 values per thread no shared memory
template <
    int group_size,
    typename index_t,
    typename value_t> __global__
void sub_warp_DTW(  // replace DTW_1024_sub_warp
    value_t * Subject,
    value_t * Dist,
    index_t num_entries,
    index_t num_features) {

    const index_t blid = blockIdx.x;
    const index_t thid = threadIdx.x;
    const index_t l = thid;
    const index_t lane = num_features+1;
    //const index_t WARP_SIZE = 32;
    //const index_t group_size = WARP_SIZE/(1024/lane);
    const index_t base = (32/group_size)*blid*num_features;

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
    value_t penalty_temp0;
    value_t penalty_temp1;

        // Init shared memeory for right column
    //for (index_t l = thid; l < lane; l += blockDim.x)
    //    Subject_cache[l] = INFINITY;
    //__syncthreads();

    //index_t iter = 0;

    if (thid % group_size == 0) {
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

    //if (thid >= 8 && thid < 16) {
    //    subject_value0 = Subject[base+num_features+16*(l-8)-1];
    //const value_t subject_value0 = Subject[base+(thid/group_size)*num_features + 32*(l-(group_size*(thid/group_size)))-1];

    //const value_t subject_value0 = base+l == 0 ? 0 : Subject[base+32*l-1];
    //const value_t subject_value1 = Subject[base+32*l-0];
    //const value_t subject_value2 = Subject[base+32*l+1];

    const value_t subject_value0 = base+(thid/group_size)*num_features + (l%group_size) == 0 ? 0 : Subject[base+(thid/group_size)*num_features + 32*(l%group_size)-1];
    //const value_t subject_value0 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)-1];
    const value_t subject_value1 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)-0];
    const value_t subject_value2 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+1];
    const value_t subject_value3 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+2];
    const value_t subject_value4 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+3];
    const value_t subject_value5 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+4];
    const value_t subject_value6 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+5];
    const value_t subject_value7 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+6];
    const value_t subject_value8 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+7];
    const value_t subject_value9 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+8];
    const value_t subject_value10 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+9];
    const value_t subject_value11 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+10];
    const value_t subject_value12 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+11];
    const value_t subject_value13 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+12];
    const value_t subject_value14 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+13];
    const value_t subject_value15 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+14];

    const value_t subject_value16 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+15];
    const value_t subject_value17 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+16];
    const value_t subject_value18 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+17];
    const value_t subject_value19 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+18];
    const value_t subject_value20 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+19];
    const value_t subject_value21 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+20];
    const value_t subject_value22 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+21];
    const value_t subject_value23 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+22];
    const value_t subject_value24 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+23];
    const value_t subject_value25 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+24];
    const value_t subject_value26 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+25];
    const value_t subject_value27 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+26];
    const value_t subject_value28 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+27];
    const value_t subject_value29 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+28];
    const value_t subject_value30 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+29];
    const value_t subject_value31 = Subject[base+(thid/group_size)*num_features + 32*(l%group_size)+30];
        //if (blid == 0) Dist[2*l] = subject_value0;
        //if (blid == 0) Dist[2*l+1] = subject_value1;
    index_t counter = 1;
    value_t query_value = INFINITY;
    value_t new_query_value = cQuery[thid%group_size];
    if (thid % group_size == 0) query_value = new_query_value;
    if (thid % group_size == 0) penalty_here1 = 0; // (query_value - subject_value0)*(query_value - subject_value0);
    new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);

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
    penalty_here31 = (query_value-subject_value31) * (query_value-subject_value31) + min(penalty_here30, min(penalty_here31, penalty_temp0));

    query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
    if (thid % group_size == 0) query_value = new_query_value;
    new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
    counter++;

    penalty_diag = penalty_left;
    penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here31, 1, 32);

    //penalty_left = INFINITY;

    if (thid % group_size == 0)  penalty_left = INFINITY;

   const index_t group_id = thid%group_size;

    for (index_t k = 3; k < lane+group_size-1; k++) {
        const index_t i = k-l%group_size;
            //outside = k <= l || i >= lane;

            //const value_t residue = outside ? INFINITY : cQuery[i-1]-subject_value;
            //const value_t residue = outside ? INFINITY : query_value-subject_value;
            //if (thid == 0 && iter == 0 && k == 2) penalty_temp = INFINITY; else
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
            penalty_here31 = (query_value-subject_value31) * (query_value-subject_value31) + min(penalty_here30, min(penalty_here31, penalty_temp0));

            //if (blid == 0 && thid == 31 && iter == 0) Dist[2*(k-1)] = penalty_here0;
            //if (blid == 0 && thid == 31 && iter == 0) Dist[2*(k-1)+1] = penalty_here1;
            //if (blid == 0 && iter == 0) Dist[64*(k-1)+2*thid] = penalty_here0;
            //if (blid == 0 && iter == 0) Dist[64*(k-1)+2*thid+1] = penalty_here1;

            //if (counter%32 == 0 && counter > 1) new_query_value = cQuery[i+2*thid-1];
        if (counter%group_size == 0) new_query_value = cQuery[i+2*(thid%group_size)-1];
        query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
        if (!group_id) query_value = new_query_value;
        new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
        counter++;

            // save the right column
            //if (!outside && thid == 31 && iter < lane/WARP_SIZE-1) Subject_cache[i] = penalty_here; // TO DO: replace this by shhffles
            //if (iter < lane/WARP_SIZE-1 && thid == 31 && k>l) Subject_cache[i] = penalty_here31;

            // shuffle the penalty
        penalty_diag = penalty_left;
        penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here31, 1, 32);
            //if (thid == 0 && !outside) penalty_left = Subject_cache[i+1]; // TO DO: replace by shuffles
        //    if (iter > 0 && thid == 0 && k>l) penalty_left = Subject_cache[i+1];
            //if (iter && thid == 0) penalty_left = Subject_cache[i+1];
            //if (!iter && thid == 0) penalty_left = INFINITY;
        if (!group_id) penalty_left = INFINITY;

            //if (thid == 0 && k>l) if (iter > 0) penalty_left = Subject_cache[i+1]; else penalty_left = INFINITY;
            //if (thid == 0 && k>l && iter > 0) penalty_left = Subject_cache[i+1];
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
        penalty_here31 = (query_value-subject_value31) * (query_value-subject_value31) + min(penalty_here30, min(penalty_here31, penalty_temp0));

    if(thid % group_size == group_size-1) Dist[(32/group_size)*blid+thid/group_size] = penalty_here31;
}


#endif
