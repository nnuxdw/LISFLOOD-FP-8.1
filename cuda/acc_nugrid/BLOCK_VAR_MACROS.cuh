#pragma once

#ifndef THREADS_PER_BLOCK_MRA
    // MUST BE EITHER 64, 256 OR 1024
    #define THREADS_PER_BLOCK_MRA 256 

    #if   THREADS_PER_BLOCK_MRA == 64
        #define SHARED_MEMORY_BLOCK_DIM 8
        #define LVL_SINGLE_BLOCK 3

    #elif THREADS_PER_BLOCK_MRA == 256
        #define SHARED_MEMORY_BLOCK_DIM 16
        #define LVL_SINGLE_BLOCK 4

    #elif THREADS_PER_BLOCK_MRA == 1024
        #define SHARED_MEMORY_BLOCK_DIM 32
        #define LVL_SINGLE_BLOCK 5

    #endif
#endif

#ifndef NUM_VARS
    #define NUM_VARS 4
#endif

#ifndef NUM_CHILDREN
   #define NUM_CHILDREN 4
#endif

#ifndef SIGNIFICANT
    #define SIGNIFICANT true
#endif

#ifndef INSIGNIFICANT
    #define INSIGNIFICANT false
#endif