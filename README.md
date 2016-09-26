# DataBase-Systems-Implementation-Part-2-3

In this part of the project, we re-implemented the probe routine in C using Intel’s SSE instruction set. Our code more closely
matched the code fragments from the reference document.
Our program should be invokable using the same interface as before, but we only support fanouts of 5, 9, and 17 (i.e., 4, 8, 
16 keys, respectively). The output from parts 1 and 2 should be identical. While our method is fully general, we also hardcoded
one particular tree structure, namely the one with fanouts of 9, 5, and 9 on successive levels (see the reference document).
Hardcoding means that certain overheads such as dereferencing function pointers or checking the fanout at each level can be 
avoided during probes. The hardcoded probe function for a given key has no if or while or for statements. Almost all of the 
statements in the probe routine use SSE intrinsics. The hardcoded version incorporates two additional optimizations. First, 
the root node of the index is explicitly loaded into register variables so that it is not re-read from the array for each 
search. Second, probes should are done four at a time. To load and broadcast the 4 keys from the input to register
variables, we use:
__m128i k = _mm_load_si128((__m128i*) &probe_keys[i]);
register __m128i k1 = _mm_shuffle_epi32(k, _MM_SHUFFLE(0,0,0,0));
register __m128i k2 = _mm_shuffle_epi32(k, _MM_SHUFFLE(1,1,1,1));
register __m128i k3 = _mm_shuffle_epi32(k, _MM_SHUFFLE(2,2,2,2));
register __m128i k4 = _mm_shuffle_epi32(k, _MM_SHUFFLE(3,3,3,3));
These 4 keys are processed in an unrolled fashion: access level 1 for all four, then access level 2 for all four etc. 
This optimization allows work to happen for one probe while another has a data access stall, for example.
Again, we implement the method in three distinct phases to facilitate timing measurements. In phase 1, the index is built, and
the probes loaded into an array of integers. In phase 2, the index is probed using the more advanced method and the range 
identifier of the match is written to an array of output values. In phase 3, the output value array is written to stdout.
