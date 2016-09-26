#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include "p2random.h"
#include "tree.h"


#include <malloc.h>
#include <limits.h>
#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#define BILLION 1000000000L
void print128_num(__m128i var)
{
    uint32_t *val = (uint32_t*) &var;
    printf("Numerical: %i %i %i %i \n", 
           val[0], val[1], val[2], val[3]);
}

int main(int argc, char* argv[]) {
        // parsing arguments
        assert(argc > 3);
        size_t num_keys = strtoull(argv[1], NULL, 0);
        size_t num_probes = strtoull(argv[2], NULL, 0);
        
        size_t num_levels = (size_t) argc - 3;
        size_t* fanout = malloc(sizeof(size_t) * num_levels);
        assert(fanout != NULL);
        for (size_t i = 0; i < num_levels; ++i) {
                fanout[i] = strtoull(argv[i + 3], NULL, 0);
                assert(fanout[i]==5 || fanout[i] ==9 || fanout[i]==17);
                assert(fanout[i] >= 2 && fanout[i] <= 17);
        }

        
        
        struct timespec time1,time2;
        double dt2;

        clock_gettime(CLOCK_MONOTONIC,&time1);

        // building the tree index
        rand32_t* gen = rand32_init((uint32_t) time(NULL));
        
        assert(gen != NULL);
        /*for(int i=0;i<num_keys;i++)
            delimiter[i]=i;*/
        int32_t* delimiter = generate_sorted_unique(num_keys, gen);
        assert(delimiter != NULL);
        Tree* tree = build_index(num_levels, fanout, num_keys, delimiter);
        //free(delimiter);
        //free(fanout);
        if (tree == NULL) {
                //free(gen);
                exit(EXIT_FAILURE);
        }

        // generate probes
        int32_t* probe = generate(num_probes, gen);
        /*int32_t* probe =  (int32_t *)malloc(sizeof(int32_t) * num_probes);
        for(int i=0;i<num_probes;i++)
        {
        	probe[i] = (int32_t)i+25;
        }*/
        clock_gettime(CLOCK_MONOTONIC,&time2);
        uint64_t time_indexLoad = BILLION*(time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec);


        assert(probe != NULL);
        //free(gen);
        uint32_t* result_old_probe = malloc(sizeof(uint32_t) * num_probes);
        uint32_t* result_probe_generic = malloc(sizeof(uint32_t) * num_probes);
        uint32_t* result_probe_hardcoded = malloc(sizeof(uint32_t) * num_probes);
        //assert(result != NULL);

        // perform index probing (Phase 2)
        //old probe
        clock_gettime(CLOCK_MONOTONIC,&time1);
        for (size_t i = 0; i < num_probes; ++i) {
                    result_old_probe[i] = probe_index(tree, probe[i]);
        }
        clock_gettime(CLOCK_MONOTONIC,&time2);
        uint64_t time_probing = BILLION*(time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec);
        
        //new probe with Intel SSE Instructions
        clock_gettime(CLOCK_MONOTONIC,&time1);
        for (size_t i = 0; i < num_probes; ++i) {
                    result_probe_generic[i] = probe_index_generic(tree, probe[i]);
                    if(result_probe_generic[i]!=result_old_probe[i])
                        printf("Not Same\n");
        }
        clock_gettime(CLOCK_MONOTONIC,&time2);
        uint64_t time_generic_probing = BILLION*(time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec);

        //if tree 9-5-9 then hardcoded probe
        bool hardcoded=false;
        clock_gettime(CLOCK_MONOTONIC,&time1);
        //printf("Tree Details:%zu %zu %zu %zu",tree->num_levels,tree->node_capacity[0],tree->node_capacity)
        if(tree->num_levels==3)
        {
            if(tree->node_capacity[0]==8 && tree->node_capacity[1]==4 && tree->node_capacity[2]==8)
            {  
                result_probe_hardcoded=probe_index_959(tree,probe,num_probes);
                hardcoded=true;
            }
        }
        clock_gettime(CLOCK_MONOTONIC,&time2);
        uint64_t time_hardcoded_probing = BILLION*(time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec);

        
        // output results
        clock_gettime(CLOCK_MONOTONIC,&time1);
        for (size_t i = 0; i < num_probes; ++i) {
                fprintf(stdout, "%d %u\n", probe[i], result_probe_generic[i]);
        }
        clock_gettime(CLOCK_MONOTONIC,&time2);
        uint64_t time_output = BILLION*(time2.tv_sec - time1.tv_sec) + (time2.tv_nsec - time1.tv_nsec);

        printf("\n\nTime Phase 1(Creating index and initializing probe): %llu nanoseconds  \n", (long long unsigned int)time_indexLoad);
        printf("Time Phase 2(Performing probes - old way) : %llu seconds \n",(long long unsigned int)time_probing);
        printf("Time Phase 2(Performing probes - with Intel SSE Instructions) : %llu nanoseconds \n",(long long unsigned int)time_generic_probing);
        if(hardcoded==true)
            printf("Time Phase 2(Performing probes - hardcoded for 9 5 9 tree) : %llu nanoseconds \n",(long long unsigned int)time_hardcoded_probing);
        printf("Time Phase 3(Printing result to stdout) : %llu nanoseconds  \n\n\n",(long long unsigned int)time_output);
        

        // cleanup and exit
        //free(result_old_probe);
        //free(result_probe_generic);
        //free(result_probe_hardcoded);
        //free(probe);
        //cleanup_index(tree);
        return EXIT_SUCCESS;
}
