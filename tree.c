#include "tree.h"
#include <x86intrin.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
extern int posix_memalign(void** memptr, size_t alignment, size_t size);
size_t alignment = 16;

Tree* build_index(size_t num_levels, size_t fanout[], size_t num_keys, int32_t key[]) {
        // return null pointer for invalid tree configuration
        size_t min_num_keys = 1;
        for (size_t i = 0; i < num_levels - 1; ++i) {
                min_num_keys *= fanout[i];
        }
        size_t max_num_keys = min_num_keys * fanout[num_levels - 1] - 1;
        if (num_keys < min_num_keys || num_keys > max_num_keys) {
                fprintf(stderr, "Error: incorrect number of keys, min %zu, max %zu\n", min_num_keys, max_num_keys);
                return NULL;
        }

        // initialize the tree index
        Tree* tree = malloc(sizeof(Tree));
        assert(tree != NULL);
        tree->num_levels = num_levels;
        tree->node_capacity = malloc(sizeof(size_t) * num_levels);
        assert(tree->node_capacity != NULL);
        for (size_t i = 0; i < num_levels; ++i) {
                tree->node_capacity[i] = fanout[i] - 1;
        }
        tree->key_array = malloc(sizeof(int32_t*) * num_levels);
        assert(tree->key_array != NULL);
        size_t* key_count = malloc(sizeof(size_t) * num_levels);
        assert(key_count != NULL);
        size_t* array_capacity = malloc(sizeof(size_t) * num_levels);
        assert(array_capacity != NULL);
        for (size_t i = 0; i < num_levels; ++i) {
                size_t size = sizeof(int32_t) * tree->node_capacity[i];         // allocate one node per level
                int error = posix_memalign((void**) &(tree->key_array[i]), alignment, size);
                assert(error == 0);
                key_count[i] = 0;
                array_capacity[i] = tree->node_capacity[i];     // array_capacity[i] is always a multiple of node_capacity[i]
        }


        // insert sorted keys into index
        for (size_t i = 1; i < num_keys; ++i) {
                assert(key[i - 1] < key[i]);
        }
        for (size_t i = 0; i < num_keys; ++i) {
                size_t level = num_levels - 1;
                while (key_count[level] == array_capacity[level])
                        level -= 1;
                tree->key_array[level][key_count[level]] = key[i];
                key_count[level] += 1;
                while (level < num_levels - 1) {
                        level += 1;
                        size_t new_capacity = array_capacity[level] + tree->node_capacity[level];
                        size_t size = sizeof(int32_t) * new_capacity;           // allocate one more node
                        int32_t* new_array = NULL;
                        int error = posix_memalign((void**) &new_array, alignment, size);
                        assert(error == 0);
                        memcpy(new_array, tree->key_array[level], sizeof(int32_t) * key_count[level]);
                        free(tree->key_array[level]);
                        tree->key_array[level] = new_array;
                        array_capacity[level] = new_capacity;
                }
        }

        // pad with INT32_MAXs
        for (size_t i = 0; i < num_levels; ++i) {
                for (size_t j = key_count[i]; j < array_capacity[i]; ++j)
                        tree->key_array[i][j] = INT32_MAX;
                key_count[i] = array_capacity[i];
        }

        // print the tree
         for (size_t i = 0; i < num_levels; ++i) {
                 printf("Level %zu:", i);
                 for (size_t j = 0; j < key_count[i]; ++j)
                         printf(" %d", tree->key_array[i][j]);
                 printf("\n");
         }

        /*printf("number of levels is : %zu", num_levels);
        for (size_t i = 0; i < num_levels; ++i) {
                tree->node_capacity[i] = fanout[i] - 1;
                printf("node capacity for %zu level is : %zu\n", i,tree->node_capacity[i]);
                for(int ind = 0;ind < tree->node_capacity[i];ind++)
                {
                        printf("%zu ",tree->key_array[i][ind]);
                } 
                printf("\n");
        }*/

        //free(array_capacity);
        //free(key_count);
        return tree;
}
void print128_num_tree(__m128i var)
{
    uint32_t *val = (uint32_t*) &var;
    printf("Key: %i %i %i %i\n", 
           val[0],val[1],val[2],val[3]);
}
void printKey(__m128i var)
{
    uint32_t *val = (uint32_t*) &var;
    printf("Key: %i ", 
           val[0]);
}
uint32_t probe_index_generic(Tree* tree, int32_t probe_key) 
{
	//load and broadcast key to register
	__m128i k = _mm_cvtsi32_si128((int) probe_key);
	register __m128i key = _mm_shuffle_epi32(k , _MM_SHUFFLE(0,0,0,0));
	size_t result = 0;
    for (size_t level = 0; level < tree->num_levels; ++level) {
    		//printf("Level %zu ",level);
    		if(tree->node_capacity[level]==4)
    		{
    			__m128i del_ABCD = _mm_load_si128((__m128i*)&tree->key_array[level][result<<2]);
    			__m128i cmp = _mm_cmpgt_epi32(key,del_ABCD);
    			int32_t mask = _mm_movemask_ps((__m128)cmp);
    			int32_t r = _bit_scan_forward(mask ^ 0x1FF);
    			result = (result<<2)+result+r;
    		}
    		else if(tree->node_capacity[level]==8)
    		{
    			__m128i del_ABCD = _mm_load_si128((__m128i*)&tree->key_array[level][result<<3]);
        		__m128i del_EFGH = _mm_load_si128((__m128i*)&tree->key_array[level][(result<<3)+4]);
        		__m128i cmp_ABCD = _mm_cmpgt_epi32(key,del_ABCD);
                __m128i cmp_EFGH = _mm_cmpgt_epi32(key,del_EFGH);
        		__m128i cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
                cmp_A_to_H = _mm_packs_epi16(cmp_A_to_H, _mm_setzero_si128());
                int32_t mask = _mm_movemask_epi8(cmp_A_to_H);
                int32_t r = _bit_scan_forward((int32_t) (mask ^ 0x1FFFF));
                result= (result<<3)+result+r;

    		}
    		else if(tree->node_capacity[level]==16)
    		{
    			__m128i del_ABCD = _mm_load_si128((__m128i*)&tree->key_array[level][result<<4]);
        		__m128i del_EFGH = _mm_load_si128((__m128i*)&tree->key_array[level][(result<<4)+4]);
        		__m128i del_IJKL = _mm_load_si128((__m128i*)&tree->key_array[level][(result<<4)+8]);
        		__m128i del_MNOP = _mm_load_si128((__m128i*)&tree->key_array[level][(result<<4)+12]);
        		
        		__m128i cmp_ABCD = _mm_cmpgt_epi32(key,del_ABCD);
        		__m128i cmp_EFGH = _mm_cmpgt_epi32(key,del_EFGH);
        		__m128i cmp_IJKL = _mm_cmpgt_epi32(key,del_IJKL);
        		__m128i cmp_MNOP = _mm_cmpgt_epi32(key,del_MNOP);

        		__m128i cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
        		__m128i cmp_I_to_P = _mm_packs_epi32(cmp_IJKL, cmp_MNOP);
        		__m128i cmp_A_to_P = _mm_packs_epi32(cmp_A_to_H, cmp_I_to_P);
        		int32_t mask = _mm_movemask_epi8(cmp_A_to_P);
        		int32_t r = _bit_scan_forward((int32_t) (mask ^ 0x1FFFFFFFF));
        		result= (result<<4)+result+r;
    		}
    		
    }
    
    return (uint32_t) result;
}
uint32_t* probe_index_959(Tree* tree, int32_t* probe_key,size_t n) 
{
        //loaded root node
        __m128i del_ABCD = _mm_load_si128((__m128i*)tree->key_array[0]);
        tree->key_array[0] += 4;
        __m128i del_EFGH = _mm_load_si128((__m128i*)tree->key_array[0]);
        int num_probes = (int)(n);
        int i=0;
        uint32_t* result = malloc(sizeof(uint32_t) * num_probes);

        while(num_probes>=0)
        {
                num_probes = num_probes - 4;
                
                //load 4 keys and broadcast
                __m128i key = _mm_load_si128((__m128i*)probe_key);
                probe_key +=4;
                register __m128i k1 = _mm_shuffle_epi32(key , _MM_SHUFFLE(0,0,0,0));
                register __m128i k2 = _mm_shuffle_epi32(key , _MM_SHUFFLE(1,1,1,1));
                register __m128i k3 = _mm_shuffle_epi32(key , _MM_SHUFFLE(2,2,2,2));
                register __m128i k4 = _mm_shuffle_epi32(key , _MM_SHUFFLE(3,3,3,3));

                //Access level 0
                //Key 0 - root level - fanout 9  
         		 __m128i cmp_ABCD = _mm_cmpgt_epi32(k1,del_ABCD);
                __m128i cmp_EFGH = _mm_cmpgt_epi32(k1,del_EFGH);         
                __m128i cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
                cmp_A_to_H = _mm_packs_epi16(cmp_A_to_H, _mm_setzero_si128());
                int32_t msk_0 = _mm_movemask_epi8(cmp_A_to_H);
                int32_t r_00 = _bit_scan_forward((int32_t) (msk_0 ^ 0x1FFFF));
               
               //Key 1 -root level access
                cmp_ABCD = _mm_cmpgt_epi32(k2,del_ABCD);
                cmp_EFGH = _mm_cmpgt_epi32(k2,del_EFGH);          
                cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
                cmp_A_to_H = _mm_packs_epi16(cmp_A_to_H, _mm_setzero_si128());
                msk_0 = _mm_movemask_epi8(cmp_A_to_H);
                int32_t r_10 = _bit_scan_forward((int32_t) (msk_0 ^ 0x1FFFF));

                //Key 2 root level access
                cmp_ABCD = _mm_cmpgt_epi32(k3,del_ABCD);
                cmp_EFGH = _mm_cmpgt_epi32(k3,del_EFGH);          
                cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
                cmp_A_to_H = _mm_packs_epi16(cmp_A_to_H, _mm_setzero_si128());
                msk_0 = _mm_movemask_epi8(cmp_A_to_H);
                int32_t r_20 = _bit_scan_forward((int32_t) (msk_0 ^ 0x1FFFF));

                //Key 3 root level access
                cmp_ABCD = _mm_cmpgt_epi32(k4,del_ABCD);
                cmp_EFGH = _mm_cmpgt_epi32(k4,del_EFGH);          
                cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
                cmp_A_to_H = _mm_packs_epi16(cmp_A_to_H, _mm_setzero_si128());
                msk_0 = _mm_movemask_epi8(cmp_A_to_H);
                int32_t r_30 = _bit_scan_forward((int32_t) (msk_0 ^ 0x1FFFF));

                //Access level 1
                //Key 0 level 1
                __m128i lvl_1 = _mm_load_si128((__m128i*)&tree->key_array[1][r_00 << 2]);
                __m128i cmp_1 = _mm_cmpgt_epi32(k1,lvl_1);
                int32_t msk_1 = _mm_movemask_ps( (__m128)cmp_1);
                int32_t r_01 = _bit_scan_forward(msk_1 ^ 0x1FF);
                r_01 += (r_00 << 2) + r_00;

                //Key 1 level 1
                lvl_1 = _mm_load_si128((__m128i*)&tree->key_array[1][r_10 << 2]);
                cmp_1 = _mm_cmpgt_epi32(k2,lvl_1);
                msk_1 = _mm_movemask_ps( (__m128)cmp_1);
                int32_t r_11 = _bit_scan_forward(msk_1 ^ 0x1FF);
                r_11 += (r_10 << 2) + r_10;

                //Key 2 level 1
                lvl_1 = _mm_load_si128((__m128i*)&tree->key_array[1][r_20 << 2]);
                cmp_1 = _mm_cmpgt_epi32(k3,lvl_1);
                msk_1 = _mm_movemask_ps( (__m128)cmp_1);
                int32_t r_21 = _bit_scan_forward(msk_1 ^ 0x1FF);
                r_21 += (r_20 << 2) + r_20;

                //Key 3 level 1
                lvl_1 = _mm_load_si128((__m128i*)&tree->key_array[1][r_30 << 2]);
                cmp_1 = _mm_cmpgt_epi32(k4,lvl_1);
                msk_1 = _mm_movemask_ps( (__m128)cmp_1);
                int32_t r_31 = _bit_scan_forward(msk_1 ^ 0x1FF);
                r_31 += (r_30 << 2) + r_30;

                // access level 2 of the index (9-way)
                //Key 0 level 2
                __m128i lvl_2_A = _mm_load_si128((__m128i*)&tree->key_array[2][ r_01 << 3]);
                __m128i lvl_2_B = _mm_load_si128((__m128i*)&tree->key_array[2][(r_01 << 3) + 4]);
                __m128i cmp_2_A = _mm_cmpgt_epi32(k1,lvl_2_A);
                __m128i cmp_2_B = _mm_cmpgt_epi32(k1,lvl_2_B);
                __m128i cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
                cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
                int32_t msk_2 = _mm_movemask_epi8(cmp_2);
                int32_t r_02 = _bit_scan_forward(msk_2 ^ 0x1FFFF);
                r_02 += (r_01 << 3) + r_01;
                result[i++]=r_02;

                //Key 1 access level 2 of the index (9-way)
                lvl_2_A = _mm_load_si128((__m128i*)&tree->key_array[2][ r_11 << 3]);
                lvl_2_B = _mm_load_si128((__m128i*)&tree->key_array[2][(r_11 << 3) + 4]);
                cmp_2_A = _mm_cmpgt_epi32(k2,lvl_2_A);
                cmp_2_B = _mm_cmpgt_epi32(k2,lvl_2_B);
                cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
                cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
                msk_2 = _mm_movemask_epi8(cmp_2);
                int32_t r_12 = _bit_scan_forward(msk_2 ^ 0x1FFFF);
                r_12 += (r_11 << 3) + r_11;
                result[i++]=r_12;

                //Key 2- key2 level 2 of the index (9-way)
                lvl_2_A = _mm_load_si128((__m128i*)&tree->key_array[2][ r_21 << 3]);
                lvl_2_B = _mm_load_si128((__m128i*)&tree->key_array[2][(r_21 << 3) + 4]);
                cmp_2_A = _mm_cmpgt_epi32(k3,lvl_2_A);
                cmp_2_B = _mm_cmpgt_epi32(k3,lvl_2_B);
                cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
                cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
                msk_2 = _mm_movemask_epi8(cmp_2);
                int32_t r_22 = _bit_scan_forward(msk_2 ^ 0x1FFFF);
                r_22 += (r_21 << 3) + r_21;
                result[i++]=r_22;

                //Key 3- key3 access level 2 of the index (9-way)
                lvl_2_A = _mm_load_si128((__m128i*)&tree->key_array[2][ r_31 << 3]);
                lvl_2_B = _mm_load_si128((__m128i*)&tree->key_array[2][(r_31 << 3) + 4]);
                cmp_2_A = _mm_cmpgt_epi32(k4,lvl_2_A);
                cmp_2_B = _mm_cmpgt_epi32(k4,lvl_2_B);
                cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
                cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
                msk_2 = _mm_movemask_epi8(cmp_2);
                int32_t r_32 = _bit_scan_forward(msk_2 ^ 0x1FFFF);
                r_32 += (r_31 << 3) + r_31;
                result[i++]=r_32;
        }
        
        return result;
}

uint32_t probe_index(Tree* tree, int32_t probe_key) {
        size_t result = 0;
        for (size_t level = 0; level < tree->num_levels; ++level) {
                size_t offset = result * tree->node_capacity[level];
                size_t low = 0;
                size_t high = tree->node_capacity[level];
                //printf("Tree Node Capacity Level\n");
                //printf("%zu",level);
                //printf("\n");
                //printf("%zu",tree->node_capacity[level]);
                while (low != high) {
                        size_t mid = (low + high) / 2;
                        if (tree->key_array[level][mid + offset] < probe_key)
                                low = mid + 1;
                        else
                                high = mid;
                }
                size_t k = low;       // should go to child k
                result = result * (tree->node_capacity[level] + 1) + k;
                //printf("Result at level %zu:%zu \n",level,result);
        }
        return (uint32_t) result;
}
void cleanup_index(Tree* tree) {
        free(tree->node_capacity);
        for (size_t i = 0; i < tree->num_levels; ++i)
                free(tree->key_array[i]);
        free(tree->key_array);
        free(tree);
}
