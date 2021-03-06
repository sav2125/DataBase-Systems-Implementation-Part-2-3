#ifndef TREE_H_
#define TREE_H_

#include <stdint.h>
#include <stdlib.h>

typedef struct {
        size_t num_levels;
        size_t* node_capacity;
        int32_t** key_array;
} Tree;

Tree* build_index(size_t num_levels, size_t fanout[], size_t num_keys, int32_t key[]);
uint32_t probe_index(Tree* tree, int32_t probe_key);
uint32_t* probe_index_959(Tree* tree, int32_t* probe_key,size_t no_of_keys);
uint32_t probe_index_generic(Tree* tree, int32_t probe_key);
void cleanup_index(Tree* tree);

#endif
