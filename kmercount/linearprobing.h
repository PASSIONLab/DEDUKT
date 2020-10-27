
#ifndef LINEARPROBING_H
#define LINEARPROBING_H


#include "common_gpu.h"
// #define keyType unsigned long long int
// #define keyType unsigned int
// struct KeyValue 
// {
//     keyType key;
//     uint32_t value;
// };


// inline cudaError_t checkCuda(cudaError_t result, int s){
//     if (result != cudaSuccess) {
//     fprintf(stderr, "CUDA Runtime Error in line : %s - %d\n", cudaGetErrorString(result), s);
//     // assert(result == cudaSuccess);
//   }
//   return result;
// }

KeyValue* create_hashtable_GPU(int rank);

void insert_hashtable(KeyValue* hashtable, keyType* kvs, uint32_t num_kvs, int rank);

void lookup_hashtable(KeyValue* hashtable, KeyValue* kvs, uint32_t num_kvs);

void delete_hashtable(KeyValue* hashtable, const KeyValue* kvs, uint32_t num_kvs);

std::vector<KeyValue> iterate_hashtable(KeyValue* hashtable);

void destroy_hashtable(KeyValue* hashtable, int rank);

uint64_t * getKmers_GPU( char *seq, int klen, int nproc, int *owner_counter, int rank, int BUFF_LEN);
#endif
