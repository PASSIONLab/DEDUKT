// #pragma once

// #define keyType uint32_t

#define keyType unsigned long long int
// #define keyType unsigned int
struct KeyValue 
{
    keyType key;
    uint32_t value;
};

inline cudaError_t checkCuda(cudaError_t result, int s){

  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error in line : %s - %d\n", cudaGetErrorString(result), s);
    // assert(result == cudaSuccess);
  }
  return result;
}


inline void cuda_timer_start(cudaEvent_t start){
	checkCuda(cudaEventRecord(start), __LINE__);
}
inline void cuda_timer_stop(cudaEvent_t start, cudaEvent_t stop, float &mili){
	checkCuda(cudaEventRecord(stop), __LINE__);
    cudaEventSynchronize(stop);
    checkCuda(cudaEventElapsedTime(&mili, start, stop), __LINE__);
    cudaDeviceSynchronize();
}


const uint64_t kHashTableCapacity = 128 * 1024 * 1024;

const uint64_t kNumKeyValues = kHashTableCapacity / 2;//4168188;

const keyType kEmpty = 0;


KeyValue* create_hashtable_GPU(int rank);

void insert_hashtable(KeyValue* hashtable, keyType* kvs, uint32_t num_kvs, int rank);

void lookup_hashtable(KeyValue* hashtable, KeyValue* kvs, uint32_t num_kvs);

void delete_hashtable(KeyValue* hashtable, const KeyValue* kvs, uint32_t num_kvs);

std::vector<KeyValue> iterate_hashtable(KeyValue* hashtable);

void destroy_hashtable(KeyValue* hashtable, int rank);

uint64_t * getKmers_GPU( char *seq, int klen, int nproc, int *owner_counter, int rank);
