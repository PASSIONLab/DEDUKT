// #define keyType uint32_t

#ifndef COMMON_GPU_H
#define COMMON_GPU_H
#include <cuda_runtime_api.h>
#include <cuda.h>

typedef unsigned long long int uint64_cu; 
typedef unsigned long long int keyType; 

struct KeyValue 
{
    keyType key;
    uint32_t value;
};

const keyType kEmpty = 0;

const uint64_t kHashTableCapacity = 128 * 1024 * 1024;

const uint64_t kNumKeyValues = kHashTableCapacity / 2;//4168188;

inline cudaError_t checkCuda(cudaError_t result, int s){
	if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error in line : %s - %d\n", cudaGetErrorString(result), s);
    // assert(result == cudaSuccess);
  }
  return result;
}

inline void cuda_timer_start(cudaEvent_t start){
	cudaEventRecord(start);
}
inline void cuda_timer_stop(cudaEvent_t start, cudaEvent_t stop, float &mili){
	cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&mili, start, stop);
    cudaDeviceSynchronize();
}

#define BIG_CONSTANT(x) (x)

// 32 bit Murmur3 hash
inline __device__ keyType cuda_murmur3_32(keyType k) //TODO:: for now just changed keytype
{
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    k *= 0xc2b2ae35;
    k ^= k >> 16;
    return k & (kHashTableCapacity-1);
}

inline __device__ keyType cuda_murmur3_64( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k & (kHashTableCapacity-1);
}


#endif