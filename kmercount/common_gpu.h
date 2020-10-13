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

const keyType max64 = 18446744073709551615;

const uint64_t kHashTableCapacity = 2 * 128 * 1024 * 1024;

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
    return k;// & (kHashTableCapacity-1);
}

inline __device__ keyType cuda_murmur3_64( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;// & (kHashTableCapacity-1);
}


inline __device__ keyType cuda_MurmurHash3_x64_128(const void* key, const uint32_t len, const uint32_t seed)
{
  const uint8_t * data = (const uint8_t*)key;
  const uint32_t nblocks = len / 16;
  int32_t i;

  uint64_t h1 = seed;
  uint64_t h2 = seed;

  uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

  //----------
  // body

  // const uint64_t * blocks = (const uint64_t *)(data);

  // // nblocks is zero for now
  // for(i = 0; i < nblocks; i++)
  // {
  //   uint64_t k1 = getblock(blocks,i*2+0);
  //   uint64_t k2 = getblock(blocks,i*2+1);

  //   int8_t r = 31;
  //   keytype x = k1;
  //   k1 *= c1; k1  = (x << r) | (x >> (64 - r)); //ROTL64(k1,31); 
  //   k1 *= c2; h1 ^= k1;

  //   r = 27;
  //   x= h1;
  //   h1 = (x << r) | (x >> (64 - r));//ROTL64(h1,27); 
  //   h1 += h2; h1 = h1*5+0x52dce729;

  //   r = 33;
  //   x = k2;
  //   k2 *= c2; k2  = (x << r) | (x >> (64 - r));//ROTL64(k2,33); 
  //   k2 *= c1; h2 ^= k2;

  //   r = 31;
  //   x = h2;
  //   h2 = (x << r) | (x >> (64 - r));//ROTL64(h2,31); 
  //   h2 += h1; h2 = h2*5+0x38495ab5;
  // }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  uint64_t k1 = 0;
  uint64_t k2 = 0;

  switch(len & 15)
  {
  case 15: k2 ^= (uint64_t)(tail[14]) << 48;
  case 14: k2 ^= (uint64_t)(tail[13]) << 40;
  case 13: k2 ^= (uint64_t)(tail[12]) << 32;
  case 12: k2 ^= (uint64_t)(tail[11]) << 24;
  case 11: k2 ^= (uint64_t)(tail[10]) << 16;
  case 10: k2 ^= (uint64_t)(tail[ 9]) << 8;
  case  9: k2 ^= (uint64_t)(tail[ 8]) << 0;
  k2 *= c2; k2  =  (k2 << 33) | (k2 >> (64 - 33)) ;//ROTL64(k2,33); 
  k2 *= c1; h2 ^= k2;

  case  8: k1 ^= (uint64_t)(tail[ 7]) << 56;
  case  7: k1 ^= (uint64_t)(tail[ 6]) << 48;
  case  6: k1 ^= (uint64_t)(tail[ 5]) << 40;
  case  5: k1 ^= (uint64_t)(tail[ 4]) << 32;
  case  4: k1 ^= (uint64_t)(tail[ 3]) << 24;
  case  3: k1 ^= (uint64_t)(tail[ 2]) << 16;
  case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;
  case  1: k1 ^= (uint64_t)(tail[ 0]) << 0;
    k1 *= c1; k1  = (k1 << 31) | (k1 >> (64 - 31)) ;//ROTL64(k1,31); 
    k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len; h2 ^= len;

  h1 += h2;
  h2 += h1;

  // h1 = fmix64(h1);

  // /* regular murmur64 */
  keyType k  = h1;
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  h1 = k;

  // h2 = fmix64(h2);

  k  = h2;
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  h2 = k;

  h1 += h2;
  h2 += h1;

  // ((uint64_t*)out)[0] = h1;
  // ((uint64_t*)out)[1] = h2;

  return h1;
}

#endif