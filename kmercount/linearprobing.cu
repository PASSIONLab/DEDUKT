#include "stdio.h"
#include "stdint.h"
#include "vector"
#include "linearprobing.h"
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

// #include "nvbio/basic/timer.h"
// #include "nvbio/basic/console.h"
// #include "nvbio/basic/bloom_filter.h"
// #include "nvbio/basic/numbers.h"
// #include "nvbio/basic/primitives.h"
// #include "nvbio/basic/vector.h"
// #include "nvbio/basic/cuda/ldg.h"
// #include "nvbio/basic/omp.h"
#include <stdio.h>
#include <stdlib.h>

#define BIG_CONSTANT(x) (x)

cudaError_t checkCuda(cudaError_t result, int s){

  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error in line : %s - %d\n", cudaGetErrorString(result), s);
    // assert(result == cudaSuccess);
  }
  return result;
}

// 32 bit Murmur3 hash
__device__ keyType murmur3_32(keyType k) //TODO:: for now just changed keytype
{
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    k *= 0xc2b2ae35;
    k ^= k >> 16;
    return k & (kHashTableCapacity-1);
}

__device__ keyType murmur3_64( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k & (kHashTableCapacity-1);
}

// Hash funs and variables for NVBIO Bloom Filter 
static const uint32_t FILTER_K = 7;
static const uint32_t ITEMS_PER_THREAD = 100;
/*
//To avoid conflicts of Hash functions between diBella and nvbio
namespace BF_hash{

// IN:: copied from numbers.h in nvbio library
// compute a simple 32-bit hash

NVBIO_FORCEINLINE NVBIO_HOST_DEVICE 
uint32_t hash(uint32_t a)
{
    a = (a+0x7ed55d16) + (a<<12);
    a = (a^0xc761c23c) ^ (a>>19);
    a = (a+0x165667b1) + (a<<5);
    a = (a+0xd3a2646c) ^ (a<<9);
    a = (a+0xfd7046c5) + (a<<3);
    a = (a^0xb55a4f09) ^ (a>>16);
    return a;
}

/// Thomas Wang's 32 bit Mix Function: http://www.cris.com/~Ttwang/tech/inthash.htm
///
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE 
uint32_t hash2(uint32_t key)
{
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}

/// Thomas Wang's 64 bit Mix Function: http://www.cris.com/~Ttwang/tech/inthash.htm
///
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE 
uint64_t hash(uint64_t key)
{
    key += ~(key << 32);
    key ^=  (key >> 22);
    key += ~(key << 13);
    key ^=  (key >> 8);
    key +=  (key << 3);
    key ^=  (key >> 15);
    key += ~(key << 27);
    key ^=  (key >> 31);
    return key;
}

/// simple 64-bit hash
///
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE 
uint64_t hash2(uint64_t key)
{
    return (key >> 32) ^ key;
}

/// elf 64-bit hash
///
NVBIO_FORCEINLINE NVBIO_HOST_DEVICE 
uint64_t hash3(uint64_t key)
{
    uint32_t hash = 0u;

    #if defined(__CUDA_ARCH__)
    #pragma unroll
    #endif
    for (uint32_t i = 0; i < 8; ++i)
    {
        hash = (hash << 4) + ((key >> (i*8)) & 255u); // shift/mix

        // get high nybble
        const uint32_t hi_bits = hash & 0xF0000000;
        if (hi_bits != 0u)
            hash ^= hi_bits >> 24; // xor high nybble with second nybble

        hash &= ~hi_bits; // clear high nybble
    }
    return hash;
}
}//namespace BF_hash


struct hash_functor1
{
    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
    uint64_t operator() (const uint64_t kmer) const { return BF_hash::hash( kmer ); }
};
struct hash_functor2
{
    NVBIO_FORCEINLINE NVBIO_HOST_DEVICE
    uint64_t operator() (const uint64_t kmer) const { return BF_hash::hash2( kmer ); }
};
*/
// Create a hash table. For linear probing, this is just an array of KeyValues
KeyValue* create_hashtable_GPU(int rank) 
{
    int count, devId;
    cudaGetDeviceCount(&count);
    int gpuID = rank % count;
    cudaSetDevice(gpuID);
    // Allocate memory
    KeyValue* hashtable;
    cudaMalloc(&hashtable, sizeof(KeyValue) * kHashTableCapacity);

    // Initialize hash table to empty
    // static_assert(kEmpty == 0xffffffff, "memset expected kEmpty=0xffffffff");
    // static_assert(kEmpty == 0x1ffffffffu, "memset expected kEmpty=0xffffffff");
    // cudaMemset(hashtable, 0x1ffffffffu, sizeof(KeyValue) * kHashTableCapacity);
    cudaMemset(hashtable,  0, sizeof(KeyValue) * kHashTableCapacity);

    return hashtable;
}

// typedef nvbio::bloom_filter<2,hash_functor1,hash_functor2,uint32_t*> bloom_filter_type;

// Insert the key/values in kvs into the hashtable
__global__ void gpu_hashtable_insert(KeyValue* hashtable, 
    const keyType* kvs, 
    unsigned int numkvs)
{
    unsigned int threadid = blockIdx.x*blockDim.x + threadIdx.x;
   
    if (threadid < numkvs)
    {
        keyType new_key = kvs[threadid];//.key;
        keyType slot = murmur3_32(new_key);
        // if(filter.has(new_key))
        {

            while (true)
            {
                keyType old_key = atomicCAS(&hashtable[slot].key, kEmpty, new_key);
                      
                if (old_key == kEmpty || old_key == new_key)
                {
                    atomicAdd(&hashtable[slot].value,1);
                    return;
                }

                slot = (slot + 1) & (kHashTableCapacity-1);
            }
        }
    }
}

//TODO:: make it templated
 


// populate a nvbio bloom filter in parallel
// template <typename bloom_filter_type, typename T>
// __global__ void populate_kernel(const uint32_t N, const uint32_t* keys, bloom_filter_type filter)
// {
//     const uint32_t i = threadIdx.x + blockIdx.x * blockDim.x;
//     if (i < N)
//     {
//         filter.insert( keys[i] );
//     }
// }

std::vector<KeyValue> insert_hashtable(KeyValue* pHashTable, const keyType* keys, uint32_t num_keys, int rank)
{
    // Map MPI ranks to GPUs
    int count, devId;
    cudaGetDeviceCount(&count);
    int gpuID = rank % count;
    cudaSetDevice(gpuID);
    cudaGetDevice(&devId);
    // printf("\n FROnProcs %d: rank %d mapped to %d\n", nproc, rank, devId);

    // Copy the keyvalues to the GPU
    // Create events for GPU timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    /*------------------------
     Copy kmers to GPU      
    ------------------------*/
    keyType* device_keys;
    cudaEventRecord(start);
    
    cudaMalloc(&device_keys, sizeof(keyType) * num_keys); 
    cudaMemcpy(device_keys, keys, sizeof(keyType) * num_keys, cudaMemcpyHostToDevice);
    
    // cudaEventRecord(stop);
    // cudaEventSynchronize(stop);
    // float milliseconds = 0;
    // cudaEventElapsedTime(&milliseconds, start, stop);
    // float seconds = milliseconds / 1000.0f;
    // printf("    GPU memcopy: %d items in %f ms (%f million keys/second)\n", 
    //     num_keys, milliseconds, num_keys / (double)seconds / 1000000.0f);

    /*----------------------------
    Build NVBIO Bloom Filter       
    ------------------------------*/
    // build a set of 1M random integers
    const uint32_t N = 1;//1000000;
    // nvbio::vector<nvbio::host_tag,uint32_t> h_vector( N );
    // fill it up
    // for (uint32_t i = 0; i < N; ++i)
    //     h_vector[i] = rand();
    // copy it to the device
    // nvbio::vector<nvbio::device_tag,uint32_t> d_vector = h_vector;
    // construct an empty Bloom filter

    const uint32_t filter_words = N;  // let's use 32-bits per input item;
                                    // NOTE: this is still a lot less than 1-bit for each
                                    // of the 4+ billion integers out there...
    // nvbio::vector<nvbio::device_tag,uint32_t> d_filter_storage( filter_words, 0u );
    // bloom_filter_type d_filter( filter_words * 32, plain_view( d_filter_storage ) );

    int b = 128;
    int g= (N + (b - 1)) / b;
    // populate_kernel<<<g,128>>>( N, device_keys, d_filter );
   
    /*----------------------------
    CUDA call: Insert kmers to HT       
    ------------------------------*/

    // Have CUDA calculate the thread block size
    int mingridsize;
    int threadblocksize;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_hashtable_insert, 0, 0);

    // cudaEventRecord(start);

    // Insert all the keys into the hash table
    int gridsize = ((uint32_t)num_keys + threadblocksize - 1) / threadblocksize;
    gpu_hashtable_insert<<<gridsize, threadblocksize>>>(pHashTable, device_keys, (uint32_t)num_keys);

    cudaEventRecord(stop);

    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    float seconds = milliseconds / 1000.0f;
    printf("    GPU inserted %d items in %f ms (%f million keys/second)\n", 
        num_keys, milliseconds, num_keys / (double)seconds / 1000000.0f);

    std::vector<KeyValue> keys1;
    keys1.resize(kHashTableCapacity);

    cudaMemcpy(keys1.data(), pHashTable, sizeof(KeyValue) * kHashTableCapacity, cudaMemcpyDeviceToHost); 

    cudaFree(device_keys);
    return keys1;

}


__global__ void gpu_parse_kmer(char *seq, char *kmers, int klen, unsigned int seq_len){
    
    unsigned int tId = threadIdx.x;
    unsigned int laneId = tId & 31;
    unsigned int gId = (blockIdx.x * blockDim.x + tId);
    unsigned int wId = (gId % 128 ) >> 5;

    int nWarp = blockDim.x / 32;
    int per_block_seq_len = (seq_len + (gridDim.x - 1)) / gridDim.x;
    int per_warp_seq_len = (per_block_seq_len + (nWarp - 1)) / nWarp;
    int st_char_block = blockIdx.x * per_block_seq_len; //first char this block should read
    int st_char_warp = st_char_block + wId * per_warp_seq_len;
    int end_char = seq_len - klen + 1;

   for(int i = st_char_warp ; i < st_char_warp + per_warp_seq_len ; i++) {
        int kk = i;
        if(i <= end_char && laneId < klen)
            kmers[i * klen + laneId] = seq[kk + laneId]; // for (int j = 0; j < KMER_LENGTH; ++j,++kk)                
    }
}

__global__ void gpu_parseKmerNFillupBuff(char *seq, char *kmers, int klen, unsigned int seq_len,
    uint64_t* outgoing, int *owner_counter, int nproc, int p_buff_len){
    
    unsigned int tId = threadIdx.x;
    unsigned int laneId = tId & 31;
    unsigned int gId = (blockIdx.x * blockDim.x + tId);
    unsigned int wId = (gId % 128 ) >> 5;

    int nWarp = blockDim.x / 32;
    int per_block_seq_len = (seq_len + (gridDim.x - 1)) / gridDim.x;
    int per_warp_seq_len = (per_block_seq_len + (nWarp - 1)) / nWarp;
    int st_char_block = blockIdx.x * per_block_seq_len; //first char this block should read
    int st_char_warp = st_char_block + wId * per_warp_seq_len;
    int nKmer = seq_len - klen + 1;
    int warpSize = 32;

   for(int i = st_char_warp ; i < st_char_warp + per_warp_seq_len && (i + laneId) <= nKmer ; i+=warpSize) {
        keyType longs = 0; //GPU CAS support this for 64 bit
        bool endOfRead = false;
        for (int k = 0; k < klen; ++k) {
            char s =  seq[i + k + laneId];
            if(s == 'a') { endOfRead = true; break;}
            int j = k % 32;
            int l = k/32;
            // assert(s != '\0');
            size_t x = ((s) & 4) >> 1;
            longs |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer
        }
        if(endOfRead) continue;
        keyType owner = murmur3_64(longs) & (nproc - 1); // remove & with HTcapacity in func
        int old_count = atomicAdd(&owner_counter[owner],1); 

        if(old_count >= nKmer) return;
        if(old_count >= p_buff_len * 2 ) {
            printf("Overflow!! MISSION ABORT!!\n");
        }
        outgoing[owner * p_buff_len + old_count]=longs; //hash (longs)      
    }
}

uint64_t * getKmers_GPU(char *seq, int klen, int nproc, int *owner_counter, int rank){

    // printf("CHANGE\n");
    int count, devId;
        // Map MPI ranks to GPUs

    cudaGetDeviceCount(&count);
    int gpuID = rank % count;
    cudaSetDevice(gpuID);
    cudaGetDevice(&devId);
    printf("\n nProcs %d: rank %d mapped to %d\n", nproc, rank, devId);

    // nproc = 4;
    unsigned int seq_len = strlen(seq);
    unsigned int n_kmers =  seq_len - klen + 1;
    printf("nKMER with gaps %d\n", n_kmers );
    // printf("nkmers %u \n", n_kmers );

    char *d_kmers, *d_seq;
    uint64_t *d_outgoing, *d_outOverflowBuff;
    int *d_owner_counter; 

    // Create events for GPU timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    // checkCuda (cudaMalloc(&d_kmers, n_kmers * klen * sizeof(char*)), __LINE__);
    checkCuda (cudaMalloc(&d_outgoing, n_kmers*2 * sizeof(uint64_t*)), __LINE__);  // giving 2x space to each node 
    checkCuda (cudaMalloc(&d_seq, seq_len * sizeof(char*)), __LINE__);
    checkCuda (cudaMemcpy(d_seq, seq, seq_len * sizeof(char*) , cudaMemcpyHostToDevice), __LINE__);
    
    checkCuda (cudaMalloc(&d_owner_counter, nproc * sizeof(int)), __LINE__);
    cudaMemset(d_owner_counter,  0, sizeof(int) * nproc);

    int p_buff_len = ((n_kmers * 2) + nproc - 1)/nproc;
    int b = 128;
    int g = (seq_len + (b - 1)) / b;
    int per_block_seq_len = (seq_len + (g - 1)) / g;
    // gpu_parse_kmer<<<g, b>>>(d_seq, d_kmers, klen, seq_len);

    gpu_parseKmerNFillupBuff<<<g, b>>>(d_seq, d_kmers, klen, seq_len, d_outgoing, d_owner_counter, nproc, p_buff_len);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    float seconds = milliseconds / 1000.0f;
    // printf("    GPU parseNPack: n %f ms (%f million kmers/second)\n", 
    //      milliseconds, n_kmers / (double)seconds / 1000000.0f);


    // // char *h_kmers = (char *) malloc ( n_kmers * klen * sizeof(char*));
    // uint64_t *h_outgoing = (uint64_t *) malloc ( n_kmers * 2 * sizeof(uint64_t));

    // checkCuda (cudaMemcpy(h_outgoing, d_outgoing, n_kmers * 2 * sizeof(uint64_t) , cudaMemcpyDeviceToHost), __LINE__); 
    checkCuda (cudaMemcpy(owner_counter, d_owner_counter, nproc * sizeof(int) , cudaMemcpyDeviceToHost), __LINE__); 
   
    uint64_t total_counter = 0;
    printf("\nFrom GPU: "); 
    for (int i = 0; i < nproc; ++i)
    {
        total_counter += owner_counter[i];
       printf(" %d", owner_counter[i]);
    }
    printf(", Total: %d \n ", total_counter);
    // cudaFree(d_kmers);
    // cudaFree(d_outgoing);
    cudaFree(d_seq);
    cudaFree(d_owner_counter);
    return d_outgoing;
}


// Lookup keys in the hashtable, and return the values
__global__ void gpu_hashtable_lookup(KeyValue* hashtable, KeyValue* kvs, unsigned int numkvs)
{
    unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadid < kHashTableCapacity)
    {
        uint32_t key = kvs[threadid].key;
        uint32_t slot = murmur3_32(key);

        while (true)
        {
            if (hashtable[slot].key == key)
            {
                kvs[threadid].value = hashtable[slot].value;
                return;
            }
            if (hashtable[slot].key == kEmpty)
            {
                kvs[threadid].value = kEmpty;
                return;
            }
            slot = (slot + 1) & (kHashTableCapacity - 1);
        }
    }
}

void lookup_hashtable(KeyValue* pHashTable, KeyValue* kvs, uint32_t num_kvs)
{
    // Copy the keyvalues to the GPU
    KeyValue* device_kvs;
    cudaMalloc(&device_kvs, sizeof(KeyValue) * num_kvs);
    cudaMemcpy(device_kvs, kvs, sizeof(KeyValue) * num_kvs, cudaMemcpyHostToDevice);

    // Have CUDA calculate the thread block size
    int mingridsize;
    int threadblocksize;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_hashtable_insert, 0, 0);

    // Create events for GPU timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    // Insert all the keys into the hash table
    int gridsize = ((uint32_t)num_kvs + threadblocksize - 1) / threadblocksize;
    // gpu_hashtable_insert << <gridsize, threadblocksize >> > (pHashTable, device_kvs, (uint32_t)num_kvs);

    cudaEventRecord(stop);

    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    float seconds = milliseconds / 1000.0f;
    printf("    GPU lookup %d items in %f ms (%f million keys/second)\n",
        num_kvs, milliseconds, num_kvs / (double)seconds / 1000000.0f);

    cudaFree(device_kvs);
}

// Delete each key in kvs from the hash table, if the key exists
// A deleted key is left in the hash table, but its value is set to kEmpty
// Deleted keys are not reused; once a key is assigned a slot, it never moves
__global__ void gpu_hashtable_delete(KeyValue* hashtable, const KeyValue* kvs, unsigned int numkvs)
{
    unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadid < kHashTableCapacity)
    {
        uint32_t key = kvs[threadid].key;
        uint32_t slot = murmur3_32(key);

        while (true)
        {
            if (hashtable[slot].key == key)
            {
                hashtable[slot].value = kEmpty;
                return;
            }
            if (hashtable[slot].key == kEmpty)
            {
                return;
            }
            slot = (slot + 1) & (kHashTableCapacity - 1);
        }
    }
}

void delete_hashtable(KeyValue* pHashTable, const KeyValue* kvs, uint32_t num_kvs)
{
    // Copy the keyvalues to the GPU
    KeyValue* device_kvs;
    cudaMalloc(&device_kvs, sizeof(KeyValue) * num_kvs);
    cudaMemcpy(device_kvs, kvs, sizeof(KeyValue) * num_kvs, cudaMemcpyHostToDevice);

    // Have CUDA calculate the thread block size
    int mingridsize;
    int threadblocksize;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_hashtable_insert, 0, 0);

    // Create events for GPU timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);

    // Insert all the keys into the hash table
    int gridsize = ((uint32_t)num_kvs + threadblocksize - 1) / threadblocksize;
    gpu_hashtable_delete<< <gridsize, threadblocksize >> > (pHashTable, device_kvs, (uint32_t)num_kvs);

    cudaEventRecord(stop);

    cudaEventSynchronize(stop);

    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    float seconds = milliseconds / 1000.0f;
    printf("    GPU delete %d items in %f ms (%f million keys/second)\n",
        num_kvs, milliseconds, num_kvs / (double)seconds / 1000000.0f);

    cudaFree(device_kvs);
}

// Iterate over every item in the hashtable; return non-empty key/values
__global__ void gpu_iterate_hashtable(KeyValue* pHashTable, KeyValue* kvs, uint32_t* kvs_size)
{
    unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadid < kHashTableCapacity) 
    {
        if (pHashTable[threadid].key != kEmpty) 
        {
            uint32_t value = pHashTable[threadid].value;
            if (value != kEmpty)
            {
                uint32_t size = atomicAdd(kvs_size, 1);
                kvs[size] = pHashTable[threadid];
            }
        }
    }
}

std::vector<KeyValue> iterate_hashtable(KeyValue* pHashTable)
{
    uint32_t* device_num_kvs;
    cudaMalloc(&device_num_kvs, sizeof(uint32_t));
    cudaMemset(device_num_kvs, 0, sizeof(uint32_t));

    KeyValue* device_kvs;
    cudaMalloc(&device_kvs, sizeof(KeyValue) * kNumKeyValues);

    int mingridsize;
    int threadblocksize;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_iterate_hashtable, 0, 0);

    int gridsize = (kHashTableCapacity + threadblocksize - 1) / threadblocksize;

    gpu_iterate_hashtable<<<gridsize, threadblocksize>>>(pHashTable, device_kvs, device_num_kvs);

    uint32_t num_kvs;
    cudaMemcpy(&num_kvs, device_num_kvs, sizeof(uint32_t), cudaMemcpyDeviceToHost);

    std::vector<KeyValue> kvs;
    kvs.resize(num_kvs);

    cudaMemcpy(kvs.data(), device_kvs, sizeof(KeyValue) * num_kvs, cudaMemcpyDeviceToHost);

    cudaFree(device_kvs);
    cudaFree(device_num_kvs);

    return kvs;
}



// Free the memory of the hashtable
void destroy_hashtable(KeyValue* pHashTable, int rank)
{
    cudaFree(pHashTable);
}