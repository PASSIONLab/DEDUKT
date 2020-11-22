#include "stdio.h"
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include "stdint.h"
#include "vector"
#include "linearprobing.h"
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
using namespace std;
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

	if (threadid < numkvs){
		keyType new_key = kvs[threadid];//.key;
		keyType slot = cuda_murmur3_64(new_key) & (kHashTableCapacity-1);

		while (true){
			keyType old_key = atomicCAS(&hashtable[slot].key, kEmpty, new_key);

			if (old_key == kEmpty || old_key == new_key){
				atomicAdd(&hashtable[slot].value,1);
				return;
			}
			slot = (slot + 1) & (kHashTableCapacity-1);
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

// void insert_hashtable(KeyValue* pHashTable, const keyType* keys, uint32_t num_keys, int rank)

void insert_hashtable(KeyValue* pHashTable, vector <keyType>& recvbuf, uint32_t nkmers, int rank, int nprocs, int *recvcnt, int p_buff_len)
{

	// Map MPI ranks to GPUs
	int count, devId;
	cudaGetDeviceCount(&count);
	int gpuID = rank % count;
	cudaSetDevice(gpuID);
	keyType *d_kmers;
	// int p_buff_len = ((nkmers * BUFF_SCALE) + nprocs - 1)/nprocs;

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start);
	size_t recvdKmers = 0;
	
	cudaMalloc(&d_kmers, sizeof(keyType) * nkmers); 

	for(uint64_t i= 0; i < nprocs ; ++i) {
		cudaMemcpy(d_kmers + recvdKmers, &recvbuf[i * p_buff_len], sizeof(keyType) * recvcnt[i], cudaMemcpyHostToDevice); 
		recvdKmers += recvcnt[i];	
	}

	int b = 128;
	int g = (b - 1) / b;

	int mingridsize;
	int threadblocksize;
	cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_hashtable_insert, 0, 0);

	int gridsize = ((uint32_t)nkmers + threadblocksize - 1) / threadblocksize;
	
	gpu_hashtable_insert<<<gridsize, threadblocksize>>>(pHashTable, d_kmers, (uint32_t)nkmers);

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);

	cudaFree(d_kmers);
	
	return ;
}

// void insert_bloomfilter(KeyValue* pHashTable, keyType* d_kmers, uint32_t nkmers, int rank)
// {

// 	// Map MPI ranks to GPUs
// 	int count, devId;
// 	cudaGetDeviceCount(&count);
// 	int gpuID = rank % count;
// 	cudaSetDevice(gpuID);
// 	cudaGetDevice(&devId);
// 	// printf("\n FROnProcs %d: rank %d mapped to %d\n", nproc, rank, devId);

// 	// Copy the keyvalues to the GPU
// 	// Create events for GPU timing
// 	cudaEvent_t start, stop;
// 	cudaEventCreate(&start);
// 	cudaEventCreate(&stop);

// 	/*------------------------
// 	  Copy kmers to GPU      
// 	  ------------------------*/
// 	// keyType* device_keys;
// 	cudaEventRecord(start);
// 	size_t recvdKmers = 0;
// 	cudaMalloc(&d_kmers, sizeof(keyType) * totrecv); 

// 	for(uint64_t i= 0; i < nprocs ; ++i) {
// 		cudaMemcpy(d_kmers + recvdKmers, &recvbuf[i * p_buff_len], sizeof(keyType) * recvcnt[i], cudaMemcpyHostToDevice); 
// 		recvdKmers += recvcnt[i];	
// 	}
// 	// checkCuda( cudaMalloc(&device_keys, sizeof(keyType) * num_keys), __LINE__); 
// 	// checkCuda( cudaMemcpy(device_keys, keys, sizeof(keyType) * num_keys, cudaMemcpyHostToDevice), __LINE__); 

// 	// cudaEventRecord(stop);
// 	// cudaEventSynchronize(stop);
// 	// float milliseconds = 0;
// 	// cudaEventElapsedTime(&milliseconds, start, stop);
// 	// float seconds = milliseconds / 1000.0f;
// 	// printf("    GPU memcopy: %d items in %f ms (%f million keys/second)\n", 
// 	//     num_keys, milliseconds, num_keys / (double)seconds / 1000000.0f);

// 	/*----------------------------
// 	  Build NVBIO Bloom Filter       
// 	  ------------------------------*/
// 	// build a set of 1M random integers
// 	const uint32_t N = 1;//1000000;
// 	// nvbio::vector<nvbio::host_tag,uint32_t> h_vector( N );
// 	// fill it up
// 	// for (uint32_t i = 0; i < N; ++i)
// 	//     h_vector[i] = rand();
// 	// copy it to the device
// 	// nvbio::vector<nvbio::device_tag,uint32_t> d_vector = h_vector;
// 	// construct an empty Bloom filter

// 	const uint32_t filter_words = N;  // let's use 32-bits per input item;
// 	// NOTE: this is still a lot less than 1-bit for each
// 	// of the 4+ billion integers out there...
// 	// nvbio::vector<nvbio::device_tag,uint32_t> d_filter_storage( filter_words, 0u );
// 	// bloom_filter_type d_filter( filter_words * 32, plain_view( d_filter_storage ) );

// 	int b = 128;
// 	int g= (N + (b - 1)) / b;
// 	// populate_kernel<<<g,128>>>( N, device_keys, d_filter );

// 	/*----------------------------
// 	  CUDA call: Insert kmers to HT       
// 	  ------------------------------*/

// 	// Have CUDA calculate the thread block size
// 	int mingridsize;
// 	int threadblocksize;
// 	cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, gpu_hashtable_insert, 0, 0);

// 	// cudaEventRecord(start);

// 	// Insert all the keys into the hash table

// 	int gridsize = ((uint32_t)nkmers + threadblocksize - 1) / threadblocksize;
// 	gpu_hashtable_insert<<<gridsize, threadblocksize>>>(pHashTable, d_kmers, (uint32_t)nkmers);

// 	cudaEventRecord(stop);

// 	cudaEventSynchronize(stop);

// 	// float milliseconds = 0;
// 	// cudaEventElapsedTime(&milliseconds, start, stop);
// 	// float seconds = milliseconds / 1000.0f;
// 	// printf("    GPU inserted %d items in %f ms (%f million keys/second)\n", 
// 	//     num_keys, milliseconds, num_keys / (double)seconds / 1000000.0f);

// 	// std::vector<KeyValue> h_pHashTable;
// 	// h_pHashTable.resize(kHashTableCapacity);

// 	// cudaMemcpy(h_pHashTable.data(), pHashTable, sizeof(KeyValue) * kHashTableCapacity, cudaMemcpyDeviceToHost); 

// 	cudaFree(device_keys);
// 	return ;//h_pHashTable;

// }


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
	unsigned int laneId = tId & (blockDim.x - 1);
	unsigned int gId = (blockIdx.x * blockDim.x + tId);
	int per_block_seq_len = blockDim.x;//(seq_len + (gridDim.x - 1)) / gridDim.x;
	int st_char_block = blockIdx.x * per_block_seq_len; //first char this block should read
	int nKmer = seq_len - klen + 1; //last char is 'a'

	//*parse and compress kmer*
	for(int i = st_char_block + laneId; i  < (st_char_block + per_block_seq_len) && i < nKmer ; i+=blockDim.x) {
		keyType longs = 0; //GPU CAS support this for 64 bit
		bool validKmer = true;

		for (int k = 0; k < klen; ++k) {
			char s =  seq[i + k ];
			if(s == 'a' || s == 'N')  {
				validKmer = false; break;
			}
			int j = k % 32;
			size_t x = ((s) & 4) >> 1;
			longs |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer
		}
		
		//*populate outgoing buffer*
		if(validKmer ) {
			// keyType myhash = cuda_murmur3_64(longs); // remove & with HTcapacity in func
			keyType myhash = cuda_MurmurHash3_x64_128((const void *)&longs, 8, 313);// & (nproc - 1);
			double range = static_cast<double>(myhash) * static_cast<double>(nproc);
			size_t owner = range / max64; // dest proc rank

			int old_count = atomicAdd(&owner_counter[owner],1); 

			if(old_count > p_buff_len) {
				printf("Overflow!! MISSION ABORT!!\n"); return;
			}
			outgoing[owner * p_buff_len + old_count]=longs; //hash (longs)   
		}   
	}
}

void getKmers_GPU(string seq, vector<uint64_t> & h_outgoing, int klen, int nproc, vector<int> & owner_counter, int rank, int BUFF_SCALE){

	int count, devId;
	char *d_kmers, *d_seq;
	uint64_t *d_outOverflowBuff, *d_outgoing;
	int *d_owner_counter; 

	// Map MPI ranks to GPUs
	cudaGetDeviceCount(&count);
	int gpuID = rank % count;
	cudaSetDevice(gpuID);
	cudaGetDevice(&devId);

	unsigned int seq_len = seq.length();
	if(seq_len < klen) return;// h_outgoing;
	unsigned int n_kmers =  seq_len - klen + 1;

	//* Create events for GPU timing*
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float milliseconds = 0;

	cuda_timer_start(start);

	//* CUDA mallocs *
	checkCuda (cudaMalloc(&d_outgoing, n_kmers * BUFF_SCALE * sizeof(uint64_t)), __LINE__);  // giving 2x space to each node 
	checkCuda (cudaMalloc(&d_seq, seq_len * sizeof(char)), __LINE__);
	checkCuda (cudaMalloc(&d_owner_counter, nproc * sizeof(int)), __LINE__);

	checkCuda (cudaMemcpy(d_seq, &seq[0], seq_len * sizeof(char) , cudaMemcpyHostToDevice), __LINE__);
	cudaMemset(d_outgoing,  0, n_kmers * BUFF_SCALE * sizeof(uint64_t));
	cudaMemset(d_owner_counter,  0, sizeof(int) * nproc);

	int p_buff_len = ((n_kmers * BUFF_SCALE) + nproc - 1)/nproc;
	int b = 128;
	int g = (seq_len + (b - 1)) / b;
	int per_block_seq_len = (seq_len + (g - 1)) / g;

	//*CUDA call *
	gpu_parseKmerNFillupBuff<<<g, b>>>(d_seq, d_kmers, klen, seq_len, d_outgoing, d_owner_counter, nproc, p_buff_len);

	// h_outgoing = (uint64_t *) malloc ( n_kmers * BUFF_SCALE * sizeof(uint64_t));
	checkCuda (cudaMemcpy(&(h_outgoing[0]), d_outgoing, n_kmers * BUFF_SCALE * sizeof(uint64_t) , cudaMemcpyDeviceToHost), __LINE__); 
	checkCuda (cudaMemcpy(owner_counter.data(), d_owner_counter, nproc * sizeof(int) , cudaMemcpyDeviceToHost), __LINE__); 

	uint64_t total_counter = 0;
	// printf("GPU ParseNPack: outgoing buufers: ");
	for (int i = 0; i < nproc; ++i) {   
		total_counter += owner_counter[i];    
		// printf(" %d", owner_counter[i]);
	}
	// cout << endl;
	cudaFree(d_seq);
	cudaFree(d_outgoing);
	cudaFree(d_owner_counter);

	cuda_timer_stop(start, stop, milliseconds);
	return;// h_outgoing;
}


// double Exchange_GPUsupermers(vector<keyType>& outgoing, vector<unsigned char> &len_smers, int *sendcnt, int *recvcnt, int n_kmers);
// double GPUtime=0;
// uint64_t *sendbuf_GPU = NULL;
// int * owner_counter; 
// int nkmers_thisBatch = 0;
// int nkmers_processed = 0;
// size_t nkmersAllseq_thisBatch = 0;
std::unordered_map<uint64_t,uint64_t> kcounter_cpu; 
// keyType *device_keys;
// double tot_alltoallv_GPU = 0;
// size_t nkmers_gpu = 0;

// void insert_hashtable_cpu(keyType* mykmers_GPU, int nkeys){

// 	for (uint32_t i = 0; i < nkeys; i++){
// 		uint64_t longs =  mykmers_GPU[i];
// 		auto found = kcounter_cpu.find(longs);// == kcounter_cpu.end() )
// 		if(found != kcounter_cpu.end())
// 			found->second += 1;
// 		else 
// 			kcounter_cpu.insert({longs,1});
// 	}

// 	uint64_t totalPairs = 0, HTsize= 0;
// 	for ( auto it = kcounter_cpu.begin(); it!= kcounter_cpu.end(); ++it ){
// 		if(it->second > 0){
// 			HTsize++;
// 			totalPairs += it->second;
// 		}
// 	}

// 	size_t allrank_hashsize = 0, allrank_totalPairs = 0;
// 	CHECK_MPI( MPI_Reduce(&HTsize, &allrank_hashsize,  1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
// 	CHECK_MPI( MPI_Reduce(&totalPairs, &allrank_totalPairs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );

// 	size_t allrank_kmersthisbatch = 0;
// 	size_t allrank_kmersprocessed = 0;
// 	// CHECK_MPI( MPI_Reduce(&nkmers_thisBatch, &allrank_kmersthisbatch, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
// 	// CHECK_MPI( MPI_Reduce(&nkmers_sofar, &allrank_kmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
// 	if(myrank == 0)
// 		cout << allrank_hashsize << ", #kmers from CPU_HT: " << allrank_totalPairs << endl;

// }


// int getKmers_test(char * seq) {

// 	int * tmp_counter = new int[nprocs];
// 	unsigned int seq_len = strlen(seq);
// 	unsigned int n_kmers =  seq_len - KMER_LENGTH + 1;

// 	for (int i = 0; i < nprocs; ++i)	tmp_counter[i] = 0;

// 	std::vector<uint64_t> kmers_compressed;
// 	int na = 0;

// 	for(int i = 0; i < seq_len - KMER_LENGTH + 1; ++i) {
// 		int kk = i;
// 		uint64_t longs = 0;
// 		bool validKmer = true;
// 		for (int k = 0; k < KMER_LENGTH; ++k,++kk)      {
// 			char s = seq[kk];
// 			if(s == 'a' || s == 'N') { 
// 				na++; 
// 				validKmer = false; break;
// 			}
// 			int j = k % 32;
// 			int l = k/32;
// 			// assert(s != '\0'); 
// 			size_t x = ((s) & 4) >> 1;
// 			longs |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer
// 		}
// 		if(validKmer){
// 			// keyType	owner = murmur3_64(longs);
// 			uint64_t myhash_1 = tmp_MurmurHash3_x64_128((const void *)&longs, 8 , 313); 
// 			double range = static_cast<double>(myhash_1) * static_cast<double>(nprocs);
// 			size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());
// 			tmp_counter[owner]++;
// 			kmers_compressed.push_back(longs);    
// 		}
// 	}
// 	for (int i = 0; i < nprocs; ++i)
// 	{
// 		cout << tmp_counter[i] << " ";
// 	}

// 	// cout << "Similar CPU Parse & pack: Total kmers: " <<  kmers_compressed.size() << " - " << na << endl;// if(i < 5)// toString_custm(kmers_compressed[i]);  
// 	return kmers_compressed.size(); 
// }


// double GPU_DealWithInMemoryData(KeyValue * pHashTable, vector<keyType> & mykmers_GPU, int pass, struct bloom * bm){

// 	// std::vector<keyType> insert_kvs = populate_GPUarray(mykmers_GPU); 
// 	double t_count = MPI_Wtime();
// 	// Insert items into the hash table
// 	// insert_hashtable(pHashTable, mykmers_GPU.data(), (uint32_t)mykmers_GPU.size(), myrank);
// 	insert_hashtable_cpu(mykmers_GPU.data(), (uint32_t)mykmers_GPU.size());
// 	size_t num_keys = nkmers_gpu;
// 	insert_hashtable(pHashTable, device_keys, num_keys, myrank);

// 	double t_count1 = MPI_Wtime() - t_count;

// 	// correctness check
// 	std::vector<KeyValue> h_pHashTable(kHashTableCapacity);
// 	// h_pHashTable.resize(kHashTableCapacity);
// 	cudaMemcpy(h_pHashTable.data(), pHashTable, sizeof(KeyValue) * kHashTableCapacity, cudaMemcpyDeviceToHost); 

// 	uint64_t HTsize = 0, totalPairs= 0;  
// 	for (int i = 0; i < kHashTableCapacity; ++i)
// 	{
// 		// if (hosthash[i].value > 1 && hosthash[i].value < 8) { //kEmpty 
// 		if (h_pHashTable[i].value > 0) { //kEmpty  		
// 			HTsize++;
// 			totalPairs += h_pHashTable[i].value;
// 		}
// 	}

// 	size_t allrank_hashsize = 0, allrank_totalPairs = 0;
// 	CHECK_MPI( MPI_Reduce(&HTsize, &allrank_hashsize,  1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
// 	CHECK_MPI( MPI_Reduce(&totalPairs, &allrank_totalPairs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );

// 	size_t allrank_kmersthisbatch = 0;
// 	size_t allrank_kmersprocessed = 0;
// 	CHECK_MPI( MPI_Reduce(&nkmers_thisBatch, &allrank_kmersthisbatch, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
// 	CHECK_MPI( MPI_Reduce(&nkmers_processed, &allrank_kmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
// 	// cout << myrank << " local GPU HT size " << HTsize << " pair " << totalPairs << endl;
// 	if(myrank == 0){
// 		cout << "\n GPU HTsize: " 
// 			<< allrank_hashsize << ", #kmers from GPU_HT: " << allrank_totalPairs 
// 			<< " ideal #kmers " << allrank_kmersprocessed << endl;
// 	}
// 	return t_count1;
// }

// // double tot_alltoallv_GPU = 0;
// double GPU_Exchange(vector<keyType> & mykmers_GPU, int pass, int nprocs)
// {
// 	MPI_Pcontrol(1,"Exchange");
// 	double tot_exch_time = MPI_Wtime();
// 	double performance_report_time = 0.0; 
// 	// mpicomm mcomm;
// 	// size_t bytesperkmer = Kmer::numBytes();
// 	// size_t bytesperentry = bytesperkmer + (pass == 2 ? sizeof(ReadId) + sizeof(PosInRead) : 0);
// 	int * sendcnt = new int[nprocs];
// 	int * sdispls = new int[nprocs];
// 	int * rdispls = new int[nprocs];
// 	int * recvcnt = new int[nprocs];

// 	int size = nprocs;
// 	for (int i=0; i < nprocs; i++) {
// 		// outgoing[i].clear(); // CPU parseNPack is populating this..remove when that call is removed
// 		sendcnt[i] = owner_counter[i];
// 	}
// 	free(owner_counter);

// 	CHECK_MPI( MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD) );  // share the request counts

// 	int64_t totsend = accumulate(sendcnt, sendcnt+nprocs, static_cast<int64_t>(0));
// 	if (totsend < 0) { cerr << myrank << " detected overflow in totsend calculation, line" << __LINE__ << endl; }
// 	int64_t totrecv = accumulate(recvcnt, recvcnt+nprocs, static_cast<int64_t>(0));
// 	if (totrecv < 0) { cerr << myrank << " detected overflow in totrecv calculation, line" << __LINE__ << endl; }
// 	DBG("totsend=%lld totrecv=%lld\n", (lld) totsend, (lld) totrecv);

// 	int n_kmers = nkmersAllseq_thisBatch;
// 	int p_buff_len = ((n_kmers * BUFF_SCALE) + nprocs - 1)/nprocs;

// 	for (int i=0; i < nprocs; i++) {
// 		sdispls[i] = i * p_buff_len;
// 		rdispls[i] = i * p_buff_len;
// 	}

// 	uint64_t* recvbuf = (uint64_t*) malloc(n_kmers * BUFF_SCALE * sizeof(uint64_t)); 

// 	double exch_time = MPI_Wtime(); //so weird
// 	for (int i = 0; i < COMM_ITER; ++i)
// 		CHECK_MPI( MPI_Alltoallv(sendbuf_GPU, sendcnt, sdispls, MPI_UINT64_T, recvbuf, recvcnt, rdispls, MPI_UINT64_T, MPI_COMM_WORLD) );	
// 	exch_time = (MPI_Wtime()-exch_time) /COMM_ITER;
// 	tot_alltoallv_GPU += exch_time;

// 	****** Performance reporting ******
// 	// performance_report_time = MPI_Wtime();
// 	// const int SND=0, RCV=1;
// 	// int64_t local_counts[2];
// 	// local_counts[SND] = totsend;
// 	// local_counts[RCV] = totrecv;

// 	// int64_t global_mins[2]={0,0};
// 	// CHECK_MPI( MPI_Reduce(&local_counts, &global_mins, 2, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD) );
// 	// int64_t global_maxs[2]={0,0};
// 	// CHECK_MPI( MPI_Reduce(&local_counts, &global_maxs, 2, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD) );
// 	// int64_t global_sums[2] = {0,0};
// 	// CHECK_MPI( MPI_Reduce(&local_counts, &global_sums, 2, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );

// 	// double global_min_time = 0.0;
// 	// CHECK_MPI( MPI_Reduce(&exch_time, &global_min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
// 	// double global_max_time = 0.0;
// 	// CHECK_MPI( MPI_Reduce(&exch_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
// 	// double global_sum_time = 0.0;
// 	// CHECK_MPI( MPI_Reduce(&exch_time, &global_sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );

// 	// int bytedata = 8;//sizeof(uint64_t);
// 	// // serial_printf("KmerMatch:%s sent min %lld bytes, sent max %lld bytes, sent avg %lld bytes, recv min %lld bytes, \
// 	// // 	recv max %lld bytes, recv avg %lld bytes, in min %.3f s, max %.3f s, avg %.3f s\n", __FUNCTION__, global_mins[SND]*bytedata, \
// 	// // 	global_maxs[SND]*bytedata, (global_sums[SND]/nprocs)*bytedata, global_mins[RCV]*bytedata, global_maxs[RCV]*bytedata, (global_sums[RCV]/nprocs)*bytedata, \
// 	// // 	global_min_time, global_max_time, global_sum_time/nprocs);
// 	// performance_report_time = MPI_Wtime()-performance_report_time;
// 	/**************/

// 	nkmers_gpu = 0;
// 	   for(uint64_t i= 0; i < nprocs ; ++i) {
// 			for(uint64_t j= 0; j <  recvcnt[i] ; ++j)				 
// 				mykmers_GPU.push_back(recvbuf[i * p_buff_len + j]);
// 	}

// 	cudaMalloc(&device_keys, sizeof(keyType) * totrecv); 

// 	for(uint64_t i= 0; i < nprocs ; ++i) {
// 		cudaMemcpy(device_keys + nkmers_gpu, &recvbuf[i * p_buff_len], sizeof(keyType) * recvcnt[i], cudaMemcpyHostToDevice); 
// 		nkmers_gpu += recvcnt[i];	
// 	}

// 	if(totsend > 0)  free(sendbuf_GPU);
// 	if(totrecv > 0)  free(recvbuf);

// 	DeleteAll(rdispls, sdispls, recvcnt, sendcnt);
// 	tot_exch_time=MPI_Wtime()-tot_exch_time-performance_report_time;
// 	return tot_exch_time;
// }

// size_t GPU_ParseNPack(vector<string> & seqs, int pass, size_t offset, size_t endoffset, int nprocs)
// {
// 	size_t nreads = endoffset;// seqs.size();
// 	size_t nskipped = 0;
// 	size_t maxsending = 0, kmersthisbatch = 0;
// 	size_t memoryThreshold = (MAX_ALLTOALL_MEM / nprocs) * 2; // 2x any single rank
// 	std::string all_seqs;
// 	uint64_t all_seq_size = 0;
	
// 	int * recvcnt = new int[nprocs];
// 	int * sendcnt = new int[nprocs];
// 	memset(sendcnt, 0, sizeof(int) * nprocs);
// 	nkmers_thisBatch = 0;
// 	nkmersAllseq_thisBatch = 0;
	
// 	for(size_t i=offset; i< nreads; ++i){
// 		if (seqs[i].length() <= KMER_LENGTH) continue;
// 		all_seq_size += seqs[i].length();
// 	}
// 	if(all_seq_size == 0) return nreads;

// 	for(size_t i=offset; i < nreads; ++i){	
// 		size_t found  = seqs[i].length();
		
// 		if (seqs[i].length() <= KMER_LENGTH) continue;

// 		nkmers_thisBatch += seqs[i].length() - KMER_LENGTH + 1;
// 		// int approx_maxsending = (4 * nkmers_thisBatch + nprocs - 1)/nprocs;
// 		all_seqs += seqs[i] + "a";

// 		// if (approx_maxsending * bytesperentry >= memoryThreshold || (nkmers_thisBatch + seqs[i].length()) * bytesperentry >= MAX_ALLTOALL_MEM) { 
// 		// 	nreads = i+1; // start with next read
// 		// 	break;
// 		// }
// 	}
// 	owner_counter = (int*) malloc (nprocs * sizeof(int)) ;
// 	memset(owner_counter, 0, nprocs * sizeof(int));

// 	sendbuf_GPU = getKmers_GPU(all_seqs, KMER_LENGTH, nprocs, owner_counter, myrank, BUFF_SCALE);
// 	nkmersAllseq_thisBatch = all_seqs.length() - KMER_LENGTH + 1;
// 	// getKmers_test(all_seqs);
// 	nkmers_processed += nkmers_thisBatch;

// 	return nreads;
// }

// Lookup keys in the hashtable, and return the values
__global__ void gpu_hashtable_lookup(KeyValue* hashtable, KeyValue* kvs, unsigned int numkvs)
{
	unsigned int threadid = blockIdx.x * blockDim.x + threadIdx.x;
	if (threadid < kHashTableCapacity)
	{
		uint32_t key = kvs[threadid].key;
		uint32_t slot = cuda_murmur3_32(key) & (kHashTableCapacity-1);

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
		uint32_t slot = cuda_murmur3_32(key) & (kHashTableCapacity-1);

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
