#include "Kmer.hpp"
#include "Friends.h"
#include "MPIType.h"
#include "SimpleCount.h"
#include "FriendsMPI.h"
#include "Pack.h"
#include "supermer.h"


// int MINIMIZER_LENGTH = 5;
// int mlen = MINIMIZER_LENGTH;

__device__ keyType find_minimizer(keyType kmer, int &order, int klen, int mlen, keyType max64){
	
	keyType minimizer = max64;
	int local_order;

	for (int m = 0; m < (klen - mlen + 1); ++m){
		keyType mmer = ((kmer >> (2*(31-(m+mlen-1)))) );
		
		if( mmer < minimizer ) {
			local_order = m;
			minimizer = mmer;
		}
	}
	order += local_order;
	return minimizer;
}

__device__ keyType find_minimizer(char* kmer, int &order, int klen, int mlen){
	
	char *minimizer = "ZZZZZ";
	keyType com_mini = 0;/////////INT_MAX;
	// int local_order;
	// char mmer[5];
	
	// for (int m = 0; m < (klen - mlen); ++m){
		
	// 	for (int i = 0; i < mlen; ++i)
	// 		mmer[i] = kmer[m+i];
	// 	if( mmer < minimizer ) {
	// 		local_order = m;
	// 		minimizer = mmer;
	// 	}
	// }
	// order += local_order;

	// for (int m = 0; m < mlen; ++m) {
 //        char s = mmer[m];//minimizer[m];
 //        int j = m % 32;
 //        size_t x = ((s) & 4) >> 1;
 //        com_mini |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer 	
 //    }

	return com_mini;
}

// size_t total_kmers = 0, total_supermers = 0, tot_char = 0;



__global__ void cuda_build_supermer(char *seq, char *kmers, int klen, int mlen, unsigned int seq_len,
    keyType* outgoing, unsigned char *out_slen, int *owner_counter, int nproc, unsigned int p_buff_len, int per_block_seq_len){
    
    unsigned int tId = threadIdx.x;
    unsigned int laneId = tId & (blockDim.x - 1);
    unsigned int gId = (blockIdx.x * blockDim.x + tId);
    int window = klen - mlen + 1 ;    
    // int per_block_seq_len = blockDim.x * window;   
    // int per_block_seq_len = (seq_len + (gridDim.x - 1)) / gridDim.x;
    int st_char_block = blockIdx.x * per_block_seq_len; //first char this block should read
    int nKmer = seq_len - klen + 1; //last char is 'a'
    keyType max64 = 18446744073709551615;

    for(int i = st_char_block + laneId * window; i < (st_char_block + per_block_seq_len) && i < nKmer ; i+=blockDim.x*window) {
     // if(tId == 0)
     //    printf("GPU smer info %d %d %d %d \n", gridDim.x,  blockIdx.x, per_block_seq_len, st_char_block );  
 
        keyType longs = 0; //GPU CAS support this for 64 bit
        bool validKmer = true, inserted = false;
        int slen = klen;
        keyType comprs_Smer = 0;
        keyType comprs_Kmer = 0;

        keyType cur_mini = max64;  keyType prev_mini = cur_mini;	
        
        int order = 0, prev_order = 0;
        
        for (int w = 0; w < window && (i+ w + klen) < nKmer; ++w){ //make it linear time	
       		inserted = false; validKmer = true;
        	
        	for (int k = 0; k < klen ; ++k) {
	            char s =  seq[i + w + k ];
	            if(s == 'a' || s == 'N')  { //improvement scope..discard a chunk based on loc of N
	                validKmer = false; break;
	            }
	            int j = k % 32;
	            size_t x = ((s) & 4) >> 1;
	            comprs_Kmer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer 	
	            order = i + w;
	            // cur_kmer[k] = seq[i + w + k ];
	        }
			if(validKmer)
            {
                if(w == 0) //
                { //not initilized yet with any mini
                    // cur_mini = find_minimizer(cur_kmer, order, klen, mlen);//find_minimizer(comprs_Kmer, order, klen, mlen);
                    cur_mini = find_minimizer(comprs_Kmer, order, klen, mlen, max64);
                    comprs_Smer = comprs_Kmer; slen = klen;
                     // if(gId == 0) printf(" brand new smer %lu %lu \n", cur_mini, comprs_Smer );   
                }
                else 
                {
		  	     	
                    cur_mini = find_minimizer(comprs_Kmer, order, klen, mlen, max64);;//find_minimizer(cur_kmer, order, klen, mlen);
		        	
		        	if(prev_mini == cur_mini ){	
                        // printf("mini match  %lu %lu \n", cur_mini, comprs_Smer ); 	     
			        	char s =  seq[i + w + klen - 1];
			            int j = slen % 32; 
			            size_t x = ((s) & 4) >> 1;
			            comprs_Smer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j)));
			            slen++;
		        	}
			     	else 
                    {		                  
		            	keyType owner = cuda_murmur3_64(comprs_Smer) & (nproc - 1); // remove & with HTcapacity in func
		            	int old_count = atomicAdd(&owner_counter[owner],1); 
			            // if(gId == 0) printf("%lu %d\n", owner, old_count );
                        if(old_count >= p_buff_len )  { 
                            printf("%d %d \n", old_count, p_buff_len );
                            printf("Overflow!! MISSION ABORT!!\n");
                        }           	
                        outgoing[owner * p_buff_len + old_count] = comprs_Smer; //hash (longs) 
		            	out_slen[owner * p_buff_len + old_count] = slen;  
		            	comprs_Smer = comprs_Kmer;
                        slen = klen;
		            	// if(w == (window - 1) || (i == nKmer - 1)) inserted = true;
			        }
			    }
			    prev_mini = cur_mini;
		    	prev_order = order;
			}
			// else cur_mini == max64;
	    }
	    if(validKmer)
        {
	    	keyType owner = cuda_murmur3_64(comprs_Smer) & (nproc - 1); // remove & with HTcapacity in func
        	int old_count = atomicAdd(&owner_counter[owner],1); 
            if(old_count >= p_buff_len)   printf("Overflow!! MISSION ABORT!!\n");     
        	outgoing[owner * p_buff_len + old_count] = comprs_Smer; //hash (longs)
        	out_slen[owner * p_buff_len + old_count] = slen;   
	    }
    }
}
int buff_scale = 2;

void getSupermers_GPU(char* seq, int klen, int mlen, int nproc, int *owner_counter, 
	keyType* h_send_smers, unsigned char* h_send_slens, int n_kmers, int rank )
{
    klen= 17;
	int count, devId;
    char *d_kmers, *d_seq;
    keyType *d_supermers, *d_outOverflowBuff;
    unsigned char *d_slen;
    int *d_owner_counter; 
        
    // Map MPI ranks to GPUs
    cudaGetDeviceCount(&count);
    int gpuID = rank % count;
    cudaSetDevice(gpuID);
    cudaGetDevice(&devId);

    unsigned int seq_len = strlen(seq);
    if(seq_len < klen) return;// h_outgoing;
    // unsigned int n_kmers =  seq_len - klen + 1;

    // Create events for GPU timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float milliseconds = 0;
    
    cuda_timer_start(start);


    // CUDA mallocs
    checkCuda (cudaMalloc(&d_supermers, n_kmers * buff_scale * sizeof(keyType)), __LINE__);  // giving 2x space to each node 
    checkCuda (cudaMalloc(&d_slen, n_kmers * buff_scale * sizeof(unsigned char)), __LINE__);  // giving 2x space to each node 
    checkCuda (cudaMalloc(&d_seq, seq_len * sizeof(char)), __LINE__);
    checkCuda (cudaMalloc(&d_owner_counter, nproc * sizeof(int)), __LINE__);
    // CUDA memcopies
    checkCuda (cudaMemcpy(d_seq, seq, seq_len * sizeof(char) , cudaMemcpyHostToDevice), __LINE__);
    cudaMemset(d_supermers,  0, n_kmers * buff_scale * sizeof(keyType));
    cudaMemset(d_owner_counter,  0, sizeof(int) * nproc);

	int window = klen - mlen + 1 ;
	unsigned int p_buff_len = ((n_kmers * buff_scale) + nproc - 1)/nproc;
    int b = 128;
    int g = (seq_len + (b -1) ) / (b);// * window;
    // g = (g+window-1)/window;
    int per_block_seq_len = b * window;// ((seq_len+window-1/window) + (g - 1)) / g;
    // cout << seq_len << "  info " << g << " " << per_block_seq_len << endl;

    // Kernel call
    cuda_build_supermer<<<g, b>>>(d_seq, d_kmers, klen, mlen, seq_len, d_supermers, d_slen, d_owner_counter, nproc, p_buff_len, per_block_seq_len);

    // h_outgoing = (keyType *) malloc ( n_kmers * buff_scale * sizeof(keyType));
    //***** copy back to CPU *****
    checkCuda (cudaMemcpy(h_send_smers, d_supermers, n_kmers * buff_scale * sizeof(keyType), cudaMemcpyDeviceToHost), __LINE__); 
    checkCuda (cudaMemcpy(h_send_slens, d_slen, n_kmers * buff_scale * sizeof(unsigned char), cudaMemcpyDeviceToHost), __LINE__); 
    checkCuda (cudaMemcpy(owner_counter, d_owner_counter, nproc * sizeof(int) , cudaMemcpyDeviceToHost), __LINE__); 
   
    size_t total_counter = 0;
    for (int i = 0; i < nproc; ++i) {   
        total_counter += owner_counter[i];    
        // printf("GPU Supermer pack: output buffer: %d %d \n", owner_counter[i], total_counter);
    }
    cudaFree(d_seq);
    cudaFree(d_supermers);
    cudaFree(d_slen);
    cudaFree(d_owner_counter);

    cuda_timer_stop(start, stop, milliseconds);
    return;
}


__global__ void cu_kcounter_smer(KeyValue* hashtable, const keyType* kvs, const unsigned char* slens,  unsigned int numkvs, int klen, keyType mask)
{
    unsigned int threadid = blockIdx.x*blockDim.x + threadIdx.x;
   
    if (threadid < numkvs){

        keyType new_smer = kvs[threadid];
        unsigned char c = slens[threadid];
        int slen = (int)c;

           
        for(int k = 0; k < (slen - klen + 1); ++k)
        {
            keyType new_key = ((new_smer) >> (2*(31-(klen+k -1)))) & mask;//kvs[threadid];//.key;
             
            keyType slot = cuda_murmur3_64(new_key) & (kHashTableCapacity-1);
            
            while (true){
                keyType old_key = atomicCAS(&hashtable[slot].key, kEmpty, new_key);
                      
                if (old_key == kEmpty || old_key == new_key) {
                    atomicAdd(&hashtable[slot].value,1);
                    break;
                }
                slot = (slot + 1) & (kHashTableCapacity-1);
            }
        }
    }
}

void kcounter_supermer_GPU(KeyValue* pHashTable, keyType* d_smers, unsigned char* d_slen, uint32_t num_keys, int klen, int rank)
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
	const uint32_t N = 1;
    /*------------------------
     Copy kmers to GPU      
    ------------------------*/
    // keyType* device_keys;
    cudaEventRecord(start);

    int b = 128;
    int g= (N + (b - 1)) / b;
   
    /*----------------------------
    CUDA call: Insert kmers to HT       
    ------------------------------*/

    // Have CUDA calculate the thread block size
    keyType mask = pow(2, 2 * klen) - 1;
    int mingridsize;
    int threadblocksize;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, cu_kcounter_smer, 0, 0);


    unsigned char * h_slens = (unsigned char *) malloc(num_keys * sizeof(unsigned char));
    checkCuda (cudaMemcpy(h_slens, d_slen, num_keys * sizeof(unsigned char), cudaMemcpyDeviceToHost), __LINE__); 
   // keyType * h_smers = (keyType *) malloc(num_keys * sizeof(keyType));
   //  checkCuda (cudaMemcpy(h_smers, d_smers, num_keys * sizeof(keyType), cudaMemcpyDeviceToHost), __LINE__); 
   
   cout << "after memcopies: ";
for (int i = 0; i < num_keys; ++i)
{
    if(i < 10 ) cout << h_slens[i] << " ";
    // if(h_slens[i] > 17)
    // cout << i << ": " << h_slens[i] << " larger ";
}


    int gridsize = ((uint32_t)num_keys + threadblocksize - 1) / threadblocksize;
    cu_kcounter_smer<<<gridsize, threadblocksize>>>(pHashTable, d_smers, d_slen, (uint32_t)num_keys, klen, mask);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaFree(d_smers);
    cudaFree(d_slen);
    return ;//h_pHashTable;

}


	