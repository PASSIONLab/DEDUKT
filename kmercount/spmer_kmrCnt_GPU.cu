#include "Kmer.hpp"
#include "Friends.h"
#include "MPIType.h"
#include "SimpleCount.h"
#include "FriendsMPI.h"
#include "Pack.h"
#include "spmer_kmrCnt.h"

__device__ keyType find_minimizer(keyType kmer, int klen, int mlen, keyType max64){

	keyType minimizer = max64;
	keyType mask = pow(2, 2 * mlen) - 1;

	for (int m = 0; m < (klen - mlen ); ++m){
		keyType mmer =  (kmer >> (2*(31-(mlen+m -1)))) & mask;

		if( mmer < minimizer ) 
			minimizer = mmer;
	}
	return minimizer;
}

__global__ void cuda_build_supermer(char *seq, char *kmers, int klen, int mlen, unsigned int seq_len,
		keyType* outgoing, unsigned char *out_slen, int *owner_counter, int nproc, unsigned int p_buff_len, 
		int per_block_seq_len, int window, int rank){

	unsigned int tId = threadIdx.x;
	unsigned int laneId = tId & (blockDim.x - 1);
	unsigned int gId = (blockIdx.x * blockDim.x + tId);
	int st_char_block = blockIdx.x * per_block_seq_len; //first char this block should read
	int nKmer = seq_len - klen + 1; //last char is 'a'
	keyType max64 = 18446744073709551615;

	bool validKmer = true;
	int slen = klen;
	keyType comprs_Smer = 0;
	keyType comprs_Kmer = 0;
	int owner = -1;
	int old_count=-1;

	keyType cur_mini = max64;  keyType prev_mini = cur_mini;    

	//****First kmer of this window *****
	int i = st_char_block + laneId * window;

	if(i <  nKmer) {

		comprs_Kmer = 0;
		for (int k = 0; k < klen ; ++k) {
			char s =  seq[i + k ];
			if(s == 'a' || s == 'N')  { 
				// w += klen-1;
				validKmer = false; break; //FIXIT can have 'n'
			}
			int j = k % 32;
			size_t x = ((s) & 4) >> 1;
			// comprs_Kmer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer  

			switch(s) { //redefined
				case 'A': comprs_Kmer |= ((x + (x^1)) << (2*(31-j)));  break;
				case 'C': comprs_Kmer |= ((x + (x^0)) << (2*(31-j)));  break;
				case 'G': comprs_Kmer |= ((x + (x^3)) << (2*(31-j)));  break;
				case 'T': comprs_Kmer |= ((x + (x^2)) << (2*(31-j)));  break;
			}
		}
		if(validKmer){
			cur_mini = find_minimizer(comprs_Kmer, klen, mlen, max64);
			prev_mini = cur_mini; 
			// owner = cuda_murmur3_64(cur_mini) & (nproc - 1); // remove & with HTcapacity in func
			// keyType myhash = cuda_murmur3_64(cur_mini); // remove & with HTcapacity in func
			keyType myhash = cuda_MurmurHash3_x64_128((const void *)&cur_mini, 8, 313);// & (nproc - 1);
			double range = static_cast<double>(myhash) * static_cast<double>(nproc);
			owner = range / max64;

			old_count = atomicAdd(&owner_counter[owner],1); 
			outgoing[owner * p_buff_len + old_count] = comprs_Kmer; //hash (longs)
			out_slen[owner * p_buff_len + old_count] = klen;  

		}
		comprs_Smer = comprs_Kmer;
		slen = klen;

		int c = st_char_block + (laneId * window);

		for(int w = 1; w < window && (c+w) < nKmer ; w++) {

			validKmer = true;
			comprs_Kmer = 0;
			// if ((i + klen-1) > nKmer) return;
			for (int k = 0; k < klen ; ++k) {
				char s =  seq[c + w + k ];
				if(s == 'a' || s == 'N')  { 
					// w += klen-1;
					validKmer = false; break;
				}
				int j = k % 32;
				size_t x = ((s) & 4) >> 1;
				// comprs_Kmer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer  
				switch(s) { //redefined
					case 'A': comprs_Kmer |= ((x + (x^1)) << (2*(31-j)));  break;
					case 'C': comprs_Kmer |= ((x + (x^0)) << (2*(31-j)));  break;
					case 'G': comprs_Kmer |= ((x + (x^3)) << (2*(31-j)));  break;
					case 'T': comprs_Kmer |= ((x + (x^2)) << (2*(31-j)));  break;
				}
			}  

			if(validKmer){ 

				cur_mini = find_minimizer(comprs_Kmer, klen, mlen, max64);

				if(prev_mini == cur_mini){ 
					// printf("mini match  %lu %lu \n", cur_mini, comprs_Smer );         
					char s =  seq[c + w + klen - 1];
					int j = slen % 32; 
					size_t x = ((s) & 4) >> 1;
					// comprs_Smer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j)));
					switch(s) { //redefined
						case 'A': comprs_Smer |= ((x + (x^1)) << (2*(31-j)));  break;
						case 'C': comprs_Smer |= ((x + (x^0)) << (2*(31-j)));  break;
						case 'G': comprs_Smer |= ((x + (x^3)) << (2*(31-j)));  break;
						case 'T': comprs_Smer |= ((x + (x^2)) << (2*(31-j)));  break;
					}
					slen++;
				}
				else 	{ 

					if(owner > -1 && old_count > -1)
					{
						outgoing[owner * p_buff_len + old_count] = comprs_Smer; //hash (longs) 
						out_slen[owner * p_buff_len + old_count] = slen;                          
					}
					//* new supermer */
					slen = klen;
					comprs_Smer = comprs_Kmer;
					// owner = cuda_murmur3_64(cur_mini) & (nproc - 1); // remove & with HTcapacity in func
					keyType myhash = cuda_MurmurHash3_x64_128((const void *)&cur_mini, 8, 313);
					// keyType myhash = cuda_murmur3_64(cur_mini); // remove & with HTcapacity in func
					double range = static_cast<double>(myhash) * static_cast<double>(nproc);
					owner = range / max64;

					old_count = atomicAdd(&owner_counter[owner],1); 
					if(old_count > p_buff_len )  { 
						printf("Overflow!! MISSION ABORT!!\n"); return;
					}               
					outgoing[owner * p_buff_len + old_count] = comprs_Smer; //hash (longs) 
					out_slen[owner * p_buff_len + old_count] = slen;  
				}
				prev_mini = cur_mini;
			}
		}   
		if(old_count > -1 && owner > -1) {
			outgoing[owner * p_buff_len + old_count] = comprs_Smer; //hash (longs)
			out_slen[owner * p_buff_len + old_count] = slen; 
		}
	}         
}

void getSupermers_GPU(string seq, int klen, int mlen, int nproc, int *owner_counter, 
		vector<keyType>& h_send_smers, vector<unsigned char>& h_send_slens, int n_kmers, int rank, int BUFF_SCALE )
{

	int count, devId;
	char *d_kmers, *d_seq;
	keyType *d_supermers, *d_outOverflowBuff;
	unsigned char *d_slen;
	int *d_owner_counter; 

	//* Map MPI ranks to GPUs */
	cudaGetDeviceCount(&count);
	int gpuID = rank % count;
	cudaSetDevice(gpuID);
	cudaGetDevice(&devId);

	unsigned int seq_len = seq.length();
	if(seq_len < klen) return;// h_outgoing;

	//* Create events for GPU timing */
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float milliseconds = 0;

	cuda_timer_start(start);

	// CUDA mallocs
	checkCuda (cudaMalloc(&d_supermers, n_kmers * BUFF_SCALE * sizeof(keyType)), __LINE__);  // giving 2x space to each node 
	checkCuda (cudaMalloc(&d_slen, n_kmers * BUFF_SCALE * sizeof(unsigned char)), __LINE__);  // giving 2x space to each node 
	checkCuda (cudaMalloc(&d_seq, seq_len * sizeof(char)), __LINE__);
	checkCuda (cudaMalloc(&d_owner_counter, nproc * sizeof(int)), __LINE__);
	// CUDA memcopies
	checkCuda (cudaMemcpy(d_seq, &seq[0], seq_len * sizeof(char) , cudaMemcpyHostToDevice), __LINE__);
	cudaMemset(d_supermers,  0, n_kmers * BUFF_SCALE * sizeof(keyType));
	cudaMemset(d_owner_counter,  0, sizeof(int) * nproc);

	int window = 32 - klen;// - mlen + 1 ;

	unsigned int p_buff_len = ((n_kmers * BUFF_SCALE) + nproc - 1)/nproc;

	int b = 128;
	int g = (seq_len + (b*window - 1) ) / (b*window); ;//(seq_len + (b -1) ) / (b);// * window;
	int per_block_seq_len = b * window;// ((seq_len+window-1/window) + (g - 1)) / g;

	// Kernel call
	cuda_build_supermer<<<g, b>>>(d_seq, d_kmers, klen, mlen, seq_len, d_supermers, d_slen, d_owner_counter, nproc, p_buff_len, per_block_seq_len, window, rank);

	//* Memcopy to CPU */
	checkCuda (cudaMemcpy(h_send_smers.data(), d_supermers, n_kmers * BUFF_SCALE * sizeof(keyType), cudaMemcpyDeviceToHost), __LINE__); 
	checkCuda (cudaMemcpy(h_send_slens.data(), d_slen, n_kmers * BUFF_SCALE * sizeof(unsigned char), cudaMemcpyDeviceToHost), __LINE__); 
	checkCuda (cudaMemcpy(owner_counter, d_owner_counter, nproc * sizeof(int) , cudaMemcpyDeviceToHost), __LINE__); 

	// size_t total_counter = 0;
	// cout << rank << " smer distribution: ";
	// for (int i = 0; i < nproc; ++i) {   
	//     total_counter += owner_counter[i];    
	//     cout << owner_counter[i] << " "; 
	//     // printf("GPU Supermer pack: output buffer: %d %d \n", owner_counter[i], total_counter);
	// }
	// cout << endl;

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

		//*kmers of the supermer*
		for(int k = 0; k < (slen - klen + 1); ++k){
			
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

void GPU_SP_buildCounter(KeyValue* pHashTable, vector<keyType> &recvbuf, vector<unsigned char> &recvbuf_len,
		int * recvcnt, uint32_t totrecv, int klen, int rank, int p_buff_len)
{
	// Map MPI ranks to GPUs
	int count, devId;
	cudaGetDeviceCount(&count);
	int gpuID = rank % count;
	cudaSetDevice(gpuID);
	// cudaGetDevice(&devId);

	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start);

	int b = 128;
	// Have CUDA calculate the thread block size
	keyType mask = pow(2, 2 * klen) - 1;
	int mingridsize;
	int threadblocksize;
	cudaOccupancyMaxPotentialBlockSize(&mingridsize, &threadblocksize, cu_kcounter_smer, 0, 0);
	unsigned char * d_slens;
	keyType *d_smers;
	// unsigned char * h_slens = (unsigned char *) malloc(num_keys * sizeof(unsigned char));
	// checkCuda (cudaMemcpy(h_slens, d_slen, num_keys * sizeof(unsigned char), cudaMemcpyDeviceToHost), __LINE__); 
	
	checkCuda( cudaMalloc(&d_smers, sizeof(keyType) * totrecv), __LINE__); 
	checkCuda( cudaMalloc(&d_slens, sizeof(unsigned char) * totrecv), __LINE__);
	
	size_t num_keys = 0;

	for(uint64_t i= 0; i < nprocs ; ++i) {
		if(totrecv > 0) {
			checkCuda( cudaMemcpy(d_smers + num_keys, &recvbuf[i * p_buff_len], sizeof(keyType) * recvcnt[i], cudaMemcpyHostToDevice), __LINE__); 
			checkCuda( cudaMemcpy(d_slens + num_keys, &recvbuf_len[i * p_buff_len], sizeof(unsigned char) * recvcnt[i], cudaMemcpyHostToDevice), __LINE__); 
		}
		num_keys += recvcnt[i];	
	}

	int gridsize = ((uint32_t)num_keys + threadblocksize - 1) / threadblocksize;
	cu_kcounter_smer<<<gridsize, threadblocksize>>>(pHashTable, d_smers, d_slens, (uint32_t)num_keys, klen, mask);

	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaFree(d_smers);
	cudaFree(d_slens);
	return ;//h_pHashTable;

}

double tot_GPUsmer_alltoallv = 0;

double Exchange_GPUsupermers(vector<keyType> &outgoing, vector<unsigned char> &len_smers, 
	vector<keyType> &recvbuf, vector<unsigned char> &recvbuf_len,
	int *sendcnt, int *recvcnt, int nkmers,  int * owner_counter)
{
	double tot_exch_time = MPI_Wtime();

	// int * sendcnt = new int[nprocs];
	int * sdispls = new int[nprocs];
	int * rdispls = new int[nprocs];
	// int * recvcnt = new int[nprocs];

	uint64_t totsend = 0, totrecv = 0;
	for (int i=0; i < nprocs; i++) {
		sendcnt[i] = owner_counter[i];
		totsend += sendcnt[i];
	}
	free(owner_counter);

	CHECK_MPI( MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD) );  // share the request counts

	// cout << "recv count " ;
	for (int i=0; i < nprocs; i++) {
		totrecv += recvcnt[i];
		// cout << recvcnt[i] << " "; 
	}
	// cout << endl;

	// int64_t totsend = accumulate(sendcnt, sendcnt+nprocs, static_cast<int64_t>(0));
	// if (totsend < 0) { cerr << myrank << " detected overflow in totsend calculation, line" << __LINE__ << endl; }
	// int64_t totrecv = accumulate(recvcnt, recvcnt+nprocs, static_cast<int64_t>(0));
	// if (totrecv < 0) { cerr << myrank << " detected overflow in totrecv calculation, line" << __LINE__ << endl; }
	// DBG("totsend=%lld totrecv=%lld\n", (lld) totsend, (lld) totrecv);

	int p_buff_len = ((nkmers * BUFF_SCALE) + nprocs - 1)/nprocs;

	for (int i=0; i < nprocs; i++) {
		sdispls[i] = i * p_buff_len;
		rdispls[i] = i * p_buff_len;
	}

	// uint64_t* recvbuf = (uint64_t*) malloc(nkmers * BUFF_SCALE * sizeof(uint64_t)); 
	// unsigned char* recvbuf_len = (unsigned char*) malloc(nkmers * BUFF_SCALE * sizeof(unsigned char)); 

	double exch_time = MPI_Wtime();
	for (int i = 0; i < COMM_ITER; ++i)
	{
		CHECK_MPI( MPI_Alltoallv(&(outgoing[0]), sendcnt, sdispls, MPI_UINT64_T, &(recvbuf[0]), recvcnt, rdispls, MPI_UINT64_T, MPI_COMM_WORLD) );
		CHECK_MPI( MPI_Alltoallv(&(len_smers[0]), sendcnt, sdispls, MPI_UNSIGNED_CHAR, &(recvbuf_len[0]), recvcnt, rdispls, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD) );
	}
	exch_time = (MPI_Wtime() - exch_time)/COMM_ITER;
	tot_GPUsmer_alltoallv += exch_time;

	double performance_report_time = 0;//perf_reporting(exch_time, totsend, totrecv);

	// checkCuda( cudaMalloc(&d_recv_smers, sizeof(keyType) * totrecv), __LINE__); 
	// checkCuda( cudaMalloc(&d_recv_slens, sizeof(unsigned char) * totrecv), __LINE__);
	// size_t num_keys = 0;

	// for(uint64_t i= 0; i < nprocs ; ++i) {
	// 	if(totrecv > 0) {
	// 		checkCuda( cudaMemcpy(d_recv_smers + num_keys, &recvbuf[i * p_buff_len], sizeof(keyType) * recvcnt[i], cudaMemcpyHostToDevice), __LINE__); 
	// 		checkCuda( cudaMemcpy(d_recv_slens + num_keys, &recvbuf_len[i * p_buff_len], sizeof(unsigned char) * recvcnt[i], cudaMemcpyHostToDevice), __LINE__); 
	// 	}
	// 	num_keys += recvcnt[i];	
	// }

	// if(totsend > 0)  {free(outgoing); free(len_smers);}
	// if(totrecv > 0)  {free(recvbuf); free(recvbuf_len);}

	delete(rdispls); delete(sdispls); 
	tot_exch_time=MPI_Wtime()-tot_exch_time; //-performance_report_time;

	return tot_exch_time;
}


