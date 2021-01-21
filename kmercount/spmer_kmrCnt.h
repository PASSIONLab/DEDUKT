
#ifndef SP_KC_H
#define SP_KC_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include "common.h"
#include "common_gpu.h"
using namespace std;

typedef unsigned long long int keyType;


uint64_t murmur3_64 (uint64_t k);
int murmur3_32 (int k);
uint64_t find_minimizer (uint64_t kmer, int mlen);

size_t spmer_kmerCount(vector<string> seqs, size_t offset, size_t endoffset, 
			int klen, int mlen);

void getSupermers_GPU (string seq, int klen, int mlen, int nproc,
		       int *owner_counter, vector <keyType> &h_send_smers,
		       vector <unsigned char> &h_send_slens, int n_kmers, int rank,
		       int BUFF_LEN);
void GPU_SP_buildCounter(KeyValue* pHashTable, vector<keyType> &recvbuf, vector<unsigned char> &recvbuf_len,
		int * recvcnt, uint32_t num_keys, int klen, int rank, int p_buff_len);
double Exchange_GPUsupermers(vector<keyType> &outgoing, vector<unsigned char> &len_smers, 
	vector<keyType> &recvbuf, vector<unsigned char> &recvbuf_len,
	int *sendcnt, int *recvcnt, int nkmers,  int * owner_counter);


#endif
