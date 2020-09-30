
#ifndef SUPERMER_H
#define SUPERMER_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include "common_gpu.h"
using namespace std;


typedef unsigned long long int uint64_cu; 
typedef unsigned long long int keyType; 
// struct KeyValue 
// {
//     keyType key;
//     uint32_t value;
// };




size_t build_supermer(vector<string> seqs, size_t offset); 
void getSupermers_GPU(char* seq, int klen, int mlen, int nproc, int *owner_counter, 
	keyType* h_send_smers, unsigned char* h_send_slens, int n_kmers, int rank );
void kcounter_supermer_GPU(KeyValue* pHashTable, keyType* d_smers, unsigned char* d_slen, uint32_t num_keys, int klen, int rank);
#endif
