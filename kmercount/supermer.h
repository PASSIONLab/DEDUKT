
#ifndef SUPERMER_H
#define SUPERMER_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include "common_gpu.h"
using namespace std;


typedef unsigned long long int keyType; 


uint64_t murmur3_64(uint64_t k);
int murmur3_32(int k);

uint64_t find_minimizer(uint64_t kmer, int mlen);
size_t build_supermer(vector<string> seqs, size_t offset, size_t endoffset, int KMER_LEN,  int MINIMIZER_LENGTH); 
void getSupermers_GPU(char* seq, int klen, int mlen, int nproc, int *owner_counter, keyType* h_send_smers, unsigned char* h_send_slens, int n_kmers, int rank );
void kcounter_supermer_GPU(KeyValue* pHashTable, keyType* d_smers, unsigned char* d_slen, uint32_t num_keys, int klen, int rank);

#endif
