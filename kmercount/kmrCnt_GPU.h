
#ifndef KC_GPU_H
#define KC_GPU_H

// #include "Friends.h"
#include "common.h"
#include "common_gpu.h"
// #include "Kmer.hpp"
#include "Friends.h"
#include "MPIType.h"
// #include "SimpleCount.h"
#include "FriendsMPI.h"
// #include "Pack.h"
// #include "supermer.h"

extern KeyValue* d_hashTable;
KeyValue* create_hashtable_GPU(int rank);

void insert_hashtable(KeyValue* pHashTable, vector <keyType>& recvbuf, 
	uint32_t nkmers, int rank, int nprocs, int *recvcnt, int p_buff_len);

void lookup_hashtable(KeyValue* hashtable, KeyValue* kvs, uint32_t num_kvs);

void delete_hashtable(KeyValue* hashtable, const KeyValue* kvs, uint32_t num_kvs);

std::vector<KeyValue> iterate_hashtable(KeyValue* hashtable);

void destroy_hashtable(KeyValue* hashtable, int rank);

void getKmers_GPU(string seq, std::vector<uint64_t> & h_outgoing, int klen, 
	int nproc, vector<int> & owner_counter, int rank, int BUFF_SCALE);

double GPU_buildCounter(KeyValue *d_hashTable,  vector <keyType>& recvbuf, int pass, 
	struct bloom * bm, int totrcv, int *recvcnt, int p_buff_len);

double GPU_Exchange(vector <keyType> &recvbuf, vector <uint64_t> &h_outgoing, int pass, 
	int nproc, int n_kmers, std::vector<int> & owner_counter, size_t &nRecvdKmers, 
	int * recvcnt);

#endif
