
#ifndef LINEARPROBING_H
#define LINEARPROBING_H

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


KeyValue* create_hashtable_GPU(int rank);

void insert_hashtable(KeyValue* pHashTable, vector <keyType>& recvbuf, uint32_t nkmers, int rank, int nprocs, int *recvcnt, int p_buff_len);

void lookup_hashtable(KeyValue* hashtable, KeyValue* kvs, uint32_t num_kvs);

void delete_hashtable(KeyValue* hashtable, const KeyValue* kvs, uint32_t num_kvs);

std::vector<KeyValue> iterate_hashtable(KeyValue* hashtable);

void destroy_hashtable(KeyValue* hashtable, int rank);

void getKmers_GPU(string seq, std::vector<uint64_t> & h_outgoing, int klen, int nproc, vector<int> & owner_counter, int rank, int BUFF_SCALE);
size_t KC_GPU(vector<string> & seqs, int pass, size_t offset, size_t endoffset, int nproc);
// double GPU_buildCounter(KeyValue * pHashTable, vector<keyType> & mykmers_GPU, int pass, struct bloom * bm);
// double GPU_Exchange(vector<keyType> & mykmers_GPU, int pass, int nproc);
#endif
