// #pragma once

// #define keyType uint32_t

#define keyType unsigned long long int
// #define keyType unsigned int
struct KeyValue 
{
    keyType key;
    uint32_t value;
};

// struct HashTable_GPU
// {
	const uint64_t kHashTableCapacity = 128 * 1024 * 1024;

	const uint64_t kNumKeyValues = kHashTableCapacity / 2;//4168188;

	const keyType kEmpty = 0;

	// HashTable_GPU(uint64_t HT_cap, uint64_t size){
	// 	kHashTableCapacity = HT_cap;
	// 	kNumKeyValues = kHashTableCapacity / 2 ;
	// }
	
// };


// void init_HT_GPU()
// {
// 	const uint64_t kHashTableCapacity = 128 * 1024 * 1024;

// 	const uint64_t kNumKeyValues = kHashTableCapacity / 2;

// 	const keyType kEmpty = 0;//0x1ffffffffu; //0xffffffff;
	
// };



KeyValue* create_hashtable_GPU(int rank);

std::vector<KeyValue> insert_hashtable(KeyValue* hashtable, const keyType* kvs, uint32_t num_kvs, int rank);

void lookup_hashtable(KeyValue* hashtable, KeyValue* kvs, uint32_t num_kvs);

void delete_hashtable(KeyValue* hashtable, const KeyValue* kvs, uint32_t num_kvs);

std::vector<KeyValue> iterate_hashtable(KeyValue* hashtable);

void destroy_hashtable(KeyValue* hashtable, int rank);

uint64_t * getKmers_GPU( char *seq, int klen, int nproc, int rank);
