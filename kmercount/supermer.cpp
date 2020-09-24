// #include <stdio.h>
// #include <iostream>
// #include <vector>
// #include <string>

#include "Kmer.hpp"
#include "supermer.h"
// #include "KmerIterator.hpp"
// #include "Deleter.h"
#include "ParallelFASTQ.h"
#include "Friends.h"
#include "MPIType.h"
#include "SimpleCount.h"
#include "FriendsMPI.h"
#include "Pack.h"
// #include "../common/Buffer.h"
// #include "../common/hash_funcs.h"

using namespace std;


#define BIG_CONSTANT(x) (x)
// 32 bit Murmur3 hash
int murmur3_32( int k )
{
  	k ^= k >> 16;
	k *= 0x85ebca6b;
	k ^= k >> 13;
	k *= 0xc2b2ae35;
	k ^= k >> 16;

  	return k & (nprocs-1);
}

int murmur3_32_tmp( int k )
{
  	k ^= k >> 16;
	k *= 0x85ebca6b;
	k ^= k >> 13;
	k *= 0xc2b2ae35;
	k ^= k >> 16;

  	return k;
}


string decompress_smer(uint64_t c_kmer, int slen){
  
  size_t i,j;
  char *s = new char[slen+1];
  for (i = 0; i < slen; i++) {
    j = i % 32;
    // l = i / 32;

    switch(( (c_kmer) >> (2*(31-j))) & 0x03 ) {
      case 00: s[i] = 'A';  break;
      case 01: s[i] = 'C';  break;
      case 02: s[i] = 'G';  break;
      case 03: s[i] = 'T';  break;  	
    }
  }
  s[slen] = '\0'; 
  string tmp = s;
  return tmp;
}

string decompress_smer(uint64_t c_kmer, unsigned char char_slen){
  
  size_t i,j;
  int slen = (int)char_slen;
  char *s = new char[slen+1];
  for (i = 0; i < slen; i++) {
    j = i % 32;
    // l = i / 32;

    switch(( (c_kmer) >> (2*(31-j))) & 0x03 ) {
      case 00: s[i] = 'A';  break;
      case 01: s[i] = 'C';  break;
      case 02: s[i] = 'G';  break;
      case 03: s[i] = 'T';  break;  	
      // case 0x00: *s = 'A'; ++s; break;
      // case 0x01: *s = 'C'; ++s; break;
      // case 0x02: *s = 'G'; ++s; break;
      // case 0x03: *s = 'T'; ++s; break;  
    }
  }
  s[slen] = '\0'; 
  string tmp = s;
  return tmp;
}

uint64_t compress_smer(string cur_kmer, int len){
	char *s = new char[len];
	strcpy(s, cur_kmer.c_str()); 
	uint64_t c_kmer = 0;
	for (int k = 0; k < len && s[k] != '\0'; ++k)
   	{
		int j = k % 32;
		int l = k/32;
		// assert(s != '\0'); 
		size_t x = ((*s) & 4) >> 1;
		c_kmer |= ((x + ((x ^ (*s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer
		s++;
	}
	return c_kmer;
}

uint64_t compress_kmer(string cur_kmer){
	char *s = new char[KMER_LENGTH];
	strcpy(s, cur_kmer.c_str()); 
	uint64_t c_kmer = 0;
	for (int k = 0; k < KMER_LENGTH; ++k)
   	{
		int j = k % 32;
		int l = k/32;
		// assert(s != '\0'); 
		size_t x = ((*s) & 4) >> 1;
		c_kmer |= ((x + ((x ^ (*s & 2)) >>1)) << (2*(31-j))); //make it longs[] to support larger kmer
		s++;
	}
	return c_kmer;
}


int MINIMIZER_LENGTH = 5;

string find_minimizer(string kmer, int &order){
	
	string minimizer = "ZZZZZZZ";
	int local_order;
	for (int m = 0; m < (KMER_LENGTH - MINIMIZER_LENGTH); ++m)
	{
		string mmer = kmer.substr(m, MINIMIZER_LENGTH);
		if( mmer < minimizer ) {
			local_order = m;
			minimizer = mmer;
		}
	}
	order += local_order;
	return minimizer;
}

size_t total_kmers = 0, total_supermers = 0, tot_char = 0;
int supermer_bins[20];

// void parse_supermer_N_build_kmercounter(vector<uint64_t> c_smers, vector<unsigned char> len_smers){
size_t nkmers_sofar = 0;
void parse_supermer_N_build_kmercounter(vector<uint64_t> c_smer, vector<unsigned char> len_smers, int* recvcnt, int p_buff_len){
   
    std::unordered_map<std::string,uint64_t> kcounter_cpu; 
    int count = 0;
    std::vector<uint64_t> kmers_compressed;
    for (int i = 0; i < c_smer.size(); ++i)
	{
		uint64_t cur_csmer = c_smer[i];// recvbuf[i * p_buff_len + j];
		unsigned char c = len_smers[i];// len_smers[i * p_buff_len + j];
		int slen = (int)c;
   		string d_smer = decompress_smer(cur_csmer, c);
   		
	    for(int k = 0; k < slen - KMER_LENGTH + 1; ++k) {

	    	string cur_kmer = d_smer.substr(k, KMER_LENGTH);
	    	count++;
	    	auto found = kcounter_cpu.find(cur_kmer);// == kcounter_cpu.end() )
	    	if(found != kcounter_cpu.end())
	    		found->second += 1;
	    	else kcounter_cpu.insert({cur_kmer,1}); 
		}
	}

    uint64_t totkmer = 0, nHTsize= 0;
	for ( auto it = kcounter_cpu.begin(); it!= kcounter_cpu.end(); ++it )
	{
		// if(it->second > 1)
		{
		nHTsize++;
		totkmer += it->second;
		}
	}
	
    std::cout << "\ncounter from supermers - HTsize: " << nHTsize 
    	<< " totkmers from HT: " << totkmer << std::endl;
    return;// kmers;
}

std::unordered_map<std::string,uint64_t> kcounter_cpu; 
void parse_supermer_N_build_kmercounter(uint64_t* recvbuf, unsigned char* len_smers, int* recvcnt, int p_buff_len){
   
   std::vector<uint64_t> kmers_compressed;
   // for (int i = 0; i < c_smers.size(); ++i)
   	for(uint64_t i= 0; i < nprocs ; ++i) 
   	{
		for(uint64_t j= 0; j <  recvcnt[i] ; ++j)
   		{ 
   			uint64_t cur_csmer = recvbuf[i * p_buff_len + j];
   			unsigned char c = len_smers[i * p_buff_len + j];
   			int slen = (int)c;
	   		string d_smer = decompress_smer(cur_csmer, c);
	   		// if(j < 5)
	   		// 	cout << myrank << ": " << d_smer << " " << len_smers[j] << endl;
	   		
		    for(int k = 0; k < slen - KMER_LENGTH + 1; ++k) {

		    	string cur_kmer = d_smer.substr(k, KMER_LENGTH);
		    	auto found = kcounter_cpu.find(cur_kmer);// == kcounter_cpu.end() )
		    	if(found != kcounter_cpu.end())
		    		found->second += 1;
		    	else kcounter_cpu.insert({cur_kmer,1}); 
			}
		}
	}

    uint64_t totalPairs = 0, HTsize= 0;
	for ( auto it = kcounter_cpu.begin(); it!= kcounter_cpu.end(); ++it )
	{
		// if(it->second > 1)
		{
		HTsize++;
		totalPairs += it->second;
		}
	}
	
	size_t allrank_hashsize = 0, allrank_totalPairs = 0;
    CHECK_MPI( MPI_Reduce(&HTsize, &allrank_hashsize,  1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
    CHECK_MPI( MPI_Reduce(&totalPairs, &allrank_totalPairs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
  
	size_t allrank_kmersthisbatch = 0;
	size_t allrank_kmersprocessed = 0;
    // CHECK_MPI( MPI_Reduce(&nkmers_thisBatch, &allrank_kmersthisbatch, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
    CHECK_MPI( MPI_Reduce(&nkmers_sofar, &allrank_kmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );

    if(myrank == 0)
    std::cout << "rank: " << myrank << ", Smer+len exch- : HTsize: " << allrank_hashsize 
    	<< " #kmers from HT: " << allrank_totalPairs << ", ideal #kmers: " << allrank_kmersprocessed << std::endl;
    return;// kmers;
}

std::unordered_map<std::string,uint64_t> kcounter_cpu_concat; 
void concat_parse_supermer_N_build_kmercounter(uint64_t* recvbuf, unsigned char* len_smers, int* recvcnt, int* recvcnt_nsmers, int p_buff_len)
{
   string all_supermers;// = NULL;;

   	for(uint64_t i= 0; i < nprocs ; ++i) 
   	{
		for(uint64_t j= 0; j <  recvcnt[i] ; ++j)
   		{
   			uint64_t cur_csmer = recvbuf[i * p_buff_len + j];

	   		all_supermers += decompress_smer(cur_csmer, 32);
   		}
    }
    // cout << "after rv array len " << loc << " " << all_supermers.length() << endl;
    int st_smer = 0;
   	for(uint64_t i= 0; i < nprocs ; ++i) 
   	{
		for(uint64_t j= 0; j <  recvcnt_nsmers[i] ; ++j)
   		{		
	   		unsigned char c = (len_smers[i * p_buff_len + j]);
	   	    int slen = (int)c;
	   		string cur_supermer = all_supermers.substr(st_smer, slen);

	   		for (int k = 0; k < slen - KMER_LENGTH + 1; ++k)
	   		{
	   		    string cur_kmer = cur_supermer.substr(k, KMER_LENGTH);
			    auto found = kcounter_cpu_concat.find(cur_kmer);// == kcounter_cpu_concat; .end() )
			    if(found != kcounter_cpu_concat.end())
			    	found->second += 1;
			    else kcounter_cpu_concat.insert({cur_kmer,1}); 	  
	   		}
			st_smer += slen;
		}
	}

    uint64_t totalPairs = 0, HTsize= 0;
	for ( auto it = kcounter_cpu_concat.begin(); it!= kcounter_cpu_concat.end(); ++it )
	{
		// if(it->second > 1)
		{
		HTsize++;
		totalPairs += it->second;
		}
	}
	
	size_t allrank_hashsize = 0, allrank_totalPairs = 0;
    CHECK_MPI( MPI_Reduce(&HTsize, &allrank_hashsize,  1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
    CHECK_MPI( MPI_Reduce(&totalPairs, &allrank_totalPairs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
  
	size_t allrank_kmersthisbatch = 0;
	size_t allrank_kmersprocessed = 0;
    // CHECK_MPI( MPI_Reduce(&nkmers_thisBatch, &allrank_kmersthisbatch, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
    CHECK_MPI( MPI_Reduce(&nkmers_sofar, &allrank_kmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );

    if(myrank == 0)
    std::cout << "rank: " << myrank << ", Smer+len exch- : HTsize: " << allrank_hashsize 
    	<< " #kmers from concated HT: " << allrank_totalPairs << ", ideal #kmers: " << allrank_kmersprocessed << std::endl;
    return;// kmers;
}

double tot_exch_time_smer = 0; 
uint64_t* Exchange_supermers(uint64_t* outgoing, unsigned char* len_smers, 
	uint64_t* recvbuf, unsigned char* recvbuf_len, int *sendcnt, int *recvcnt, int n_kmers)
{

	MPI_Pcontrol(1,"Exchange");
	double tot_exch_time = MPI_Wtime();
	double performance_report_time = 0.0; 

	// int * sendcnt = new int[nprocs];
	int * sdispls = new int[nprocs];
	int * rdispls = new int[nprocs];
	// int * recvcnt = new int[nprocs];
    
    int size = nprocs;
    
	CHECK_MPI( MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD) );  // share the request counts

	int64_t totsend = accumulate(sendcnt, sendcnt+nprocs, static_cast<int64_t>(0));
	if (totsend < 0) { cerr << myrank << " detected overflow in totsend calculation, line" << __LINE__ << endl; }
	int64_t totrecv = accumulate(recvcnt, recvcnt+nprocs, static_cast<int64_t>(0));
	if (totrecv < 0) { cerr << myrank << " detected overflow in totrecv calculation, line" << __LINE__ << endl; }
	DBG("totsend=%lld totrecv=%lld\n", (lld) totsend, (lld) totrecv);
       
    int p_buff_len = ((n_kmers * 2) + nprocs - 1)/nprocs;

    for (int i=0; i < nprocs; i++) {
        sdispls[i] = i * p_buff_len;
        rdispls[i] = i * p_buff_len;
        // printf("GPU COMM: rank %d: %d %d %d %d \n", myrank, sendcnt[i], sdispls[i], recvcnt[i], rdispls[i] );
    }
    // int *d_recvbuf;
    // checkCuda (cudaMalloc(&d_recvbuf, n_kmers * 2 * sizeof(uint64_t*)), __LINE__);

	double exch_time = MPI_Wtime();

	CHECK_MPI( MPI_Alltoallv(outgoing, sendcnt, sdispls, MPI_UINT64_T, recvbuf, recvcnt, rdispls, MPI_UINT64_T, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Alltoallv(len_smers, sendcnt, sdispls, MPI_UNSIGNED_CHAR, recvbuf_len, recvcnt, rdispls, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD) );
	
	exch_time = MPI_Wtime() - exch_time;
	tot_exch_time_smer += exch_time;

	// cout << "Smer exch- totsend: " << totsend << ", exch (smers+len) time: " 
	// << exch_time << ", total: " << tot_exch_time_smer << endl;
	// CHECK_MPI( MPI_Alltoallv(d_outgoing, sendcnt, sdispls, MPI_LONG, d_recvbuf, recvcnt, rdispls, MPI_LONG, MPI_COMM_WORLD) );
	// cout << "CPU alltoallv() + cudaMemcpy(future): " <<  MPI_Wtime() - exch_time1 << endl;
	// cout << "GPU direct alltoallv(): " <<  MPI_Wtime() - exch_time1 << endl;

	// /******* Performance reporting *******/
	// performance_report_time = MPI_Wtime();

	// double global_min_time = 0.0;
	// CHECK_MPI( MPI_Reduce(&exch_time, &global_min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	// double global_max_time = 0.0;
	// CHECK_MPI( MPI_Reduce(&exch_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );

	// serial_printf("KmerMatch:%s exchange iteration %d pass %d: sent min %lld bytes, sent max %lld bytes, recv min %lld bytes, recv max %lld bytes, in min %.3f s, max %.3f s\n",
	// 	__FUNCTION__, exchange_iter, pass, global_mins[SND], global_maxs[SND], global_mins[RCV], global_maxs[RCV], global_min_time, global_max_time);
	// performance_report_time = MPI_Wtime()-performance_report_time;
	/*************************************/
   
    // checkCuda (cudaMemcpy(recvbuf, d_recvbuf, n_kmers * 2 * sizeof(uint64_t), cudaMemcpyDeviceToHost), __LINE__); 

	DBG("DeleteAll: recvcount=%lld, sendct=%lld\n", (lld) recvcnt, (lld) sendcnt);
	// DeleteAll(rdispls, sdispls, recvcnt, sendcnt);
	delete(rdispls); delete(sdispls); 
	// c_smers.clear();
    //serial_printf("exchanged totsend=%lld, totrecv=%lld, pass=%d\n", (lld) totsend, (lld) totrecv, pass);
    // exchange_iter++;
    tot_exch_time=MPI_Wtime()-tot_exch_time-performance_report_time;
	MPI_Pcontrol(-1,"Exchange");
	return recvbuf;
}


double tot_exch_time_smer_concat = 0; 
void Exchange_concat_supermers(uint64_t* outgoing, unsigned char* len_smers,
	uint64_t *recvbuf, unsigned char *recvbuf_len, int *sendcnt, int *sendcnt_nsmers, int *recvcnt, 
	int *recvcnt_nsmers, int n_kmers)
{

	MPI_Pcontrol(1,"Exchange");
	double tot_exch_time = MPI_Wtime();
	double performance_report_time = 0.0; 

	int * sdispls = new int[nprocs];
	int * rdispls = new int[nprocs];
    
    CHECK_MPI( MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD) );  // share the request counts
    CHECK_MPI( MPI_Alltoall(sendcnt_nsmers, 1, MPI_INT, recvcnt_nsmers, 1, MPI_INT, MPI_COMM_WORLD) );  // share the request counts
    
   
    int p_buff_len = (n_kmers + nprocs - 1)/nprocs * 4;
    int sd = 0, rv =0;
    for (int i=0; i < nprocs; i++) {
        sdispls[i] = i * p_buff_len;
        rdispls[i] = i * p_buff_len;
    }
	
	int64_t totsend = accumulate(sendcnt, sendcnt+nprocs, static_cast<int64_t>(0));
	if (totsend < 0) { cerr << myrank << " detected overflow in totsend calculation, line" << __LINE__ << endl; }
	int64_t totrecv = accumulate(recvcnt, recvcnt+nprocs, static_cast<int64_t>(0));
	if (totrecv < 0) { cerr << myrank << " detected overflow in totrecv calculation, line" << __LINE__ << endl; }
	DBG("totsend=%lld totrecv=%lld\n", (lld) totsend, (lld) totrecv);

 	double exch_time = MPI_Wtime();

	CHECK_MPI( MPI_Alltoallv(outgoing, sendcnt, sdispls, MPI_UINT64_T, recvbuf, recvcnt, rdispls, MPI_UINT64_T, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Alltoallv(len_smers, sendcnt_nsmers, sdispls, MPI_UNSIGNED_CHAR, recvbuf_len, recvcnt_nsmers, rdispls, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD) );
	
	exch_time = MPI_Wtime() - exch_time;
	tot_exch_time_smer_concat += exch_time;

	// cout << "concated Smer exch - totsend: " << totsend << ", exch time: " 
	// << exch_time << ", total: " << tot_exch_time_smer_concat << endl;

	delete(rdispls); delete(sdispls); 

    tot_exch_time=MPI_Wtime()-tot_exch_time-performance_report_time;
	MPI_Pcontrol(-1,"Exchange");
	return;// recvbuf;
}

size_t build_concat_supermer(vector<string> seqs, size_t offset)
{
	size_t nreads = seqs.size(), max_slen = 0;
	
	size_t memoryThreshold = (MAX_ALLTOALL_MEM / nprocs) * 2;  //copied form diBella
	size_t bytesperentry = 8; // unint64 //confirm with Aydin
	int * recvcnt = new int[nprocs];
	int * sendcnt = new int[nprocs];
	int * nsmers = new int[nprocs];
	int * recvcnt_nsmers = new int[nprocs];
	memset(sendcnt, 0, sizeof(int) * nprocs);
	memset(nsmers, 0, sizeof(int) * nprocs);
	unsigned int n_kmers = 0, kmersthisbatch = 0;
	
	for(size_t i=offset; i< nreads; ++i)
	{		
		if (seqs[i].length() <= KMER_LENGTH)  // skip too short seqs
			continue;
		n_kmers += seqs[i].length() - KMER_LENGTH + 1;
	}
	nkmers_sofar += n_kmers;

	/*long 1D array. Size randomly chosen for now. Might overflow. */
	unsigned char* outgoing_lensmers = new unsigned char[n_kmers * 4]; // worst case #smers = #kmers
	int p_buff_len = (n_kmers + nprocs - 1)/nprocs * 4; //random

	// memset(supermer_bins, 0, 10 * sizeof(int));

	vector<string> long_smer(nprocs);
	for(size_t i=offset; i< nreads; ++i)
	{		
		if (seqs[i].length() <= KMER_LENGTH)  // skip too short seqs
			continue;
		kmersthisbatch += (seqs[i].length()-KMER_LENGTH+1);
		
		/* Build supermers from the current read*/	
		int order = 0, prev_order;
		string cur_kmer = seqs[i].substr(0, KMER_LENGTH);
		string minimizer = find_minimizer(cur_kmer, order);
		string prev_minimizer = minimizer;
		prev_order = order;
		int cur_slen =  KMER_LENGTH;
		
		char c_m[32]; strcpy(c_m, minimizer.c_str());
		unsigned int owner = murmur3_32(*(uint32_t *)c_m);
		outgoing_lensmers[owner * p_buff_len + nsmers[owner]++] = (unsigned char)cur_slen;
		long_smer[owner] += cur_kmer;
		// nsmers[owner]++;

		for (int c = 1; c < seqs[i].length()-KMER_LENGTH+1; ++c){ 
			
			cur_kmer = seqs[i].substr(c, KMER_LENGTH);
			order = c;
			minimizer = find_minimizer(cur_kmer, order);
			
			if(minimizer == prev_minimizer && order == prev_order) {	
				long_smer[owner] += seqs[i][c+KMER_LENGTH-1];
	    		cur_slen++;	
			}
			else {			
				outgoing_lensmers[owner * p_buff_len + nsmers[owner] - 1] = (unsigned char)cur_slen;

				cur_slen = KMER_LENGTH;
				char c_m[32]; strcpy(c_m, minimizer.c_str());
				owner = murmur3_32(*(uint32_t *)c_m);
				long_smer[owner] += cur_kmer;
				outgoing_lensmers[owner * p_buff_len + nsmers[owner]++] = (unsigned char)cur_slen;

				if( nsmers[owner] > p_buff_len)
					cout << endl << "MISSION ABORT" << endl << endl;
			}
			prev_minimizer = minimizer;
			prev_order = order;
		}
		outgoing_lensmers[owner * p_buff_len + nsmers[owner]-1] = (unsigned char)cur_slen;
	}

	size_t totsend = 0;
	for (int i = 0; i < nprocs; ++i){
		size_t arr_len = (long_smer[i].length() + 31)/32;
		sendcnt[i] = arr_len;
		totsend += arr_len;
	}

	/* compressing the every 32 char and putting in one long array */
	int loc = 0;
	uint64_t *long_csmer = new uint64_t[totsend]; 

	for (int i = 0; i < nprocs; ++i){
		size_t arr_len = (long_smer[i].length() + 31)/32;
		
		for (int j = 0; j < arr_len; j+=32) {
			uint64_t tmp = 0;//compress_smer(cur_string, 32);
			
			for (int k = 0; k < 32 && long_smer[i][j+k] != '\0'; ++k){
				char s = long_smer[i][j+k];
		    	int jj = k % 32;
		   		size_t x = ((s) & 4) >> 1;
		    	tmp |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-jj)));
			}
			long_csmer[loc++] = tmp;
		}
	}

	unsigned char *recv_slen = (unsigned char*) malloc(n_kmers * 4 * sizeof(unsigned char)); 
	uint64_t *recv_concat_smers = (uint64_t*) malloc(n_kmers * 4 * sizeof(uint64_t)); 
	
	Exchange_concat_supermers(long_csmer, outgoing_lensmers, recv_concat_smers, recv_slen, sendcnt, nsmers, recvcnt, recvcnt_nsmers, n_kmers);

	concat_parse_supermer_N_build_kmercounter(recv_concat_smers, recv_slen, recvcnt,  recvcnt_nsmers, p_buff_len);
	
	free(recv_slen); 
	free(recv_concat_smers); 
	// free(outgoing_csmers); free(outgoing_lensmers);
	delete[] sendcnt; delete[] recvcnt; delete[] nsmers; delete[] recvcnt_nsmers;
	return nreads;
}
std::unordered_map<int,uint64_t> histogram_mini; 
size_t build_supermer(vector<string> seqs, size_t offset)
{
	size_t nreads = seqs.size(), max_slen = 0;
	vector<string> supermers;
	vector<uint64_t> c_supermers; //compressed supermers
	vector<unsigned char> len_smers; //length of c_supermers 
	size_t memoryThreshold = (MAX_ALLTOALL_MEM / nprocs) * 2;  //copied form diBella
	size_t bytesperentry = 8; // unint64 //confirm with Aydin
	int * recvcnt = new int[nprocs];
	int * sendcnt = new int[nprocs];
	memset(sendcnt, 0, sizeof(int) * nprocs);
	unsigned int n_kmers = 0, kmersthisbatch = 0;
	
	for(size_t i=offset; i< nreads; ++i)
	{		
		if (seqs[i].length() <= KMER_LENGTH)  // skip too short seqs
			continue;
		n_kmers += seqs[i].length() - KMER_LENGTH + 1;
	}
	nkmers_sofar += n_kmers;

	//long 1D array. Size randomly chosen for now. Might overflow.
	uint64_t* outgoing_csmers = new uint64_t[n_kmers * 2]; // worst case #smers = #kmers
	unsigned char* outgoing_lensmers = new unsigned char[n_kmers * 2]; // worst case #smers = #kmers
	int p_buff_len = ((n_kmers * 2) + nprocs - 1)/nprocs;

	// memset(supermer_bins, 0, 10 * sizeof(int));
	int max_hit = 0;
	for(size_t i=offset; i< nreads; ++i)
	{		
		if (seqs[i].length() <= KMER_LENGTH) { // skip too short seqs
			continue;
		}

		kmersthisbatch += (seqs[i].length()-KMER_LENGTH+1);
		
		// Build supermers from the current read	
		
		unsigned int owner;
		int order = 0, prev_order;
		string cur_kmer = seqs[i].substr(0, KMER_LENGTH);
		string minimizer = find_minimizer(cur_kmer, order);
		string prev_minimizer = minimizer;
		prev_order = order;
		int cur_slen =  KMER_LENGTH;
		uint64_t compressed_supermer = compress_kmer(cur_kmer); 
		
		char c_m[32];
		strcpy(c_m, minimizer.c_str());
		owner = murmur3_32(*(uint32_t *)c_m);

		supermers.push_back(cur_kmer);
		c_supermers.push_back(compressed_supermer);
		len_smers.push_back((unsigned char)cur_slen); 		// unsigned char cc = cur_slen;
		outgoing_csmers[owner * p_buff_len + sendcnt[owner]] = compressed_supermer;
		outgoing_lensmers[owner * p_buff_len + sendcnt[owner]++] = (unsigned char)cur_slen;
		
		for (int c = 1; c < seqs[i].length()-KMER_LENGTH+1; ++c)
		{ 
			cur_kmer = seqs[i].substr(c, KMER_LENGTH);
			order = c;
			minimizer = find_minimizer(cur_kmer, order);

			if(minimizer == prev_minimizer && order == prev_order) {
							
				char s = seqs[i][c+KMER_LENGTH-1];
				supermers[supermers.size()-1] += s;			
				int k = cur_slen;
	    		int j = k % 32;
	   			size_t x = ((s) & 4) >> 1;
	    		compressed_supermer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j)));
	    		cur_slen++;
	   			// c_supermers[supermers.size()-1] = compressed_supermer;
				// len_smers[len_smers.size()-1] = (unsigned char)cur_slen;		
			}
			else {
				// this needs woek...5 diff in total kmer count from HT..probably form the very first kmer
				c_supermers[supermers.size()-1] = compressed_supermer;
				len_smers[len_smers.size()-1] = (unsigned char)cur_slen;	
				outgoing_csmers[owner * p_buff_len + sendcnt[owner] - 1] = compressed_supermer;
				outgoing_lensmers[owner * p_buff_len + sendcnt[owner] - 1] = (unsigned char)cur_slen;

				// new supermer starts
				compressed_supermer = compress_kmer(cur_kmer);
				cur_slen = KMER_LENGTH;
				supermers.push_back(cur_kmer);
				c_supermers.push_back(compressed_supermer);	
				len_smers.push_back((unsigned char)KMER_LENGTH);	
				char c_m[32];
				strcpy(c_m, minimizer.c_str());
				owner = murmur3_32(*(uint32_t *)c_m);
				outgoing_csmers[owner * p_buff_len + sendcnt[owner]] = compressed_supermer;
				outgoing_lensmers[owner * p_buff_len + sendcnt[owner]] = (unsigned char)cur_slen;
				sendcnt[owner]++;

				int hist = murmur3_32_tmp(*(uint32_t *)c_m);
				auto found = histogram_mini.find(hist);// == kcounter_cpu_concat; .end() )
			
				if(found != histogram_mini.end()){
				   	found->second += 1;
				   	max_hit = max(max_hit, (int)found->second);
				}
				else histogram_mini.insert({hist,1}); 	 
					
			}
			prev_minimizer = minimizer;
			prev_order = order;
		}
		c_supermers[supermers.size()-1] = compressed_supermer;
		len_smers[len_smers.size()-1] = (unsigned char)cur_slen;	
		outgoing_csmers[owner * p_buff_len + sendcnt[owner]-1] = compressed_supermer;
		outgoing_lensmers[owner * p_buff_len + sendcnt[owner]-1] = (unsigned char)cur_slen;
		// sendcnt[owner]++;

		int maxsending = 0;
		for (int p = 0; p < nprocs; ++p)
			maxsending = max(sendcnt[p], maxsending);
		
		if (maxsending * bytesperentry >= memoryThreshold || (kmersthisbatch + seqs[i].length()) * bytesperentry >= MAX_ALLTOALL_MEM) { 
			nreads = i+1; // start with next read
			break;
		}
	}
    uint64_t totkmer = 0, nHTsize= 0;
    max_hit = max_hit/20;
    int histbin[max_hit];
    memset(&histbin, 0, max_hit * sizeof(int));
    size_t tot_hist = 0;
	for ( auto it = histogram_mini.begin(); it!= histogram_mini.end(); ++it ){
		tot_hist+=it->second;

		cout << it->second << " ";
	}
	cout  << "\ntotl ele in hist " << tot_hist << endl;
	

	
	for (int p = 0; p < nprocs; ++p)
		cout << sendcnt[p] << " ";
	cout << endl;

	// Supermer stats
	total_supermers += supermers.size();
	// cout << "#supermers: " << c_supermers.size() << " " << supermers.size() << endl;

	unsigned char *recv_slen = (unsigned char*) malloc(n_kmers * 2 * sizeof(unsigned char)); 
	uint64_t *recv_smers = (uint64_t*) malloc(n_kmers * 2 * sizeof(uint64_t)); 
	
	Exchange_supermers(outgoing_csmers, outgoing_lensmers, recv_smers, recv_slen, sendcnt, recvcnt, n_kmers);

	// cout << "After exchange: " << endl;
	int totssmer = 0, totrsmer = 0;
	for (int i = 0; i < nprocs; ++i)
	{
		totssmer += sendcnt[i];
		totrsmer += recvcnt[i];
		// cout << "from main(): proc " << myrank << ": #sendcnt  " << sendcnt[i] << ", #recvcnt: " << recvcnt[i]  << "; ";
	}
	// cout << "\nTotal send-recv count " << totssmer << " " << totrsmer << endl;
	parse_supermer_N_build_kmercounter(recv_smers, recv_slen, recvcnt, p_buff_len);
	// parse_supermer_N_build_kmercounter(c_supermers, len_smers, recvcnt, p_buff_len);


	/*********Correctness check*******/
	for (int j = 0; j < supermers.size(); ++j)
	{
		string d_kmer = decompress_smer(c_supermers[j], len_smers[j]);
			   		
		// if(supermers[j] != d_kmer)
		// 	cout << j << ": Didnt match " << supermers[j] << " " << d_kmer << endl;
		// else
		// 	cout << j << ": Matched " << supermers[j] << " " << d_kmer << endl;
		// cout << i << ": "  << supermers[i] <<" " << c_supermers[i] << " " << decompress_smer(c_supermers[i], len_smers[i]) << " " << len_smers[i] << endl; 
		// size_t slen = supermers[j].length();
		// max_slen = std::max(max_slen, slen);
		// supermer_bins[(slen - KMER_LENGTH)]++;
		// tot_char += slen;
	}
	// supermers.clear();
 //    std::cout << std::setw(30);
	// cout << "#supermers: " << total_supermers <<", Total length of all supermers: " 
	// 	<< tot_char << ", #kmers: " << total_kmers <<", len of longest supermer: " << max_slen << endl;
	// for (int j = 0; j < 20; ++j)
	// 	cout << "#supermers between length " << j+KMER_LENGTH << " - " << j+1+KMER_LENGTH << " " << supermer_bins[j] << endl;
	free(recv_slen); free(recv_smers); free(outgoing_csmers); free(outgoing_lensmers);
	delete[] sendcnt; delete[] recvcnt;
	return nreads;
}

void getKmers_noncompress(vector<string> seqs, int offset) {
   
	std::unordered_map<std::string,uint64_t> kcounter_cpu_kmer;
	size_t nreads = seqs.size();
	for(size_t i=offset; i< nreads; ++i){
		
		if (seqs[i].length() <= KMER_LENGTH) continue;
		unsigned int seq_len = seqs[i].length();
   		
   		for(int j = 0; j < seq_len - KMER_LENGTH + 1; ++j) {

			string cur_kmer = seqs[i].substr(j, KMER_LENGTH);
	    	auto found = kcounter_cpu_kmer.find(cur_kmer);// == kcounter_cpu_kmer.end() )
	    	if(found != kcounter_cpu_kmer.end())
	    		found->second += 1;
	    	else kcounter_cpu_kmer.insert({cur_kmer,1}); 
		}
	}
	uint64_t totkmer = 0, nHTsize= 0;
	for ( auto it = kcounter_cpu_kmer.begin(); it!= kcounter_cpu_kmer.end(); ++it )
	{
		// if(it->second > 1)
		{
		nHTsize++;
		totkmer += it->second;
		}
	}

    std::cout << "\ncounter from regular kmers - HTsize: " << nHTsize 
    	<< " totkmers from HT: " << totkmer << std::endl;
    return;// kmers;
}
