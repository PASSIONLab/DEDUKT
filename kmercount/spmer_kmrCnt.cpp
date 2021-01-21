#include <cstdint>
#include <stdint.h>
#include <limits>
#include <iostream>
#include "Kmer.hpp"
#include "Friends.h"
#include "MPIType.h"
#include "SimpleCount.h"
#include "FriendsMPI.h"
#include "Pack.h"
#include "spmer_kmrCnt.h"

using namespace std;

// 32 bit Murmur3 hash
int murmur3_32( int k )
{
	k ^= k >> 16;
	k *= 0x85ebca6b;
	k ^= k >> 13;
	k *= 0xc2b2ae35;
	k ^= k >> 16;

	return k;
}

uint64_t murmur3_64( uint64_t k )
{
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xff51afd7ed558ccd);
	k ^= k >> 33;
	k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
	k ^= k >> 33;

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

// uint64_t compress_kmer(string cur_kmer){

//         uint64_t c_kmer = 0, cc_kmer = 0;
//         for (int k = 0; k < KMER_LENGTH; ++k)
//         {
//                 char s = cur_kmer[k];
//                 int j = k % 32;
//                 size_t x = (s & 4) >> 1;
//                 c_kmer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j)));
//                 // switch(s) { //redefined
//   //                    case 'A': c_kmer |= ((x + (x^1)) << (2*(31-j)));  break;
//   //                    case 'C': c_kmer |= ((x + (x^0)) << (2*(31-j)));  break;
//   //                    case 'G': c_kmer |= ((x + (x^3)) << (2*(31-j)));  break;
//   //                    case 'T': c_kmer |= ((x + (x^2)) << (2*(31-j)));  break;
//   //            }
//         }
//         // const uint64_t mask = ((int64_t) 0x3);
//  //     uint64_t endmask = 0;
//  //     if (KMER_LENGTH % 32) {
//  //       endmask = (((uint64_t) 2) << (2*(31-(KMER_LENGTH%32)) + 1)) - 1;
//  //       // k == 0 :                0x0000000000000000
//  //       // k == 1 : 2 << 61  - 1 : 0x3FFFFFFFFFFFFFFF
//  //       // k == 31: 2 << 1   - 1 : 0x0000000000000003
//  //     }
//  //     endmask = ~endmask;
//  //     c_kmer &= endmask;
//         return c_kmer;
// }


uint64_t find_minimizer(uint64_t kmer, int mlen){

	int klen = KMER_LENGTH;
	keyType minimizer = std::numeric_limits<uint64_t>::max(); 
	keyType mask = pow(2, 2*mlen) - 1;

	for (int m = 0; m < (klen - mlen ); ++m){
		keyType mmer =  (kmer >> (2*(31-(mlen+m -1)))) & mask;

		if( mmer < minimizer ) 
			minimizer = mmer;
	}
	return minimizer;
}
// obsolete
string find_minimizer(string kmer, int mlen){

	string minimizer = "ZZZZZZZ";

	for (int m = 0; m < (KMER_LENGTH - mlen); ++m)
	{
		string mmer = kmer.substr(m, mlen);
		if( mmer < minimizer ) {
			minimizer = mmer;
		}
	}
	return minimizer;
}

size_t total_kmers = 0, total_supermers = 0, tot_char = 0;
size_t nkmers_sofar = 0;

std::unordered_map<uint64_t,uint64_t> smer_kcounter; 
void parse_supermer_N_build_kmercounter(uint64_t* recvbuf, unsigned char* len_smers, int* recvcnt, int p_buff_len){

	std::vector<uint64_t> kmers_compressed;
	int klen = KMER_LENGTH;
	uint64_t mask = pow(2, 2*klen) - 1;
	for(uint64_t i= 0; i < nprocs ; ++i) 
	{
		for(uint64_t j= 0; j <  recvcnt[i] ; ++j)
		{ 
			uint64_t cur_csmer = recvbuf[i * p_buff_len + j];
			unsigned char c = len_smers[i * p_buff_len + j];
			int slen = (int)c;
			string d_smer = decompress_smer(cur_csmer, c);
			int nkmer = slen - KMER_LENGTH + 1;
			for(int k = 0; k < slen - KMER_LENGTH + 1; ++k) {
				uint64_t cur_ckmer =  (cur_csmer >> (2*(31-(klen+k -1)))) & mask; //Assume LSB
				auto found = smer_kcounter.find(cur_ckmer);
				if(found != smer_kcounter.end())
					found->second += 1;
				else smer_kcounter.insert({cur_ckmer,1}); 
			}
		}
	}

	uint64_t totalPairs = 0, HTsize= 0, totalPairs_num = 0, HTsize_num= 0;
	for ( auto it = smer_kcounter.begin(); it!= smer_kcounter.end(); ++it ){
		HTsize++; 
		totalPairs += it->second;
	}

	size_t allrank_hashsize = 0, allrank_totalPairs = 0;
	CHECK_MPI( MPI_Reduce(&HTsize, &allrank_hashsize,  1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&totalPairs, &allrank_totalPairs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );

	size_t allrank_kmersthisbatch = 0;
	size_t allrank_kmersprocessed = 0;
	// CHECK_MPI( MPI_Reduce(&nkmers_thisBatch, &allrank_kmersthisbatch, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&nkmers_sofar, &allrank_kmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );

	// std::cout << "rank: " << myrank << ", Smer - CPU HTsize: " << HTsize 
	// 	<< " #kmers from HT: " << totalPairs << ", ideal #kmers: " << allrank_kmersprocessed << std::endl;

	if(myrank == 0)
		std::cout << "Smer - CPU HTsize: " << allrank_hashsize 
			<< " #kmers from HT: " << allrank_totalPairs << ", ideal #kmers: " << allrank_kmersprocessed << std::endl;

	return;
}

double tot_exch_time_smer = 0; 
uint64_t* Exchange_supermers(uint64_t* outgoing, unsigned char* len_smers, 
		uint64_t* recvbuf, unsigned char* recvbuf_len, int *sendcnt, int *recvcnt, int n_kmers){

	double tot_exch_time = MPI_Wtime();
	
	int * sdispls = new int[nprocs];
	int * rdispls = new int[nprocs];

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
	}

	double exch_time = MPI_Wtime();

	CHECK_MPI( MPI_Alltoallv(outgoing, sendcnt, sdispls, MPI_UINT64_T, recvbuf, recvcnt, rdispls, MPI_UINT64_T, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Alltoallv(len_smers, sendcnt, sdispls, MPI_UNSIGNED_CHAR, recvbuf_len, recvcnt, rdispls, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD) );

	exch_time = MPI_Wtime() - exch_time;
	tot_exch_time_smer += exch_time;

	double performance_report_time = 0;//perf_reporting(exch_time, totsend, totrecv);
	
	delete(rdispls); delete(sdispls); 

	tot_exch_time = MPI_Wtime() - tot_exch_time - performance_report_time;

	return recvbuf;
}

/**
 * @param outgoing_csmers: outgoing mpi buffer to store compressed supermers
 * @param outgoing_lensmers: outgoing mpi buffer to store the lengths of the 
 compressed supermers
 * @param sendcnt: size of outgoing mpi buffers 
 */

size_t build_supermer(vector<string> seqs, size_t offset, size_t endoffset, uint64_t* outgoing_csmers,
		unsigned char* outgoing_lensmers, int *sendcnt, int mlen)
{
	int klen = KMER_LENGTH;
	if(myrank == 0) cout << "FIXIT " << endl;
	size_t nreads = endoffset;// seqs.size(), 
	vector<string> supermers;
	vector<unsigned char> len_smers; //length of c_supermers 
	size_t memoryThreshold = (MAX_ALLTOALL_MEM / nprocs) * 2;  //copied form diBella
	size_t bytesperentry = 8; // unint64 //confirm with Aydin

	memset(sendcnt, 0, sizeof(int) * nprocs);
	unsigned int n_kmers = 0, kmersthisbatch = 0;

	for(size_t i=offset; i< nreads; ++i){		

		if (seqs[i].length() <= KMER_LENGTH)  // skip too short seqs
			continue;
		n_kmers += seqs[i].length() - KMER_LENGTH + 1;
	}

	int p_buff_len = ((n_kmers * 2) + nprocs - 1)/nprocs;

	for(size_t i=offset; i< nreads; ++i){	

		if (seqs[i].length() <= KMER_LENGTH) continue; //skip short seqs

		kmersthisbatch += (seqs[i].length()-KMER_LENGTH+1);
		int cur_nkmer = seqs[i].length()-KMER_LENGTH+1;

		//* Build supermers from the current read */
		unsigned int owner;
		int order = 0, prev_order;
		int window = 32 - KMER_LENGTH;

		for (int c = 0; c < cur_nkmer; c+=window)
		{
			string cur_kmer = seqs[i].substr(c, KMER_LENGTH);
			uint64_t compressed_supermer = compress_kmer(cur_kmer);
			// string minimizer = find_minimizer(cur_kmer, order); 
			// string prev_minimizer = minimizer;
			// char c_m[32]; strcpy(c_m, minimizer.c_str());
			// owner = (murmur3_64(*(uint64_t *)c_m)) & (nprocs - 1);

			uint64_t c_minimizer = find_minimizer(compressed_supermer, mlen); 
			uint64_t c_prev_minimizer = c_minimizer;

			uint64_t hsh = murmur3_64(c_minimizer);
			double range = static_cast<double>(hsh) * static_cast<double>(nprocs);
			owner = range / static_cast<double>(numeric_limits<uint64_t>::max());

			int cur_slen = KMER_LENGTH;
			supermers.push_back(cur_kmer); //sanity check
			len_smers.push_back((unsigned char)cur_slen); 		
			outgoing_csmers[owner * p_buff_len + sendcnt[owner]] = compressed_supermer;
			outgoing_lensmers[owner * p_buff_len + sendcnt[owner]++] = (unsigned char)cur_slen;

			for (int w = 1; w < window && (c+w) < cur_nkmer; ++w){

				cur_kmer = seqs[i].substr(c+w, KMER_LENGTH);
				c_minimizer = find_minimizer(compress_kmer(cur_kmer), mlen);

				if(c_minimizer == c_prev_minimizer ){

					char s = seqs[i][c+w+KMER_LENGTH-1];
					supermers[supermers.size()-1] += s;		
					int k = cur_slen;
					int j = k % 32;
					size_t x = ((s) & 4) >> 1;
					compressed_supermer |= ((x + ((x ^ (s & 2)) >>1)) << (2*(31-j)));
					cur_slen++;	
				}
				else {
					//* update last supermer */
					len_smers[len_smers.size()-1] = (unsigned char)cur_slen;	
					outgoing_csmers[owner * p_buff_len + sendcnt[owner] - 1] = compressed_supermer;
					outgoing_lensmers[owner * p_buff_len + sendcnt[owner] - 1] = (unsigned char)cur_slen;

					//* new supermer starts */
					compressed_supermer = compress_kmer(cur_kmer);
					cur_slen = KMER_LENGTH;
					supermers.push_back(cur_kmer);

					len_smers.push_back((unsigned char)KMER_LENGTH);	
					hsh = murmur3_64(c_minimizer);// & (nprocs - 1);
					range = static_cast<double>(hsh) * static_cast<double>(nprocs);
					owner = range / static_cast<double>(numeric_limits<uint64_t>::max());
					outgoing_csmers[owner * p_buff_len + sendcnt[owner]] = compressed_supermer;
					outgoing_lensmers[owner * p_buff_len + sendcnt[owner]] = (unsigned char)cur_slen;
					sendcnt[owner]++;	 
				}
				// prev_minimizer = minimizer;
				c_prev_minimizer = c_minimizer;
			}
			len_smers[len_smers.size()-1] = (unsigned char)cur_slen;	
			outgoing_csmers[owner * p_buff_len + sendcnt[owner]-1] = compressed_supermer;
			outgoing_lensmers[owner * p_buff_len + sendcnt[owner]-1] = (unsigned char)cur_slen;
		}

		int maxsending = 0;
		for (int p = 0; p < nprocs; ++p)
			maxsending = max(sendcnt[p], maxsending);

		if (maxsending * bytesperentry >= memoryThreshold || (kmersthisbatch + seqs[i].length()) * bytesperentry >= MAX_ALLTOALL_MEM) { 
			nreads = i+1; // start with next read
			break;
		}
	}

	// uint64_t totkmer = 0, nHTsize= 0;

	// int totssmer = 0;
	// for (int i = 0; i < nprocs; ++i)
	// 	totssmer += sendcnt[i];

	// cout << myrank << " CPU totsmer: " <<  totssmer << ", smer distribution: avg: " << totssmer/nprocs << "; ";
	// for (int p = 0; p < nprocs; ++p)
	// 	cout << p << ": " << sendcnt[p] << ", ";
	// cout << endl;

	return nreads;
}


size_t spmer_kmerCount(vector<string> seqs, size_t offset, size_t endoffset, int klen, int mlen)
{
	KMER_LENGTH = klen;
	Kmer::set_k(klen);
	// cout << "\n\n lengths " << KMER_LENGTH << " " << mlen;
	if(myrank == 0) cout << "FIXIT " << endl;
	size_t nreads = endoffset;// seqs.size(), 
	size_t max_slen = 0;
	vector<string> supermers;
	// vector<uint64_t> c_supermers; //compressed supermers
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

	//***** Build Supermers on GPU *****
	build_supermer(seqs, offset, endoffset, outgoing_csmers, outgoing_lensmers, sendcnt, mlen);

	uint64_t totkmer = 0, nHTsize= 0;

	int totssmer = 0;
	for (int i = 0; i < nprocs; ++i)
		totssmer += sendcnt[i];

	// cout << "CPU totsmer: " <<  totssmer << ", smer distribution: avg: " << totssmer/nprocs << "; ";
	// for (int p = 0; p < nprocs; ++p)
	// 	cout << p << ": " << sendcnt[p] << ", ";
	// cout << endl;

	unsigned char *recv_slen = (unsigned char*) malloc(n_kmers * 2 * sizeof(unsigned char)); 
	uint64_t *recv_smers = (uint64_t*) malloc(n_kmers * 2 * sizeof(uint64_t)); 

	//***** Exchange supermers on CPU *****
	Exchange_supermers(outgoing_csmers, outgoing_lensmers, recv_smers, recv_slen, sendcnt, recvcnt, n_kmers);

	int totrsmer = 0;
	for (int i = 0; i < nprocs; ++i)
		totrsmer += recvcnt[i];

	//***** Parse supermers and build kcounter on CPU *****	
	parse_supermer_N_build_kmercounter(recv_smers, recv_slen, recvcnt, p_buff_len);

	/*********Correctness check*******/
	/*	for (int j = 0; j < supermers.size(); ++j)
		{
		string d_kmer = decompress_smer(c_supermers[j], len_smers[j]);		   		
		if(supermers[j] != d_kmer)
		cout << j << ": Didnt match " << supermers[j] << " " << d_kmer << endl;
		else
		cout << j << ": Matched " << supermers[j] << " " << d_kmer << endl;
		}
		*/
	free(recv_slen); free(recv_smers); 
	free(outgoing_csmers); free(outgoing_lensmers);
	delete[] sendcnt; delete[] recvcnt;
	return nreads;
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

size_t build_concat_supermer(vector<string> seqs, size_t offset, int mlen)
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
		string minimizer = find_minimizer(cur_kmer, mlen);
		string prev_minimizer = minimizer;
		prev_order = order;
		int cur_slen =  KMER_LENGTH;

		char c_m[32]; strcpy(c_m, minimizer.c_str());
		unsigned int owner = (murmur3_64(*(uint64_t *)c_m)) & (nprocs - 1);
		outgoing_lensmers[owner * p_buff_len + nsmers[owner]++] = (unsigned char)cur_slen;
		long_smer[owner] += cur_kmer;
		// nsmers[owner]++;

		for (int c = 1; c < seqs[i].length()-KMER_LENGTH+1; ++c){ 

			cur_kmer = seqs[i].substr(c, KMER_LENGTH);
			order = c;
			minimizer = find_minimizer(cur_kmer, mlen);

			if(minimizer == prev_minimizer && order == prev_order) {	
				long_smer[owner] += seqs[i][c+KMER_LENGTH-1];
				cur_slen++;	
			}
			else {			
				outgoing_lensmers[owner * p_buff_len + nsmers[owner] - 1] = (unsigned char)cur_slen;

				cur_slen = KMER_LENGTH;
				char c_m[32]; strcpy(c_m, minimizer.c_str());
				owner = (murmur3_64(*(uint64_t *)c_m)) & (nprocs - 1);
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
