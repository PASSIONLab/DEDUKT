#ifndef _FRIENDS_H_
#define _FRIENDS_H_

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <tuple>
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "Deleter.h"
#include "../common/readid_records.h"

#include <sys/stat.h>

using namespace std;

extern int reliable_max;

//#define KMERLONGS (MAX_KMER_SIZE/32)      // 32 = numbits(uint64_t)/2-  with 2 being the number of bits needed per nucleotide
static int KMER_LENGTH = MAX_KMER_SIZE-1;
static int KMER_PACKED_LENGTH = (MAX_KMER_SIZE-1+3)/4;

extern int myrank;
typedef array<int,4> ACGT;

struct filedata
{
    char filename[MAX_FILE_PATH];
    size_t filesize;
};

inline void RevCompSwitch(char & sit)
{
    switch(sit) {
        case 'G':
            sit = 'C';
            break;
        case 'A':
            sit = 'T';
            break;
        case 'T':
            sit = 'A';
            break;
        case 'C':
            sit = 'G';
            break;
        default:
            break;
    }
}

inline void Increment1stBasedOnExt(tuple<ACGT,ACGT, int> & mercnt, char c)
{
    switch(c)
    {
        case 'A' :
            (get<0>(mercnt))[0]++;
            break;
        case 'C' :
            (get<0>(mercnt))[1]++;
            break;
        case 'G' :
            (get<0>(mercnt))[2]++;
            break;
        case 'T' :
            (get<0>(mercnt))[3]++;
            break;
	// an optional default label
        // default :
	    // this might be an 'N' so don't worry about it
            //cout << c << " doesn't make sense" << endl;
    }
}

inline void Increment2ndBasedOnExt(tuple<ACGT,ACGT, int>  & mercnt, char c)
{
    switch(c)
    {
        case 'A' :
            (get<1>(mercnt))[0]++;
            break;
        case 'C' :
            (get<1>(mercnt))[1]++;
            break;
        case 'G' :
            (get<1>(mercnt))[2]++;
            break;
        case 'T' :
            (get<1>(mercnt))[3]++;
            break;
	// an optional default label
        // default :
	    // this might be an 'N' so don't worry about it
            // cout << c << " doesn't make sense" << endl;
    }
}

inline ostream & operator<<(ostream & os, uint8_t val)
{
    return os << static_cast<int>(val);
}

struct kmerpack	// the pair<Kmer::MERARR,int> used as value_type in map is not guaranteed to be contiguous in memory
{
	Kmer::MERARR arr;
	int count;

	bool operator > (const kmerpack & rhs) const
	{ return (arr > rhs.arr); }
	bool operator < (const kmerpack & rhs) const
	{ return (arr < rhs.arr); }
	bool operator == (const kmerpack & rhs) const
	{ return (arr == rhs.arr); }
};


struct ufxpack	// 38bytes for k=51
{
	Kmer::MERARR arr;	// ~128-bits=16bytes for k=51
	int count;
	READIDS reads;
	POSITIONS positions;
	/*
	int leftmin;
	int leftmax;
	int rightmin;
	int rightmax;
	char left;
	char right;
	*/
	bool operator > (const ufxpack & rhs) const
	{ return (arr > rhs.arr); }
	bool operator < (const ufxpack & rhs) const
	{ return (arr < rhs.arr); }
	bool operator == (const ufxpack & rhs) const
	{ return (arr == rhs.arr); }
};

inline void PackIntoUFX(/*const array<int,4> & leftcnt, const array<int,4> & righcnt,*/ READIDS reads, POSITIONS positions, int count, ufxpack & pack)
{
	/*
	pair<int, char> lsort[4] = {make_pair(leftcnt[0], 'A'), make_pair(leftcnt[1], 'C'), make_pair(leftcnt[2], 'G'), make_pair(leftcnt[3], 'T')};
	pair<int, char> rsort[4] = {make_pair(righcnt[0], 'A'), make_pair(righcnt[1], 'C'), make_pair(righcnt[2], 'G'), make_pair(righcnt[3], 'T')};
	sort(lsort, lsort+4);
	sort(rsort, rsort+4);

	pack.left = lsort[3].second;	// max entry guarenteed to exist
	pack.leftmax = lsort[3].first;
	pack.leftmin = lsort[2].first;

	pack.right = rsort[3].second;
	pack.rightmax = rsort[3].first;
	pack.rightmin = rsort[2].first;
	*/

	for (int i = 0; i < reliable_max; i++) {
		pack.reads[i] = reads[i];
		pack.positions[i] = positions[i];
	}
	pack.count = count;
}

#endif
