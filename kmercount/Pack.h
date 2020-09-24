#ifndef _PACK_H_
#define _PACK_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <vector>
#include <sstream>
#include <limits>
#include <array>
#include <map>
#include <tuple>
#include <locale>
#include "Kmer.hpp"
#include "Friends.h"
#include "FriendsMPI.h"

using namespace std;

extern int nprocs;
extern int myrank;

#ifdef HEAVYHITTERS
//#ifndef MAXHITTERS
//#define MAXHITTERS 200000
//#endif
extern SimpleCount<Kmer, KmerHash> *heavyhitters;
extern  UFX2ReduceObj * Frequents;
#endif


// The bloom filter pass; extensions are ignored
inline size_t FinishPackPass1(vector< vector<Kmer> > & outgoing, Kmer & kmerreal)
{
    uint64_t myhash = kmerreal.hash();  // whichever one is the representative
    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());
    
#ifdef HEAVYHITTERS
    if (heavyhitters && heavyhitters->IsMember(kmerreal)) {
        assert( heavyhitters->FindIndex(kmerreal) < heavyhitters->maxsize );
        // no-op here
        // cout << kmerreal.toString() << " is high freq" << endl;
    } else {
            assert( !heavyhitters || heavyhitters->FindIndex(kmerreal) >= heavyhitters->maxsize );
#endif
            outgoing[owner].push_back(kmerreal);
#ifdef HEAVYHITTERS
    }
#endif
    return outgoing[owner].size();
}


// The hash table pass; extensions are important
inline size_t FinishPackPass2(vector< vector<Kmer> > & outgoing, vector < vector< ReadId > > & readids, vector < vector< PosInRead > > & positions, vector<vector<array<char,2>>> & extquals,
                      vector<vector<array<char,2>>> & extseqs, Kmer & kmerreal, ReadId readId, PosInRead pos)
{
    assert( kmerreal == kmerreal.rep() );
    uint64_t myhash = kmerreal.hash();  // whichever one is the representative
    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());

    size_t location = 0, maxsize = 0;
#ifdef HEAVYHITTERS
    if (heavyhitters) {
        maxsize = heavyhitters->maxsize;
        location = heavyhitters->FindIndex(kmerreal);    // maxsize if not found
    }
#endif
    if(location < maxsize) { // in words, if kmerreal is a heavy hitter...
#ifdef HEAVYHITTERS
        assert( maxsize > 0 );
        assert( heavyhitters );
        assert( heavyhitters->IsMember( kmerreal ) );
#endif
    } else {
        // Count here
#ifdef HEAVYHITTERS
        assert( !heavyhitters || ! heavyhitters->IsMember( kmerreal ) );
#endif
        //cout << myrank << ": " << kmerreal.toString() << " is NOT a heavy hitter " << endl;
        outgoing[owner].push_back(kmerreal);
        readids[owner].push_back(readId);
        positions[owner].push_back(pos);
    }
    return outgoing[owner].size();
}

inline size_t PackEndsKmer(string & seq, string & qual, int j, Kmer &kmerreal, ReadId readid, PosInRead pos, vector< vector<Kmer> > & outgoing,
		vector< vector<ReadId> > & readids, vector< vector<PosInRead> > & positions, vector<vector<array<char,2>>> & extquals,
        vector<vector<array<char,2>>> & extseqs, int pass, int lastCountedBase, int kmer_length)
{
    bool isCounted = lastCountedBase >= j + kmer_length;
    size_t procSendCount;
    assert( seq.size() >= j + kmer_length);
    assert( seq.substr(j, kmer_length).find('N') == std::string::npos );
    assert( kmerreal == Kmer(seq.c_str() + j) );
    if (pass == 1)
    {
        if (!isCounted) return 0;
        kmerreal = kmerreal.rep();
        procSendCount = FinishPackPass1(outgoing, kmerreal);
    }
    else if (pass == 2)   // otherwise we don't care about the extensions
    {
        Kmer kmertwin = kmerreal.twin();
        if(kmertwin < kmerreal)	// the real k-mer is not lexicographically smaller
        {
            kmerreal = kmertwin;
        }
        procSendCount = FinishPackPass2(outgoing, readids, positions, extquals, extseqs, kmerreal, readid, pos);
    }
    return procSendCount;
}

/* unused in longread version
size_t PackEnds(string & seq, string & qual, int j, vector< vector<Kmer> > & outgoing, vector<vector<array<char,2>>> & extquals,
                vector<vector<array<char,2>>> & extseqs, int pass, int lastCountedBase, int kmer_length)
{
    bool isCounted = lastCountedBase >= j + kmer_length;
    if (pass == 1 && !isCounted) return 0;

    string kmerextstr;
    try {
        assert(seq.size() >= j+kmer_length);
        kmerextstr = seq.substr(j, kmer_length);
    } catch (std::exception const &exc) {
        std::cerr << "Exception caught in file " << __FILE__<< " at " << __LINE__ << ": " << exc.what() << "\n";
        std::cerr << "j = " << j << ", kmer_length = " << kmer_length << ", len seq = " << seq.length() << "\n";
        std::cerr << "seq = " << seq << "\n";
        MPI_Finalize();
        exit(1);
    }
    for (auto & c: kmerextstr) c = toupper(c);	// convert all to uppercase
    size_t found=kmerextstr.find('N');
    if (found!=std::string::npos) return 0;	// if there is an 'N', toss it

    Kmer kmerreal(kmerextstr.c_str());
    return PackEndsKmer(seq, qual, j, kmerreal, outgoing, extquals, extseqs, pass, lastCountedBase, kmer_length);
}
*/

#endif
