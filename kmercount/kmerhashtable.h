/*
 * kmerhashtable.h
 *
 * was UFXextended.h
 *
 */

#ifndef KMERCOUNT_KMERMATCH_H_
#define KMERCOUNT_KMERMATCH_H_

#include <utility>

#include "../common/defines.h"
#include "../common/readid_records.h"
#include "Kmer.hpp"
#include "../common/VectorMap.hpp"

using namespace std;

typedef tuple<READIDS, POSITIONS, int> KmerCountType;
#ifdef KHASH
  typedef std::pair<Kmer::MERARR, KmerCountType> KmerValue;
  typedef khmap_t<Kmer::MERARR, KmerCountType > KmerCountsType;
#else
  typedef pair<Kmer::MERARR, KmerCountType> KmerValue;
/*
  #if UPC_ALLOCATOR
    #include "../common/upc_allocator.hpp"
    #define USE_UPC_ALLOCATOR_IN_MPI
    #ifdef NO_UPC_MEMORY_POOL
      typedef upc_allocator::simple_allocator< KmerValue > KmerAllocator;
    #else
      #include "../common/upc_memory_pool.hpp"
//      #pragma message("Using upc_allocator::SingletonMemoryPool for kmer hashtable allocations")
      typedef upc_allocator::SingletonMemoryPool< KmerValue > KmerAllocator;
    #endif
//  typedef U_MAP<Kmer::MERARR, KmerCountType, std::hash<Kmer::MERARR>, std::equal_to<Kmer::MERARR>, KmerAllocator > KmerCountsType;
    typedef VectorMap< Kmer::MERARR, KmerCountType, std::hash<Kmer::MERARR>, std::less<Kmer::MERARR>, std::equal_to<Kmer::MERARR>, KmerAllocator > KmerCountsType;
  #else // not UPC_ALLOCATOR
*/
//  typedef U_MAP< Kmer::MERARR, KmerCountType, std::hash<Kmer::MERARR> > KmerCountsType;
    typedef VectorMap< Kmer::MERARR, KmerCountType, std::hash<Kmer::MERARR>, std::less<Kmer::MERARR>, std::equal_to<Kmer::MERARR> > KmerCountsType;
//  typedef std::map< Kmer::MERARR, KmerCountType, std::less<Kmer::MERARR> > KmerCountsType;
/*  #endif*/
#endif



#endif /* KMERCOUNT_KMERMATCH_H_ */
