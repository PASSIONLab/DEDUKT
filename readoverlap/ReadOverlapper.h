/*
 * ReadOverlapper.h
 *
 *      Author: mme
 */

#ifndef READOVERLAPPER_H_
#define READOVERLAPPER_H_

#include <iostream>
#include <cstdint> // uint64_t
#include <map> // std::map
#include <unordered_map> // std::unordered_map
#include <set> // std::set
#include <utility> // std::pair
#include <tuple> // std::tuple

#include "../common/defines.h"
#include "../kmercount/kmerhashtable.h"

using namespace std;

namespace overlap {

typedef std::set<ReadOverlapPair> ReadOverlapSet;

class ReadOverlapper {

	private:
		int myrank;
		int nprocs;
		int max_num_readids;
		ReadId* read_ranges;
		KmerCountsType* kmercounts_ptr;
		ReadOverlapSet* perproc_buffers;
		ReadOverlapMap* overlapResultMap;

		inline bool nDistant(const int distance, const PosPair& a, const PosPair& b) { \
			return ((b.first - a.first >= distance) && (b.second - a.second >= distance)); \
		}
		inline bool equalPosPairs(const PosPair a, const PosPair b) { \
			return (a.first == b.first && a.second == b.second); \
		}
		static inline bool orderBefore(PosPair a, PosPair b) { \
			return (a.first < b.first && a.second < b.second);
		}
		inline ReadOverlapPair getOrderedPair(ReadId r, ReadId s, PosInRead rp, PosInRead sp);
		void printPerProcBuffers();
		//void printMap(ReadOverlapMap& returnMap, int until);
		void extractReadPairs();
		void exchangeReadNames(set<ReadId>* reqbuffers, unordered_map<ReadId,string>& readNameMap);
		void exchangeReadPairs();
		//void sortAndRemoveDuplicatePositions();
		// @deprecated
		int outputOverlapsByIdOnly(string outfilename);
		int outputOverlapsPerThread(const char* outfilename, unordered_map<ReadId,string>& readNameMap);
		int outputReadOverlapsMpi(const char* outfilename, unordered_map<ReadId,string>& readNameMap);

	public:
		ReadOverlapper(int, ReadId* readRanges, KmerCountsType* kmercounts);
		void extractAndExchangeOverlappingReads();
		int outputOverlaps(const char* filename, bool per_thread, bool cached_io, unordered_map<ReadId,string>& readNameMap);
		void filterSeeds(int min_seed_distance, int max_num_seeds);
//		int outputAndEraseOverlaps(string outfilename, unordered_map<ReadId,string>& readNameMap);
//		void filterAndExtendSeeds(ReadOverlapMap& returnMap);
		virtual ~ReadOverlapper();

	};

} /* namespace overlap */

#endif /* READOVERLAPPER_H_ */
