/*
 * Aligner.h
 *
 *      Author: mellis
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <utility>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstddef>
#include <mpi.h>
#include <algorithm>

#include "../common/mpi_common.h" // include before defines.h
#include "../common/defines.h"
#include "../common/ReadIdToProcId.h"
#include "../common/FileOpenerMPI.h"

#include "../common/ReadIdToProcId.h"
#include "../common/Buffer.h" // todo add buffer object to cmake dependencies
#include "../kmercount/ParallelFASTQ.h"
#include "../kmercount/Kmer.hpp"


extern bool _skip_algnmnt_krnl;

namespace align {

	/*
	 * fields are kernel dependent and may be null.
	 */
	struct alignment_data {
		ReadId r1_id = nullReadId, r2_id = nullReadId;
		PosInRead r1_al_start = 0, r1_al_end = 0, r2_al_start = 0, r2_al_end = 0;
		/* length of actual read */
		PosInRead length_r1 = 0, length_r2 = 0;
		/* number of shared k-mers for seed-extend alignment
		 * defaults to 1, since alignment is currently calculated using 1 k-mer seed at a time
		 * (reduction can be done during post-processing)
		 */
		short similarity = 1;
		/* extension length for seed-extend alignment */
		int score = 0;

		void reset_positions(const PosInRead sH, const PosInRead eH, const PosInRead sV, const PosInRead eV)
		{
			r1_al_start = sH;
			r1_al_end = eH;
			r2_al_start = sV;
			r2_al_end = eV;
		}

		alignment_data() {}

		alignment_data(ReadId r1, ReadId r2) : r1_id(r1), r2_id(r2) {}

		alignment_data(ReadId r1, ReadId r2, PosInRead sH, PosInRead eH, PosInRead sV, PosInRead eV,
				PosInRead lenR1, PosInRead lenR2, int sc, short sim=1)
				: r1_id(r1), r2_id(r2), r1_al_start(sV), r1_al_end(eV), r2_al_start(sH), r2_al_end(eH),
				  length_r1(lenR1), length_r2(lenR2), similarity(sim), score(sc)
		{}
	};



	/*
	 * The name substitutes (stored externally) for alignment_data.r1_id and alignment_data.r2_id
	 * should be output before this function is called (unless supplied by context).
	 */
	inline std::ostream & operator<<(std::ostream & str, const alignment_data& data) {
		return (str << data.similarity << " " << data.score << " "
				<< data.r1_al_start << " " <<  data.r1_al_end << " " << data.length_r1 << " "
				<< data.r2_al_start << " " <<  data.r2_al_end << " " << data.length_r2);
	}

	class Aligner {
		typedef std::vector< ReadOverlapPair > OverlapPairsList;
		/*
		bool operator<( const std::pair<ReadId, ReadId>& pairA, const std::pair<ReadId, ReadId>& pairB) {
			return ( (pairA.first < pairB.first)
					|| (pairA.first == pairB.first && pairA.second < pairB.second) );
		}
		*/

	private:
		bool skip_kernel;
		int kmerLen;
		int maxPerPairAlignments;
		char *bella_settings_file = NULL;
		ReadOverlapMap* readPairs;
		std::vector<std::pair<std::string,std::string>> kernelOptions;
		/*
		 * readRanges[i] contains the last and greatest ReadId of processor i
		 */
		ReadId* readRangeEnds;
		ReadId globalReadIdOffset;
		std::set< ReadId >* perProcRequestIds;
		//std::set< std::pair<ReadId,ReadId> >* prevAlignmentsSet;

		/*
		 * A map of read ID's to actual sequences
		 */
		std::unordered_map<ReadId, std::string>* readMap;
		/*
		 * A map of read ID's to names/tags
		 */
		std::unordered_map<ReadId, std::string>* readNameMap;

		inline void checkNameMapping(const ReadId read, string name);
		int checkAndOrderOverlapPair(const ReadOverlapPair&);
		void checkLocality(const ReadOverlapPair& currPair, string, string, int64_t* overlaps_stats);
		bool loadPrevAlignmentSet(std::string prevalignmentsname);
		inline void bufferReadId(const ReadId readid);
		int loadOverlapsMpi(char* inputOverlaps, bool checkExistingAlignmentIds);
		void loadOverlaps(const char* overlapfile);
		void loadOverlapsPerThread(const char*);
		void storeSequences(const std::vector<std::string>&, const std::vector<std::string>& names, ReadId);
		ReadId loadSequences(const std::vector<filenamesize> & allfiles, std::vector<std::string>& allreads, std::vector<std::string>& allnames, bool cached_io, const char* base_dir);
		//inline ReadId getGlobalIndex(ReadId localIndex);
		//inline ReadId getLocalIndex(ReadId globalIndex);
		void initReadRangesAndSequences(const std::vector<filenamesize> & allfiles, bool cachedIO, const char* base_dir);
		void exchangeSequences();
		void computeLocalAlignments(ReadOverlapMap*, std::vector<alignment_data>&, std::string outfilename);

	public:
		struct scoring_info {
			int xdrop = 7;
			int match=1;
			int mismatch=-1;
			int gap=-1;
			int static_threshold = 50;
			bool aligntoend=false;
			bool adaptive=true;
			double error_rate = 0.15;
			double chernoff_delta = 0.2;
			double relaxation_margin = 300;
			double adaptive_phi;

			void set_adaptive_threshold() {
				double p_mat = (1-error_rate) * (1-error_rate);  // match
				double p_mis = 1-p_mat;				   			// mismatch/gap
				double alpha = 1;					   			// match penalty
				double beta = 1;					   			// mismatch/gap penalty
				adaptive_phi = alpha*p_mat - beta*p_mis;
			}
		} myscoring;
		Aligner(int, int, bool, char*);
		Aligner(int kmerlength, int maxNumAlignPerPair, bool skipKernel, ReadId* readRanges, ReadOverlapMap* overlapMap);
		virtual ~Aligner();
	//	void computeAlignmentsFromOverlapsInMemory(bool cachedIO, char* inputFqFiles,
	//	       std::vector<alignment_data>& allAlignments, string outfilename);
		bool aligned_to_end(alignment_data &algnmnt) {
			int minLeft = min(algnmnt.r1_al_start, algnmnt.r2_al_start);
			int minRight = min(algnmnt.length_r1-algnmnt.r1_al_end,
								algnmnt.length_r2-algnmnt.r2_al_end);
			return (minLeft <= 0 || minRight <= 0);
		}
		void computeAlignmentsFromOverlaps(char* inputOverlaps, bool loadperthread, std::string prevAlignments, std::string outfilename,
				char* inputSequences, bool cachedFastq, bool cachedOverlaps, const char* base_dir, std::vector<alignment_data>&, bool);
		int outputAlignmentsToFile(std::vector<alignment_data>*, std::string);
		int64_t outputAlignmentsToFileWithMpi(std::vector<alignment_data>*, std::string, int64_t);
	};

	inline std::ostream & operator<<(std::ostream & str, const Aligner::scoring_info& scoring) {
		return (str << "SeqAn settings (xdrop=" << scoring.xdrop
				<< " match=" << scoring.match
				<< " mismatch=" << scoring.mismatch
				<< " gap=" << scoring.gap
				<< ") BELLA settings ("
				<< " adaptive=" << scoring.adaptive
				<< " error_rate=" << scoring.error_rate
				<< " chernoff_delta=" << scoring.chernoff_delta
				<< " relaxation_margin=" << scoring.relaxation_margin
				<< " adaptive_phi=" << scoring.adaptive_phi
				<< " static_threshold=" << scoring.static_threshold
				<< " align_to_end=" << scoring.aligntoend
				<< ")");
	}


} /* namespace align */

#endif /* ALIGNER_H_ */
