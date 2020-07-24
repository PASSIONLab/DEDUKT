/*
 * ReadOverlapper.cpp
 *
 */
#include <iostream>
#include <algorithm> // min(), ...
#include <assert.h>
#include <climits> // MAX_INT
#include <stddef.h> // offsetof()
#include <stdlib.h> // abs()
#include <string.h> // strcat
#include <sys/stat.h>
#include "../kmercount/Friends.h"

#include "ReadOverlapper.h"
#include "../common/readid_records.h"
#include "../common/MPIUtils.h"
#include "../common/ReadIdToProcId.h"
#include "../common/Buffer.h"


/*
 * Orders pairs by the first read, the second read, the position in the first read, the position in the second, in that order of priority.
 */
bool operator<( const ReadOverlapPair& a, const ReadOverlapPair& b) {
	/*
	// don't insert read pairs that only differ in positions, and that by less than k, from others already in the set
	if ( (a.readId1 == b.readId1 && a.readId2 == b.readId2)
			|| (a.readId1 == b.readId2 && a.readId2 == b.readId1) )
	{
		PosInRead ap1 = (a.readId1 == b.readId2 && a.readId2 == b.readId1)? a.posInRead2 : a.posInRead1;
		PosInRead ap2 = (a.readId1 == b.readId2 && a.readId2 == b.readId1)? a.posInRead1 : a.posInRead2;
		// if positions in read 1 differ by less than k, AND positions in read 2 differ by less than k, return false
		if (abs(ap1 - b.posInRead1) < KMER_LENGTH && abs(ap2-b.posInRead2) < KMER_LENGTH) return false;
	}
	*/

	// otherwise, order two reads by their ID's, then by their positions
	return ( a.readId1 < b.readId1 // comparing different reads, order the lesser one first (arbitrarily)
			|| (a.readId1 == b.readId1 && a.readId2 < b.readId2) // comparing different reads in the second id, order by the lesser second id
			|| (a.readId1 == b.readId1 && a.readId2 == b.readId2 // the reads are the same
				&& a.posInRead1 < b.posInRead1)
			|| (a.readId1 == b.readId1 && a.readId2 == b.readId2
				&& a.posInRead1 == b.posInRead1 && a.posInRead2 < b.posInRead2)
			);
	// implicitly returns false for read pairs that are exactly the same
}

namespace overlap {

ReadOverlapper::ReadOverlapper(int maxReadIds, ReadId* readRanges, KmerCountsType* kmercounts)
{
	CHECK_MPI( MPI_Comm_size(MPI_COMM_WORLD,&nprocs) );
	CHECK_MPI( MPI_Comm_rank(MPI_COMM_WORLD,&myrank) );
	overlapResultMap = new ReadOverlapMap();
	perproc_buffers = new set<ReadOverlapPair>[nprocs];
	max_num_readids = maxReadIds;
	read_ranges = readRanges;
	kmercounts_ptr = kmercounts;
}

/*
 * Currently returns, per Kathy's slides:
 * if r1_id odd and r1_id < r2_id +1, the entry <r1_id, r2_id>
 * if r1_id even and r1_id > r2_id +1, the entry <r1_id,r2_id>
 * else <r2_id,r1_id>
 */
ReadOverlapPair orderSymmetrically(ReadId r, ReadId s, PosInRead rp, PosInRead sp) {
	if ( (r % 2 != 0 && r < s+1)
			|| (r % 2 == 0 && r > s+1 ) )
		{  return ReadOverlapPair(r,s,rp,sp); } //make_pair( make_pair(r, s), make_pair(rp, sp) ); }
	else { return ReadOverlapPair(s,r,sp,rp); } //make_pair( make_pair(s, r), make_pair(sp, rp) ); }
}


/*
 * Currently, orders lower ReadId first
 *
 * Expects: that the ReadIds are not equal, but if they are, orders by the second readId-position pair
 */
inline ReadOverlapPair ReadOverlapper::getOrderedPair(ReadId r, ReadId s, PosInRead rp, PosInRead sp)
{
	if (r < s) { return ReadOverlapPair(r,s,rp,sp); } \
	return ReadOverlapPair(s,r,sp,rp);
}

/*
 * A read Id is never compared to itself
 *
 * Assumes
 * - perproc_buffers have all been initialized.
 * - read lists with fewer readIds than reliable_max contain non-null reads first and then nullReadIds (i.e. no non-null read Id's follow nullReadIds)
 * - kmers contained only in 1 read can be ignored here
 * - read Ids are not duplicated in the read Id lists
 */
void ReadOverlapper::extractReadPairs() {
	for (auto it = kmercounts_ptr->begin(); it != kmercounts_ptr->end(); it++) {
		const READIDS readlist = get<0>(it->second); // get the value (read list, position list, # occurrances) then get the read list from the tuple
		const POSITIONS poslist = get<1>(it->second);
		for (int i = 0; i < max_num_readids; i++) {
			if (readlist[i] == nullReadId) break;
			ASSERT( (readlist[i] > nullReadId), "i loop: readlist[" << i << "]=" << readlist[i]);
			for (int j = i+1; j <  max_num_readids; j++) { // pairs < i+1 were already compared
				if (readlist[j] == nullReadId) break;
				ASSERT( (readlist[j] > nullReadId), "j loop: readlist[" << j << "]=" << readlist[j]);
				ASSERT(readlist[i] != readlist[j],"");
				ReadOverlapPair oPair = orderSymmetrically(readlist[i], readlist[j], poslist[i], poslist[j]);// getOrderedPair( readlist[i], readlist[j], poslist[i], poslist[j] );
				ASSERT( (oPair.readId1 > nullReadId), "oPair.r1 changed, values are: r1="<<oPair.readId1<<", r2="<<oPair.readId2<<", p1="<<oPair.posInRead1<<", p2="<<oPair.posInRead2);
				int ownerOne = getOwnerProcId( read_ranges, 0, nprocs, oPair.readId1 );
				perproc_buffers[ownerOne].insert(oPair);
				//int ownerTwo = getOwnerProcId( read_ranges, 0, nprocs, oPair.readId2 );
				//if (ownerOne != ownerTwo) { perproc_buffers[ownerTwo].insert(oPair); }
			}
		}
	}
}

void ReadOverlapper::printPerProcBuffers() {
	/*
	for (int i = 0; i < nprocs; i++) {
		cout << "Rank " << myrank << "'s preproc_buffers[" << i << "] contains ";
		for (auto itr = perproc_buffers[i].begin(); itr != perproc_buffers[i].end(); itr++) {
			cout << " (" << (itr->first).first << ", " << (itr->first).second
				<< ") , (" << (itr->second).first << ", " << (itr->second).second << ")";
		}
		cout << endl;
	}
	*/
	// TODO reimplement
}

void ReadOverlapper::exchangeReadPairs() {
	// performance statistics
	const bool LOAD=0, TIME=!LOAD;
	double start_time=0.0, exch_time=0.0;
	double tot_packing_time=0.0, tot_exchange_time=0.0, tot_store_time=0.0;
	double avg_glob_packing_time=0.0, avg_glob_store_time=0.0, avg_glob_exchange_time=0.0;
	double min_glob_exchange_time=0.0, max_glob_exchange_time=0.0, tot_glob_exchange_time=0.0;
	int num_exchanges=0;
	int64_t global_mins[2]={0, 0}, global_maxs[2]={0, 0};

	// construct an MPI data type for read Id pairs with positions
	/*
	MPI_Datatype types[4] = {MPI_UINT64_T, MPI_UINT64_T, MPI_UINT16_T, MPI_UINT16_T};
	int blocks[4] = { 1, 1, 1, 1 };
	MPI_Aint uint64_ext, uint16_ext;
	MPI_Type_extent(MPI_UINT64_T, &uint64_ext); // TODO extent function was depricated...
	MPI_Type_extent(MPI_UINT16_T, &uint16_ext);
	MPI_Aint displacements[4];
	displacements[0] = static_cast<MPI_Aint>(0);
	displacements[1] = uint64_ext;
	displacements[2] = uint64_ext + uint64_ext;
	displacements[3] = uint64_ext + uint64_ext + uint16_ext;
	MPI_Datatype mpi_readpair_type;
	MPI_Type_struct(4, blocks, displacements, types, &mpi_readpair_type);
	MPI_Type_commit(&mpi_readpair_type);
	int bytesperelement;
	MPI_Type_size(mpi_readpair_type, &bytesperelement); // when using MPI custom data types, never forget to use the MPI size function correspondingly...
	*/

	Buffer sndBufferR = initBuffer(MAX_ALLTOALL_MEM);
	Buffer rcvBufferR = initBuffer(MAX_ALLTOALL_MEM);
	Buffer sndBufferP = initBuffer(MAX_ALLTOALL_MEM);
	Buffer rcvBufferP = initBuffer(MAX_ALLTOALL_MEM);

	ReadId* rcvBffrReadIds;
	PosInRead* rcvBffrPoss;
	ReadId* sndBffrReadIds;
	PosInRead* sndBffrPoss;

	int *sendcnts = new int[nprocs]();
	int *sdispls = new int[nprocs]();
	int *rcvcnts = new int[nprocs]();
	int *rdispls = new int[nprocs]();
	if (!sendcnts || !sdispls || !rcvcnts || !rdispls) {
		SDIE("Insufficient memory to complete the read-pair exchange. "
				"Try rerunning with more nodes, fatter nodes, or fewer threads per node.");
	}

	size_t bytes_per_buffer = nprocs*sizeof(int);

	int64_t sndtotcnt=0;
	int64_t prevcnt=0;

	bool moreToSend = 0, moreToRcv = 0;
	do {
		memset(sendcnts, 0, bytes_per_buffer);
		memset(sdispls, 0, bytes_per_buffer);
		memset(rcvcnts, 0, bytes_per_buffer);
		memset(rdispls, 0, bytes_per_buffer);
		resetBuffer(sndBufferR);
		resetBuffer(rcvBufferR);
		resetBuffer(sndBufferP);
		resetBuffer(rcvBufferP);

		// prepare send count and displacement buffers
		start_time = MPI_Wtime();
		sndtotcnt = 0;
		prevcnt = 0;
		for (int i = 0; i < nprocs; i++) {
			size_t numReadIdsToSend = perproc_buffers[i].size()*2; // 2 ReadIds and 2 positions per item
			sendcnts[i] = min(numReadIdsToSend, (size_t) INT_MAX-1); // INT_MAX is odd
			sdispls[i] = prevcnt;
			DBG("sending %d elements to rank %d at position sdispls[%d]=%d \n", sendcnts[i], i, i, sdispls[i]);
			prevcnt += sendcnts[i];
			sndtotcnt += sendcnts[i];
			if (sdispls[i] < 0 || sndtotcnt < 0) {
				cerr << myrank << " detected overflow in Alltoall, " << __FUNCTION__ << " line " << __LINE__ << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		tot_packing_time += MPI_Wtime() - start_time;

		// exchange send counts
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoall(sendcnts, 1, MPI_INT, rcvcnts, 1, MPI_INT, MPI_COMM_WORLD) );
		tot_exchange_time += MPI_Wtime() - start_time;

		// prepare receive and receive displacement buffers
		start_time = MPI_Wtime();
		uint64_t rcvtotcnt = 0;
		for (int i = 0; i < nprocs; i++) {
			rdispls[i] = rcvtotcnt;
			DBG("expecting %d elements from rank %d at position rdispls[%d]=%d \n", rcvcnts[i], i, i, rdispls[i]);
			rcvtotcnt += rcvcnts[i];
			if (rdispls[i] < 0 || rcvtotcnt < 0) {
				cerr << myrank << " detected overflow in Alltoall, " << __FUNCTION__ << " line " << __LINE__ << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}

		// resize receive buffers
		growBuffer(rcvBufferR, rcvtotcnt * sizeof(ReadId));
		rcvBffrReadIds = (ReadId*) getStartBuffer(rcvBufferR);
		growBuffer(rcvBufferP, rcvtotcnt * sizeof(PosInRead));
		rcvBffrPoss = (PosInRead*) getStartBuffer(rcvBufferP);
		DBG("successfully resized receive buffers of %d ReadIds and PosInReads \n", rcvtotcnt);

		// resize the send buffers
		growBuffer(sndBufferR, sndtotcnt * sizeof(ReadId));
		sndBffrReadIds = (ReadId*) getStartBuffer(sndBufferR);
		growBuffer(sndBufferP, sndtotcnt * sizeof(PosInRead));
		sndBffrPoss = (PosInRead*) getStartBuffer(sndBufferP);
		DBG("successfully resized send buffers of %d ReadIds and PosInRead \n", sndtotcnt);

		// fill send buffers
		for (int i = 0; i < nprocs; i++) {
			int index = 0;
			auto itr = perproc_buffers[i].begin();
			for (itr ; itr != perproc_buffers[i].end() && index < sendcnts[i] ; itr++)
			{
				ReadOverlapPair rpair = *itr;
				sndBffrReadIds[sdispls[i] + index] = rpair.readId1;
				sndBffrPoss[sdispls[i] + index] = rpair.posInRead1;
				index++;
				sndBffrReadIds[sdispls[i] + index] = rpair.readId2;
				sndBffrPoss[sdispls[i] + index] = rpair.posInRead2;
				index++;
			}
			perproc_buffers[i].erase(perproc_buffers[i].begin(), itr); // erases to the current itr location
			moreToSend = (!perproc_buffers[i].empty()? true : moreToSend);
			DBG("erased perproc_buffers[%d].size() is %d \n", i, perproc_buffers[i].size());
		}
		DBG("Successfully filled send buffers \n");
		tot_packing_time += MPI_Wtime() - start_time;

		// exchange data
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoallv(sndBffrReadIds, sendcnts, sdispls, MPI_UINT64_T, rcvBffrReadIds, rcvcnts, rdispls, MPI_UINT64_T, MPI_COMM_WORLD) );
#ifdef TIGHT_READS
		CHECK_MPI( MPI_Alltoallv(sndBffrPoss, sendcnts, sdispls, MPI_UINT16_T, rcvBffrPoss, rcvcnts, rdispls, MPI_UINT16_T, MPI_COMM_WORLD) );
#else
		CHECK_MPI( MPI_Alltoallv(sndBffrPoss, sendcnts, sdispls, MPI_UINT32_T, rcvBffrPoss, rcvcnts, rdispls, MPI_UINT32_T, MPI_COMM_WORLD) );
#endif
		exch_time = MPI_Wtime() - start_time;
		tot_exchange_time += exch_time;
		DBG("Successfully exchanged ReadIds and PosInReads \n");

		// report performance (not timed)
		const int SND=0, RCV=1;
		int64_t local_counts[2]={0,0};
		local_counts[SND] = sndtotcnt;
		local_counts[RCV] = rcvtotcnt;
		CHECK_MPI( MPI_Reduce(&local_counts, &global_mins, 2, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD) );
		CHECK_MPI( MPI_Reduce(&local_counts, &global_maxs, 2, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD) );
		double global_min_time = 0.0;
		CHECK_MPI( MPI_Reduce(&exch_time, &global_min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
		double global_max_time = 0.0;
		CHECK_MPI( MPI_Reduce(&exch_time, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
		serial_printf("ReadOverlapper:%s exchange iteration %d, %d bytes per ReadId, %d bytes per Position: sent min %lld IDs/Positions, sent max %lld IDs/Positions, recv min %lld IDs/Positions, recv max %lld IDs/Positions, in min %.3f s, max %.3f s\n",
				__FUNCTION__, num_exchanges, sizeof(ReadId), sizeof(PosInRead), global_mins[SND], global_maxs[SND], global_mins[RCV], global_maxs[RCV], global_min_time, global_max_time);
		//
		// loop through received data and put in map
		//
		start_time = MPI_Wtime();
		for (int i = 0; i < rcvtotcnt; i++) {
			pair<ReadId,ReadId> key = make_pair(rcvBffrReadIds[i], rcvBffrReadIds[i+1]);
			PosPair value = make_pair( rcvBffrPoss[i], rcvBffrPoss[i+1] );
			(*overlapResultMap)[ key ].push_back( value );
			i++; // intentional
		}
		tot_store_time += MPI_Wtime() - start_time;
		CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD)); // keep loops in sync
		CHECK_MPI( MPI_Allreduce(&moreToSend, &moreToRcv, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD));
		DBG("has more to send? %s \n", (moreToSend? "yes" : "no"));
		DBG("has more to receive? %s \n", (moreToRcv? "yes" : "no"));
		num_exchanges++;
	} while (moreToRcv); // || moreToSend implied in moreToRcv computation

	delete[] sendcnts;
	delete[] sdispls;
	delete[] rcvcnts;
	delete[] rdispls;
	freeBuffer(sndBufferR);
	freeBuffer(rcvBufferR);
	freeBuffer(sndBufferP);
	freeBuffer(rcvBufferP);

	DBG("waiting at barrier\n");
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	// output performance results
	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &min_glob_exchange_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &max_glob_exchange_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &tot_glob_exchange_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_exchange_time = tot_glob_exchange_time/THREADS;

	CHECK_MPI( MPI_Reduce(&tot_packing_time, &avg_glob_packing_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_packing_time /= THREADS;
	CHECK_MPI( MPI_Reduce(&tot_store_time, &avg_glob_store_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_store_time /= THREADS;

	serial_printf("ReadOverlapper:%s Average packing time: %.2f s \n", __FUNCTION__, avg_glob_packing_time);
	serial_printf("ReadOverlapper:%s Average storage time: %.2f s \n", __FUNCTION__, avg_glob_store_time);
	serial_printf("ReadOverlapper:%s Min, max, avg and total (sum) time (s) for ReadId pair exchange: %.3f (min), %.3f (max), %.3f (avg), %3f (sum)\n",
			 __FUNCTION__, min_glob_exchange_time, max_glob_exchange_time, avg_glob_exchange_time, tot_glob_exchange_time);
}

/*
 * Duplicate removal not necessary in the current version:
 *   1 position is recorded per k-mer per read.
 *   No two processors will extract k-mers from the same read.
 */
/*
void ReadOverlapper::sortAndRemoveDuplicatePositions() {
	for (auto map_itr = overlapResultMap->begin(); map_itr != overlapResultMap->end(); map_itr++)
	{
		std::sort(map_itr->second.begin(), map_itr->second.end());
		map_itr->second.erase(
			std::unique( map_itr->second.begin(), map_itr->second.end(), equalPosPairs),
			map_itr->second.end()
		);
	}
}
*/

void ReadOverlapper::extractAndExchangeOverlappingReads() {
	// performance statistics
	double start_time=0.0,
			local_alloc_time=0.0,
			local_transform_time=0.0,
			glob_avg_time=0.0;
	//
	// allocate exchange buffers
	//
	start_time = MPI_Wtime();
	DBG("Starting extract and exchange of overlapping reads with %d processors \n", nprocs);
	for (int i = 0; i < nprocs; i++) {
		ReadOverlapSet* buffer = new ReadOverlapSet();
		assert(buffer != NULL);
		perproc_buffers[i] = *buffer;
	}
	DBG("allocated %d per-processor read-overlap-pair buffers, beginning read-pair extraction \n", nprocs);
	local_alloc_time = MPI_Wtime() - start_time;
	//
	// process kmer -> [ ReadId list ] map (kmercounts) into exchange buffers
	//
	start_time = MPI_Wtime();
	extractReadPairs();
	local_transform_time = MPI_Wtime() - start_time;
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	CHECK_MPI( MPI_Reduce(&local_alloc_time, &glob_avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_avg_time /= THREADS;
	serial_printf("%s: Average time allocating per-thread buffers before computing AxA^T: %.2f s \n", __FUNCTION__, glob_avg_time);
	CHECK_MPI( MPI_Reduce(&local_transform_time, &glob_avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_avg_time /= THREADS;
	serial_printf("%s: Average time for transforming local map (kmer->[ReadId list] to ReadId->[ReadId list]): %.2f s \n", __FUNCTION__, glob_avg_time);

	//
	// exchange overlapping read pairs and store received
	//
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	DBG("beginning read-pairs exchange\n");
	exchangeReadPairs(); // timing code embedded
	// need barrier if not implicit in above (currently is)
	DBG("finished read-pairs exchange, deleting last buffers\n");
}

/* TODO: finish - not sure if we need/want it though
 * also, requires restructuring the map
void filterAndExtendSeeds() {
	ReadId lastPartner = nullReadId;
	PosInRead lastPosition;
	PosInRead currPosition;
	for (auto map_iter = overlapResultMap->begin(); map_iter != overlapResultMap->end(); map_iter++) {
		std::set<ReadTuple>* curr_set = &(map_iter->second);
		for (auto set_itr = curr_set.begin(); set_itr != curr_set.end(); set_itr++) {
			if (lastPartner != get<0>(set_itr)) {
				lastPartner = get<0>(set_itr); // correctness depends on ordered insertion
				continue;
			}
			// else

		}
	}
}
*/

void ReadOverlapper::exchangeReadNames(set<ReadId>* reqbuffers, unordered_map<ReadId,string>& readNameMap)
{
	// performance statistics
	const bool LOAD=0, TIME=!LOAD;
	double start_time=0.0, exch_time=0.0, tot_packing_time=0.0, tot_exchange_time=0.0, tot_store_time=0.0;
	double tot_glob_packing_time=0.0, tot_glob_store_time=0.0;
	double avg_glob_packing_time=0.0, avg_glob_store_time=0.0;
	double tot_glob_exchange_time=0.0, max_glob_exchange_time=0.0, min_glob_exchange_time=0.0, avg_glob_exchange_time=0.0;
	int num_exchanges=0;
	double avg_sending_sum=0.0;

	/*
	 * The maximum number of characters that can be exchanged per thread per round,
	 * due to MPI's limitation on counts and displacements (int representation).
	 * Assumes there are no more than 2^10 characters per name (uncompressed).
	 */
	size_t MAX_READS_PER_EXCHANGE_PER_THREAD = INT32_MAX/1024/THREADS;
	ASSERT(MAX_READS_PER_EXCHANGE_PER_THREAD > 0, "MAX_READS_PER_EXCHANGE is <= 0. Need to refactor read-name exchange...");

	int tot_readids_wanted=0, // restricts displacement indices to max handleable by MPI Alltoallv
		tot_cnt_owned=0,
		tot_chars_owned=0,
		tot_chars_rcvd=0;
	size_t cnts_temp=0;

	char* seqs_names_owned;
	char* rcvnames;
	char* chars_index;

	//TODO: these should be replaced with reusable buffers
	ReadId* wantedReadIds;
	ReadId* ownedReadIds;
	PosInRead* rcvseqlengths;
	PosInRead* sndseqlengths;

	int* cnts_for_wanted;
	int* displs_for_wanted_ids;
	int* cnts_for_owned;
	int* cnts_tobe_owned;
	int* displs_for_owned_ids;
	int* displs_for_owned_names;
	int* cnts_seqs_owned;
	int* dspls_seqs_owned;

	bool exch_more = 0;
	bool exch_more_global = 0;

	size_t namelen=0;
	std::string val;
	std::string* tostore;

	// initialize per thread buffers
	// naming convention:
	// 		"wanted" == this processor does not own and wants
	//		"owned"  == this processor owns and another processors wants
	// (I know, it's sad when spelling out a naming convention is necessary...)
	int bytes_per_buffer = sizeof(int)*THREADS;
	cnts_for_wanted = new int[THREADS]();			ASSERT(cnts_for_wanted != NULL,""); // TODO handle
	displs_for_wanted_ids = new int[THREADS]();		ASSERT(displs_for_wanted_ids != NULL,""); // TODO handle
	cnts_for_owned = new int[THREADS]();				ASSERT(cnts_for_owned != NULL,""); // TODO handle
	cnts_tobe_owned = new int[THREADS]();			ASSERT(cnts_for_owned != NULL,""); // TODO handle
	displs_for_owned_ids = new int[THREADS]();		ASSERT(displs_for_owned_ids != NULL,""); // TODO handle
	displs_for_owned_names = new int[THREADS]();		ASSERT(displs_for_owned_ids != NULL,""); // TODO handle
	cnts_seqs_owned = new int[THREADS]();			ASSERT(cnts_seqs_owned != NULL,""); // TODO handle
	dspls_seqs_owned = new int[THREADS]();			ASSERT(dspls_seqs_owned != NULL,""); // TODO handle

	do {
		//////////////////////////// Begin ReadId exchange ///////////////////////////////
		//
		// initialize buffer for requested number of ReadIds from each other processor
		//
		memset(cnts_for_wanted, 0, bytes_per_buffer);
		memset(cnts_for_owned, 0, bytes_per_buffer);
		memset(displs_for_wanted_ids, 0, bytes_per_buffer);
		memset(cnts_tobe_owned, 0, bytes_per_buffer);
		memset(displs_for_owned_ids, 0, bytes_per_buffer);
		memset(displs_for_owned_names, 0, bytes_per_buffer);
		memset(cnts_seqs_owned, 0, bytes_per_buffer);
		memset(dspls_seqs_owned, 0, bytes_per_buffer);

		start_time = MPI_Wtime();
		tot_readids_wanted = 0;
		for (int i = 0; i < THREADS; i++) {
			displs_for_wanted_ids[i] = tot_readids_wanted;
			cnts_temp = reqbuffers[i].size();
			exch_more = (cnts_temp > MAX_READS_PER_EXCHANGE_PER_THREAD? 1: exch_more);
			cnts_for_wanted[i] = min(cnts_temp, MAX_READS_PER_EXCHANGE_PER_THREAD );
			tot_readids_wanted += cnts_for_wanted[i];
			if (cnts_for_wanted[i] < 0 || tot_readids_wanted < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			DBG("Asking for %d reads from rank %d \n", cnts_for_wanted[i], i);
		}
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange number of ReadIds requested from/for each processor
		//
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoall(cnts_for_wanted, 1, MPI_INT, cnts_for_owned, 1, MPI_INT, MPI_COMM_WORLD) );
		exch_time = MPI_Wtime() - start_time;
		tot_exchange_time += exch_time;

		//
		// initialize buffers of ReadIds to send (requesting respective reads from other processors)
		//
		start_time = MPI_Wtime();
		wantedReadIds = (ReadId*) malloc(tot_readids_wanted * sizeof(ReadId));
		ASSERT(wantedReadIds != NULL,""); // TODO better handling
		DBG("successfully malloc'd send buffer of %d ReadIds \n", tot_readids_wanted);
		memset(wantedReadIds, 0, tot_readids_wanted * sizeof(ReadId));
		for (int i = 0; i < THREADS; i++) {
			int index = 0;
			auto itr = reqbuffers[i].begin();
			for (itr ; itr != reqbuffers[i].end() && index < cnts_for_wanted[i]; itr++)
			{
				wantedReadIds[displs_for_wanted_ids[i] + index] = (ReadId) *itr;
				index++;
			}
			ASSERT(index == cnts_for_wanted[i],"cnts_for_wanted["+to_string(i)+"]="+to_string(cnts_for_wanted[i])+", but index="+to_string(index));
			reqbuffers[i].erase(reqbuffers[i].begin(), itr); // erases to the current itr location
			DBG("erased reqbuffers[%d].size() is %zu \n", i, reqbuffers[i].size());
		}

		//
		// initialize receive and receive displacement buffers from previous alltoall
		//
		tot_cnt_owned = 0;
		for (int i = 0; i < THREADS; i++) {
			displs_for_owned_ids[i] = tot_cnt_owned;
			tot_cnt_owned += cnts_for_owned[i];
			if (cnts_for_owned[i] < 0 || tot_cnt_owned < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			DBG("Expecting %d requested ReadIds from rank %d at position rdispls[%d]=%d \n",cnts_for_owned[i], i, i, displs_for_owned_ids[i]);
		}
		ownedReadIds = (ReadId*) malloc(tot_cnt_owned * sizeof(ReadId));    		ASSERT(ownedReadIds != NULL,""); // TODO better handling
		DBG("Successfully malloc'd receive buffer of %d ReadIds \n",tot_cnt_owned);
		memset(ownedReadIds, 0, tot_cnt_owned*sizeof(ReadId));
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange ReadIds of requested names
		//
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoallv(wantedReadIds, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT64_T,
				ownedReadIds, cnts_for_owned, displs_for_owned_ids, MPI_UINT64_T, MPI_COMM_WORLD) );
		exch_time = MPI_Wtime() - start_time;
		tot_exchange_time += exch_time;
		DBG("successfully exchanged ReadIds \n");

		//
		// for each name requested by ReadID, calculate the sum of name lengths and displacements per thread
		//
		tot_chars_owned = 0;
		sndseqlengths = new PosInRead[tot_cnt_owned];
		for (int i = 0; i < THREADS; i++) {
			dspls_seqs_owned[i] += tot_chars_owned;
			for (int j = displs_for_owned_ids[i]; j < (i==THREADS-1? tot_cnt_owned : displs_for_owned_ids[i+1]); j++) { // j should index rcvbuffer directly
				ReadId currId = ownedReadIds[j];
				ASSERT(readNameMap.count(currId) > 0, "Rank "+to_string(MYTHREAD)+": readNameMap does not contain mapping for readId "+to_string(currId)+", yet was asked for it by thread "+to_string(i)+"\n");
				namelen = ( readNameMap[currId] ).length();
				DBG("name-to-send %lld has length %d \n", currId, namelen);
				sndseqlengths[j] = namelen;
				cnts_seqs_owned[i] += namelen;
			}
			tot_chars_owned += cnts_seqs_owned[i];
			if (tot_chars_owned < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		DBG("finished calculating name lengths, total %d chars \n", tot_chars_owned);
		rcvseqlengths = new PosInRead[tot_readids_wanted](); // one length for each ReadId I requested
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange name lengths and displacements
		//
		start_time = MPI_Wtime();
		// cnts and displacements stay the same since the last exchange because we're exchanging 1 length per ReadId
	#ifdef TIGHT_READS
		CHECK_MPI( MPI_Alltoallv(sndseqlengths, cnts_for_owned, displs_for_owned_ids, MPI_UINT16_T,
				rcvseqlengths, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT16_T, MPI_COMM_WORLD) );
	#else
		CHECK_MPI( MPI_Alltoallv(sndseqlengths, cnts_for_owned, displs_for_owned_ids, MPI_UINT32_T,
					rcvseqlengths, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT32_T, MPI_COMM_WORLD) );
	#endif
		tot_exchange_time += MPI_Wtime() - start_time;

		//
		// initialize receive buffers and package names for exchange
		//
		start_time = MPI_Wtime();
		seqs_names_owned = new char[tot_chars_owned];
		chars_index = seqs_names_owned;
		val = "";
		for (int i = 0; i < tot_cnt_owned; i++) {
			val = readNameMap[ ownedReadIds[i] ];
			val.copy( chars_index, sndseqlengths[i]);
			chars_index += sndseqlengths[i];
		}
		DBG("finished filling send buffers with names \n");

		// prep for reuse,
		std::memset(cnts_for_owned, 0, sizeof(int) * THREADS);

		//
		// init rcv buffers based on number of chars expected from each processor
		//
		tot_chars_rcvd = 0;
		for (int i = 0; i < THREADS; i++) {
			displs_for_owned_names[i] = tot_chars_rcvd;
			for (int j = displs_for_wanted_ids[i]; j < (displs_for_wanted_ids[i] + cnts_for_wanted[i]); j++) {
				cnts_for_owned[i] += rcvseqlengths[j];
				tot_chars_rcvd += rcvseqlengths[j];
			}
			if (tot_chars_rcvd < 0 || displs_for_owned_names[i] < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		rcvnames = (char*) malloc(tot_chars_rcvd * sizeof(char));
		DBG("successfully initialized receive buffers \n");
		memset(rcvnames, 0, tot_chars_rcvd * sizeof(char));
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange actual names
		//
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoallv(seqs_names_owned, cnts_seqs_owned, dspls_seqs_owned, MPI_CHAR,
				rcvnames, cnts_for_owned, displs_for_owned_names, MPI_CHAR, MPI_COMM_WORLD) );
		exch_time=MPI_Wtime() - start_time;
		tot_exchange_time += exch_time;
		DBG("exchanged sequence data, cleaning up \n");

		//
		// store names
		//
		start_time = MPI_Wtime();
		tostore = NULL;
		chars_index = rcvnames; // chars_index reuse //TODO what is the purpose of this?!
		DBG("storing received names \n");
		for (int i = 0; i < tot_readids_wanted; i++) { // for each sequence received
			if(rcvseqlengths[i] < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			DBG("tot_readids_wanted=%d, chars_index=%p, i=%d, rcvseqlengths[i]=%d for ReadId=%lld \n", tot_readids_wanted,
					chars_index, i, rcvseqlengths[i], wantedReadIds[i]);
			tostore = new std::string(chars_index, rcvseqlengths[i]);
			chars_index += rcvseqlengths[i];
			readNameMap[ wantedReadIds[i] ] = *tostore;
		}
		DBG("Finished storing local sequence names, doing final clean-up \n");
		tot_store_time += MPI_Wtime() - start_time;
		CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD)); // keep rounds in sync
		//
		// round clean-up
		//
		chars_index = NULL;
		delete[] rcvseqlengths; // TODO this should be a reusable buffer allocated once per exchange
		free(wantedReadIds); // TODO this should be a reusable buffer allocated once per exchange
		free(rcvnames); // TODO check
		free(ownedReadIds); // TODO this should be a reusable buffer allocated once per exchange
		//
		// anyone have anything left to exchange?
		//
		CHECK_MPI( MPI_Allreduce(&exch_more, &exch_more_global, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD));
		exch_more=0;
	} while (exch_more_global);

	//////////////////////////// END NAME TAG EXCHANGE ////////////////////////////

	//
	// do final cleanup
	//
	delete[] cnts_for_wanted;
	delete[] cnts_for_owned;
	delete[] cnts_tobe_owned;
	delete[] displs_for_owned_ids;
	delete[] displs_for_owned_names;
	delete[] displs_for_wanted_ids;
	delete[] cnts_seqs_owned;
	delete[] dspls_seqs_owned;
	delete[] sndseqlengths;
	delete[] seqs_names_owned;


	//
	// output performance statistics
	//
	CHECK_MPI( MPI_Reduce(&tot_store_time, &tot_glob_store_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_store_time = tot_glob_store_time/THREADS;

	CHECK_MPI( MPI_Reduce(&tot_packing_time, &tot_glob_packing_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_packing_time = tot_glob_packing_time/THREADS;

	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &min_glob_exchange_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &max_glob_exchange_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &tot_glob_exchange_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_exchange_time = tot_glob_exchange_time/THREADS;

	serial_printf("Finished read name exchange \n");
	serial_printf("ReadOverlapper:%s : total (sum) and average (packing) times were %.2f s and %.2f s \n",
		__FUNCTION__, tot_glob_packing_time, avg_glob_packing_time);
	serial_printf("ReadOverlapper:%s : total (sum) and average (store) times were %.2f s and %.2f s \n",
		__FUNCTION__, tot_glob_store_time, avg_glob_store_time);
	serial_printf("ReadOverlapper:%s : total (sum), average, min, and max (exchange) times (s) were: %.2f, %.2f, %.2f, %.2f \n",
		__FUNCTION__, tot_glob_exchange_time, avg_glob_exchange_time, min_glob_exchange_time, max_glob_exchange_time);
}

/*
 * Assumes:
 *   max_num_seeds >= 1
 * we don’t need to remove duplicate positions from the (ReadA,ReadB)->positions list
 * because k-mers are stored in one location and one position per read is stored for each k-mer,
 * so “duplicate” position pairs will only exist across read pairs not between read pairs
 *
 */
void ReadOverlapper::filterSeeds(int min_seed_distance, int max_num_seeds) {
	ASSERT(min_seed_distance >= 0 || max_num_seeds >= 1, "minimum seed distance or maximum number of seeds setting doesn't make sense");
	PosPair last_pos_pair;
	int num_seed_pairs = 0;
	for (auto map_itr=overlapResultMap->begin();
			map_itr != overlapResultMap->end(); map_itr++)
	{
		ASSERT( ((*map_itr).second.size() >= 1), "at least one pair of positions per overlapping read pair is expected");
		std::sort((*map_itr).second.begin(), (*map_itr).second.end(), ReadOverlapper::orderBefore);
		auto pos_itr = (*map_itr).second.begin();
		last_pos_pair = (*pos_itr); pos_itr++;
		num_seed_pairs = 1; // reset for each PosPair vector
		while (pos_itr != (*map_itr).second.end() && num_seed_pairs < max_num_seeds) // less than, using 0-based indexing
		{
			if (ReadOverlapper::nDistant(min_seed_distance, last_pos_pair, (*pos_itr) )) {
				last_pos_pair = *pos_itr;
				num_seed_pairs++;
				pos_itr++;
			}
			else { pos_itr = (*map_itr).second.erase(pos_itr); }
		}
		(*map_itr).second.erase(pos_itr, (*map_itr).second.end()); // erase the rest of them
	}
}

int ReadOverlapper::outputOverlaps(const char* filename, bool per_thread, bool cache_overlaps,
		unordered_map<ReadId,string>& readNameMap)
{
	if (per_thread || cache_overlaps) {
		char fname[MAX_FILE_PATH];
		memset(fname, '\0', sizeof(char)*MAX_FILE_PATH);
		if (cache_overlaps) { sprintf(fname, VMFILE_PATH); }
		strcat(fname, filename);
		get_rank_path(fname, myrank);
		DBG("(cache_overlaps=%d) output filename changed from %s to %s\n", cache_overlaps, filename, fname);
		// writes overlaps to per-thread file independently in parallel
		return outputOverlapsPerThread(fname, readNameMap);
	}
	// else
	return outputReadOverlapsMpi(filename, readNameMap);
}

/*
 * Writes <code>overlapResultMap<code> contents
 * (with read names stored or collected in <code>readNameMap</code>
 * to <code>outfilename</code>, independently in parallel.
 *
 * Assumes: <code>outfilename</code> is a per-thread file.
 */
int ReadOverlapper::outputOverlapsPerThread(
		const char* outfilename,
		unordered_map<ReadId,string>& readNameMap)
{
	DBG("writing read overlaps individually in parallel to %s\n", outfilename);
	double local_travs_time=0.0,
			local_io_time=0.0,
			glob_time=0.0;

	//
	// get all the read names we might not have
	//
	local_travs_time = MPI_Wtime();
	set<ReadId>* reqbuffers = new set<ReadId>[THREADS];
	for(auto map_itr = overlapResultMap->begin(); map_itr != overlapResultMap->end(); map_itr++) {
		ReadId readId1 = (ReadId) get<0>(map_itr->first);
		ReadId readId2 = (ReadId) get<1>(map_itr->first);
		if (readNameMap.count(readId1) <= 0) { // insert read 1
			 int owner = getOwnerProcId(read_ranges, 0, THREADS, readId1);
			 reqbuffers[owner].insert(readId1);
		 }
		 if (readNameMap.count(readId2) <= 0) { // insert read 2
			 int owner = getOwnerProcId(read_ranges, 0, THREADS, readId2);
			 reqbuffers[owner].insert(readId2);
		 }
	}
	local_travs_time = MPI_Wtime()-local_travs_time;

	exchangeReadNames(reqbuffers, readNameMap); // outputs its own time statistics

	local_io_time = MPI_Wtime();
	ofstream myoutputfile;
	myoutputfile.open (outfilename, std::ofstream::out | std::ofstream::trunc);
	if (!myoutputfile.is_open()) {
		DIE("Could not open %s: %s\n", outfilename, strerror(errno));
		return -1;
	}
	DBG("Successfully opened file %s\n", outfilename);
	for(auto map_itr = overlapResultMap->begin(); map_itr != overlapResultMap->end(); map_itr++) {
		ReadId readId1 = (ReadId) get<0>(map_itr->first);
		ReadId readId2 = (ReadId) get<1>(map_itr->first);
		ASSERT(readNameMap.count(readId1) > 0, "");
		ASSERT(readNameMap.count(readId2) > 0, "");
		string readName1 = readNameMap[ readId1 ];
		string readName2 = readNameMap[ readId2 ];
		myoutputfile << readId1 << ' ' << readId2 << ' ' << readName1 << ' ' << readName2;
		for(auto set_itr = map_itr->second.begin(); set_itr != map_itr->second.end(); set_itr++) {
			myoutputfile << ' ' << get<0>(*set_itr)		// position in read 1
						<< ' ' << get<1>(*set_itr);		// position in read 2
		}
		myoutputfile << '\n';
	}
	local_io_time = MPI_Wtime()-local_io_time;
	myoutputfile << flush; // not sure whether this should be included in the time
	myoutputfile.close();

	//debugging
	struct stat buffer;
	DBG("Rank %d: %s: stat %s = %d", myrank, __FUNCTION__, outfilename, (stat(outfilename, &buffer)) );

	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	delete[] reqbuffers;

	/* output timing statistics */
	CHECK_MPI( MPI_Reduce(&local_travs_time, &glob_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	serial_printf("ReadOverlapper:%s : max overlap-map traversal time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;

	CHECK_MPI( MPI_Reduce(&local_travs_time, &glob_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_time = glob_time/THREADS;
	serial_printf("ReadOverlapper:%s : avg overlap-map traversal time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;

	CHECK_MPI( MPI_Reduce(&local_io_time, &glob_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	serial_printf("ReadOverlapper:%s : max overlap output time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;

	CHECK_MPI( MPI_Reduce(&local_io_time, &glob_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_time = glob_time/THREADS;
	serial_printf("ReadOverlapper:%s : avg overlap output time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;
	/****/

	return 0;
}

/*
 * Writes <code>overlapResultMap<code> contents
 * (with read names stored or collected in <code>readNameMap</code>)
 * to shared file specified by <code>outfilename</code> via MPI I/O.
 *
 */
int ReadOverlapper::outputReadOverlapsMpi(const char* outfilename,
		unordered_map<ReadId,string>& readNameMap)
{
	DBG("writing read overlaps via MPI I/O to %s\n", outfilename);
	double local_travs_time=0.0,
		local_io_time=0.0,
		glob_time=0.0;

	MPI_File outfile;
	MPI_Status status;
	stringstream ss;
	string temp_string;
	size_t temp_length;
	Buffer out_buffer;
	size_t buffer_size;
	char* out_chars;
	int64_t my_num_bytes;
	int64_t tot_this_round;
	int64_t file_offset;
	int byte_count;
	int64_t my_offset;
	bool more_to_write;
	bool glob_more_to_write;

	// scan the overlap and read-name map for missing names
	local_travs_time = MPI_Wtime();
	set<ReadId>* reqbuffers = new set<ReadId>[THREADS];
	for(auto map_itr = overlapResultMap->begin(); map_itr != overlapResultMap->end(); map_itr++) {
		ReadId readId1 = (ReadId) get<0>(map_itr->first);
		ReadId readId2 = (ReadId) get<1>(map_itr->first);
		if (readNameMap.count(readId1) <= 0) { // insert read 1
			 int owner = getOwnerProcId(read_ranges, 0, THREADS, readId1);
			 reqbuffers[owner].insert(readId1);
		 }
		 if (readNameMap.count(readId2) <= 0) { // insert read 2
			 int owner = getOwnerProcId(read_ranges, 0, THREADS, readId2);
			 reqbuffers[owner].insert(readId2);
		 }
	}
	local_travs_time = MPI_Wtime()-local_travs_time;

	// request missing names from read owners
	exchangeReadNames(reqbuffers, readNameMap); // outputs its own time statistics

	// begin writing overlaps to file in parallel
	local_io_time = MPI_Wtime();
	if(myrank == 0) MPI_File_delete(outfilename, MPI_INFO_NULL);
	DBG("Rank 0 deleting file %s if exists\n", outfilename);
	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
	CHECK_MPI( MPI_File_open(MPI_COMM_WORLD, outfilename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile) );
	buffer_size = MAX_ALLTOALL_MEM; // arbitrary
	out_buffer = initBuffer(buffer_size);
	auto map_itr = overlapResultMap->begin();
	file_offset = 0;
	DBG("starting loop: copying data to buffer \n");
	do {
		more_to_write = 0;
		glob_more_to_write = 0;
		byte_count = 0;
		my_num_bytes = 0;
		tot_this_round = 0;
		my_offset = 0;
		resetBuffer(out_buffer);
		out_chars = getStartBuffer(out_buffer);
		// build the overlap pair output string
		for(; map_itr != overlapResultMap->end(); ) {
			ReadId readId1 = (ReadId) get<0>(map_itr->first);
			ReadId readId2 = (ReadId) get<1>(map_itr->first);
			ASSERT(readNameMap.count(readId1) > 0, "");
			ASSERT(readNameMap.count(readId2) > 0, "");
			string readName1 = readNameMap[ readId1 ];
			string readName2 = readNameMap[ readId2 ];
			ss << readId1 << ' ' << readId2 << ' ' << readName1 << ' ' << readName2;
			for(auto set_itr = map_itr->second.begin(); set_itr != map_itr->second.end(); set_itr++) {
				ss << ' ' << get<0>(*set_itr)  // position in read 1
				   << ' ' << get<1>(*set_itr);  // position in read 2
			}
			ss << '\n';
			temp_string = ss.str();
			temp_length = temp_string.length();
			if (buffer_size >= my_num_bytes + temp_length) {
				strcat(out_chars+my_num_bytes, temp_string.c_str());
				//memcpy(out_chars+my_num_bytes, temp_string.c_str(), temp_length);
				//strcatBuffer(out_buffer, temp_string.c_str());
				my_num_bytes += temp_length;
				map_itr++;
				ss.str(""); // reset the stringstream, clear() just resets the status flags
				ss.clear();
			}
			else {
				more_to_write = 1;
				break;
			}
		}
		DBG("finished loop: copying data to buffer\n");
		//assert(getLengthBuffer(out_buffer) == my_num_bytes);
		// MPI all write
		CHECK_MPI( MPI_Exscan(&my_num_bytes, &my_offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD) );
		if (my_offset < 0 || file_offset+my_offset < 0) {
			SDIE("[%s %d] Detected integer overflow while calculating MPI_Offset\n", __FUNCTION__, __LINE__);
		}
		DBG("my_num_bytes=%lld, my_offset=%lld, file_offset=%lld\n", my_num_bytes, my_offset, file_offset);
		ASSERT(sizeof(MPI_Offset) >= sizeof(int64_t), "detected potential conversion error: sizeof(MPI_Offset) NOT >= sizeof(int64_t)");
		DBG("assertion passed: sizeof(MPI_Offset) >= sizeof(int64_t). About to MPI_File_write_at_all with %lld bytes\n", my_num_bytes);
		CHECK_MPI( MPI_File_write_at_all(outfile, ((MPI_Offset) my_offset+file_offset), out_chars, my_num_bytes, MPI_CHAR, &status) );
		CHECK_MPI( MPI_Get_count(&status, MPI_CHAR, &byte_count) );
		assert( byte_count == my_num_bytes );
		CHECK_MPI( MPI_Allreduce(&more_to_write, &glob_more_to_write, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD));
		if (glob_more_to_write) {
			// calculate new offset
			CHECK_MPI(MPI_Allreduce(&my_num_bytes, &tot_this_round, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD));
			file_offset += tot_this_round;
			if (file_offset < 0) {
				SDIE("[%s %d] Detected integer overflow while updating file offset\n", __FUNCTION__, __LINE__);
			}
			DBG("%s looping again!", __FUNCTION__);
		}
	} while (glob_more_to_write);
	local_io_time = MPI_Wtime()-local_io_time;

	// begin cleanup
	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
	CHECK_MPI(MPI_File_close(&outfile));
	delete[] reqbuffers;
	freeBuffer(out_buffer);

	/* output timing statistics */
	CHECK_MPI( MPI_Reduce(&local_travs_time, &glob_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	serial_printf("ReadOverlapper:%s : max overlap-map traversal time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;

	CHECK_MPI( MPI_Reduce(&local_travs_time, &glob_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_time = glob_time/THREADS;
	serial_printf("ReadOverlapper:%s : avg overlap-map traversal time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;

	CHECK_MPI( MPI_Reduce(&local_io_time, &glob_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	serial_printf("ReadOverlapper:%s : max overlap output time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;

	CHECK_MPI( MPI_Reduce(&local_io_time, &glob_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_time = glob_time/THREADS;
	serial_printf("ReadOverlapper:%s : avg overlap output time: %.2f s \n", __FUNCTION__, glob_time);
	glob_time=0.0;
	/****/

	return 0;
}

// TODO fix: will hang if not all processors have some number of nonlocal items
/*
int ReadOverlapper::outputAndEraseOverlaps(string outfilename, unordered_map<ReadId,string>& readNameMap, ReadOverlapMap& overlapResultMap )
{
	size_t num_nonlocal_total=0,
		   num_nonlocal_temp=0;

	set<ReadId>* reqbuffers = new set<ReadId>[THREADS];
	for (int i = 0; i < THREADS; i++) {
		reqbuffers[i] = *(new set<ReadId>());
	}

	ofstream myoutputfile;
	myoutputfile.open (outfilename);
	if (!myoutputfile.is_open()) {
		DIE("Could not open %s: %s\n", outfilename, strerror(errno));
		return -1;
	}
	DBG("Successfully opened file %s\n", outfilename.c_str());

	do {
		if (num_nonlocal_temp > 0) { // shouldn't execute in the first loop by design
			CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
			exchangeReadNames(reqbuffers, readNameMap);
			myoutputfile.open (outfilename, ofstream::out | ofstream::ate | ofstream::app);
			if (!myoutputfile.is_open()) {
				DIE("Could not reopen %s: %s\n", outfilename, strerror(errno));
				return -1;
			}
		}
		num_nonlocal_total += num_nonlocal_temp;
		num_nonlocal_temp = 0;
		for(auto map_itr = overlapResultMap->begin(); map_itr != overlapResultMap->end();) {
			ReadId readId1 = (ReadId) get<0>(map_itr->first);
			ReadId readId2 = (ReadId) get<1>(map_itr->first);
			bool haveNameRead1 = (readNameMap.count(readId1) > 0);
			bool haveNameRead2 = (readNameMap.count(readId2) > 0);
			if (haveNameRead1 && haveNameRead2) {
				string readName1 = readNameMap[ readId1 ];
				string readName2 = readNameMap[ readId2 ];
				myoutputfile << readId1 << ' ' << readId2 << ' ' << readName1 << ' ' << readName2 << ' ';
				for(auto set_itr = map_itr->second.begin(); set_itr != map_itr->second.end(); set_itr++) {
					myoutputfile << get<0>(*set_itr) << ' '		// position in read 1
								 << get<1>(*set_itr) << ' ';	// position in read 2
				}
				myoutputfile << '\n';
				map_itr = overlapResultMap->erase(map_itr);
			}
			else {
				 if (!haveNameRead1) { // insert read 1
					 int owner = getOwnerProcId(read_ranges, 0, THREADS, readId1);
					 reqbuffers[owner].insert(readId1);
					 num_nonlocal_temp++;
				 }
				 if (!haveNameRead2) { // insert read 2
					 int owner = getOwnerProcId(read_ranges, 0, THREADS, readId2);
					 reqbuffers[owner].insert(readId2);
					 num_nonlocal_temp++;
				 }
				 map_itr++;
			}
		}
		myoutputfile.close();
	} while( !overlapResultMap->empty() );

	delete[] reqbuffers;

	return 0;
}
*/

/*
 * Writes <code>overlapResultMap<code> contents
 * (which excludes read names)
 * to <code>outfilename</code>, independently in parallel.
 *
 * Assumes: <code>outfilename</code> is a per-thread file.
 *
 * Doesn't output shared k-mer positions
 *
 * @deprecated
 */
int ReadOverlapper::outputOverlapsByIdOnly(string outfilename)
{
	ofstream myoutputfile;
	myoutputfile.open (outfilename);
	if (!myoutputfile.is_open()) {
		DIE("Could not open %s: %s\n", outfilename, strerror(errno));
		return -1;
	}
	DBG("Successfully opened file %s", outfilename);

	for(auto map_itr = (*overlapResultMap).begin(); map_itr != (*overlapResultMap).end(); map_itr++) {
		ReadId readId1 = (ReadId) get<0>(map_itr->first);
		ReadId readId2 = (ReadId) get<1>(map_itr->first);
		int num_shared_kmers = map_itr->second.size();
		myoutputfile << readId1 << ' ' << readId2 << ' ' << num_shared_kmers;
		myoutputfile << '\n';
	}
	myoutputfile.close();

	return 0;
}

ReadOverlapper::~ReadOverlapper() {
	if (perproc_buffers) delete[] perproc_buffers;
	if (overlapResultMap) delete overlapResultMap;
}

} /* namespace overlap */
