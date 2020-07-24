/*
 * readIdToProcId.h
 *
 *  Created on: Nov 20, 2017
 *      Author: mellis
 */

#ifndef COMMON_READIDTOPROCID_CPP_
#define COMMON_READIDTOPROCID_CPP_

#include <assert.h>

#include "ReadIdToProcId.h"


bool isValidReadId(const ReadId* readRanges, const ReadId rId) {
	return (rId > nullReadId && rId <= readRanges[THREADS-1]);
}

/*
 * Returns:
 *   -1 if the read Id is invalid
 *    0 if the read Id is not in range of <code>procId</code>
 *    1 if the read Id is in range
 */
int inThreadRange(const ReadId* readRanges, const int procId, const ReadId rId) {
	if (!isValidReadId(readRanges, rId)) return -1;
	if (procId == 0 && rId <= readRanges[procId]) return 1;
	else if (procId > 0 && rId > readRanges[procId-1] && rId <= readRanges[procId]) return 1;
	return 0;
}

/**
 * Assumes read_ranges contains the maximum (inclusive) read ID per processor.
 */
int getOwnerProcId(const ReadId* readRanges, const int firstProcInSearch, const int nProcsInSearch, const ReadId readId) {
	assert (readId > nullReadId);
	if (nProcsInSearch == 1) return firstProcInSearch;
	if (nProcsInSearch == 2) {
		if (readId > readRanges[firstProcInSearch]) return (firstProcInSearch+1);
		else return firstProcInSearch;
	}
	// else
	int half = nProcsInSearch / 2 ; // floor() implied using integer division truncation
	int midProc = firstProcInSearch + half;
	if (readId <= readRanges[midProc] && readId > readRanges[midProc-1]) return midProc;
	// recurse
	int nnprocs = half + (nProcsInSearch % 2);
	// over first half
	if (readId < readRanges[midProc]) { return getOwnerProcId(readRanges, firstProcInSearch, nnprocs, readId); }
	// else over second half
	ASSERT (readId <= readRanges[ firstProcInSearch + (nProcsInSearch-1) ], "firstProcInSearch="+ std::to_string(firstProcInSearch) +"readId=" + std::to_string(readId) + ", RHS=" + std::to_string(readRanges[ firstProcInSearch + (nProcsInSearch-1) ]) );
	return getOwnerProcId(readRanges, midProc, nnprocs, readId);
}


#endif /* COMMON_READIDTOPROCID_CPP_ */
