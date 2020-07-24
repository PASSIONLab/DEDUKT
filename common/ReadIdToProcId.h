/*
 * ReadIdToProcId.h
 *
 *  Created on: Nov 20, 2017
 *      Author: mellis
 */

#ifndef COMMON_READIDTOPROCID_H_
#define COMMON_READIDTOPROCID_H_

#include "mpi_common.h" // include before defines.h
#include "defines.h"

bool isValidReadId(const ReadId* readRanges, const ReadId rId);
int inThreadRange(const ReadId* readRanges, const int procId, const ReadId rId);
int getOwnerProcId(const ReadId* readRanges, const int firstProcInSearch, const int nProcsInSearch, const ReadId readId);




#endif /* COMMON_READIDTOPROCID_H_ */
