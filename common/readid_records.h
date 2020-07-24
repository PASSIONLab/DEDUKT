/*
 * readid_records.h
 *
 *  Created on: Jan 19, 2018
 *      Author: mellis
 */

#ifndef COMMON_READID_RECORDS_H_
#define COMMON_READID_RECORDS_H_

// #include "defines.h"

typedef std::array<ReadId, MAX_NUM_READS> READIDS;
typedef std::array<PosInRead, MAX_NUM_READS> POSITIONS; // currently records 1 position per (k-mer, read) pair



#endif /* COMMON_READID_RECORDS_H_ */
