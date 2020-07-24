#ifndef _READUFX_H_
#define _READUFX_H_

#include <zlib.h>

#include "Friends.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    FILE *f;
    gzFile gz;
    int isCachedIO;
} UFX_FILE;

FILE * OpenDebugFile(char * prefix, FILE * pFile, int myrank);
int UFXInitOpen(char * filename, int * dsize, int64_t * myshare, int nprocs, int myrank, int64_t *nEntries, int from_cached_io, int kmer_length, UFX_FILE *f);
void UFXClose(UFX_FILE *f);
//int64_t UFXRead(int dsize, char *** kmersarr, int ** counts, char ** lefts, char ** rights, int64_t requestedkmers, int dmin, double errorRate, int reuse, int myrank, int kmer_length, UFX_FILE *f);
int64_t UFXRead(int dsize, char *** kmersarr, int ** counts, READIDS ** readlists, POSITIONS ** positionlists, int64_t requestedkmers, int my_rank, int kmer_length, UFX_FILE *UFX_f, bool reuse);
void DeAllocateAll(char *** kmersarr, int ** counts, char ** lefts, char ** rights, int64_t initialread);

/* commented out since it's not currently used in the extracted standalone ufx */
//int writeNewUFX(char *filename, int *merDepths, char *kmersAndExts, int64_t nNewEntries, int my_rank, int64_t offsetInUFX, int *leftCounts, int *rightCounts, int kmer_length, UFX_FILE *f);

#ifdef __cplusplus
}
#endif
#endif
