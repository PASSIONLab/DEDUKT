#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <zlib.h>

#include "Friends.h"
#include "Kmer.hpp"
#include "readufx.h"
#include "../common/common_rankpath.h"

using namespace std;

#ifndef GZIP_EXT
#ifndef NO_GZIP
#define GZIP_EXT ".gz"
#else
#define GZIP_EXT ""
#endif
#endif

static FILE *ufx_fopen_track(const char *path, const char *mode)
{
    FILE *f = fopen(path, mode);
    if (!f) { fprintf(stderr, "Could not open ufx file %s with mode %s! %s!\n", path, mode, strerror(errno)); fflush(stderr); }
    return f;
}

inline int ufx_fclose_track(FILE *f)
{
    int status = fclose(f);
    if (status != 0) { fprintf(stderr, "Could not close ufx file! %s!\n", strerror(errno)); fflush(stderr); }
    return status;
}

FILE * OpenDebugFile(char * prefix, FILE * pFile, int my_rank)
{
    stringstream ss;
    string rank;
    ss << my_rank;
    ss >> rank;
    string ofilename = prefix;
    ofilename += rank;
    pFile  = ufx_fopen_track(ofilename.c_str(), "w");
    return pFile;
}

// input parameters: filename, dsize, nprocs, my_rank
// output: myshare, dsize
int UFXInitOpen(char * filename, int * dsize, int64_t *myshare, int nprocs, int my_rank, 
                int64_t *nEntries, int from_cached_io, int kmer_length, UFX_FILE *UFX_f)
{
    assert(UFX_f != NULL);
    UFX_f->f = NULL;
    UFX_f->gz = NULL;
    UFX_f->isCachedIO = 0;
    if (from_cached_io)
	UFX_f->isCachedIO = 1;

    if (my_rank == 0) 
        printf("kmer_length %d\n", kmer_length);

    Kmer::set_k(kmer_length);
    *dsize = sizeof(ufxpack);
    
    char my_fname[MAX_FILE_PATH] = "";
    if (UFX_f->isCachedIO) {
        snprintf(my_fname, MAX_FILE_PATH, "/dev/shm/%s%d" GZIP_EXT, filename, my_rank);
        get_rank_path(my_fname, my_rank);
    } else {
        strcpy(my_fname, filename);
    }
        
    struct stat stufx;
    if (stat(my_fname, &stufx) != 0) { fprintf(stderr, "Could not stat %s! %s\n", my_fname, strerror(errno)); return -errno; }
	
    int64_t filesize = stufx.st_size;	
    int64_t numentries = filesize / static_cast<int64_t>(*dsize);

    (*nEntries) = numentries;

    if (UFX_f->isCachedIO) {
        // need to get the number of entries from file .. additionally it is compressed now
        char entries_fname[1024];
        snprintf(entries_fname, 1024, "%s.entries", my_fname);
        FILE *f = ufx_fopen_track(entries_fname, "r");
        if (fscanf(f, "%lld\t%lld\n", (long long int*) nEntries, (long long int*) &numentries) != 2) {
            cerr << "Could not read number from " << entries_fname << endl;
            return -errno;
        }
        fclose(f);
    }

    if (my_rank == 0) {
        if (!UFX_f->isCachedIO)
            cout << "Filesize is " << filesize << " bytes, and number of records is " << numentries << endl;
        else
            cout << "Total number of records is " << (*nEntries) << endl;
    }
    if (UFX_f->isCachedIO) {
#ifndef NO_GZIP
        UFX_f->gz = gzopen(my_fname, "r");
#else
        UFX_f->f = ufx_fopen_track(my_fname, "r");
#endif
    } else {
        UFX_f->f = ufx_fopen_track(my_fname, "r+"); // read-write to support ftruncate later
    }
    if (UFX_f->gz == NULL && UFX_f->f == NULL) return -errno;

    if (!UFX_f->isCachedIO) {
        int64_t perproc = numentries / nprocs;
		int64_t begin = perproc * my_rank;
		if (my_rank == nprocs-1)
			*myshare = numentries - (nprocs-1)* (perproc);
		else		
			*myshare = perproc;
        fseek(UFX_f->f, begin * static_cast<int64_t>(*dsize ), SEEK_SET );
    } else {
        *myshare = numentries;
    }
    return 0;
}

void UFXClose(UFX_FILE *UFX_f)
{
    if (UFX_f->f) 
        ufx_fclose_track(UFX_f->f);
    if (UFX_f->gz)
        gzclose(UFX_f->gz);
    UFX_f->f = NULL;
    UFX_f->gz = NULL;
}

// inputs: f, dsize, requestedkmers, dmin
// outputs: kmersarr, counts, lefts, rights inside upacked
// returns: number of k-mers read (can be less than requestedkmers if end of file) (-1 for error)
// TODO finish refactoring for diBELLA
int64_t UFXRead(int dsize, char *** kmersarr, int ** counts, READIDS ** readlists, POSITIONS ** positionlists, int64_t requestedkmers, int my_rank, int kmer_length, UFX_FILE *UFX_f, bool reuse)
{
    if(UFX_f->f == NULL && UFX_f->gz == NULL){
        cerr << "Thread " << my_rank << ": Problem reading binary input file\n";
        return -1;
    }

    if (requestedkmers == 0)
        return 0;

    ufxpack* upacked = new ufxpack[requestedkmers];
    if (upacked == NULL) {
        cerr << "Thread " << my_rank << ": Problem allocating memory for " << requestedkmers << " UFX entries.\n";
        return -1;
    }

    int64_t totread = 0;
    if (UFX_f->f) {
        totread = fread(upacked, dsize, requestedkmers, UFX_f->f);
    } else {
        int maxRead = 64*1024*1024 / dsize; // read in batches of at most 64MB
        while (totread < requestedkmers) {
            int64_t toRead = requestedkmers - totread;
            if (toRead > maxRead) toRead = maxRead;
            int numRead = gzread(UFX_f->gz, upacked + totread, dsize * toRead) / dsize;
            if (numRead != toRead) {
                fprintf(stderr, "ERROR Thread %d: Could not read all the requested ufx entries in this batch (%lld out of %lld, total %lld) at %lld\n", my_rank, (long long int) numRead, (long long int) toRead, (long long int) requestedkmers, (long long int) totread);
                return -1;
            }
            totread += numRead;
        }
    }

    if (totread != requestedkmers) {
        fprintf(stderr, "ERROR Thread %d: Could not read the full ufx file (%lld out of %lld requested entries)\n", my_rank, (long long int) totread, (long long int) requestedkmers);
        return -1;
    }

    if(!reuse)  // OK in the last iteration too because the invariant (totread <= requestedkmers) holds
    {
        // (*kmersarr) is of type char**
        if (my_rank == 0) {
            cout << "Thread 0: Allocating memory for " << totread << " reads: " << (totread * (sizeof(char*) + kmer_length+1 + sizeof(int) + 2) ) << " bytes" << endl;
        }
        (*kmersarr) = (char**) malloc(sizeof(char*) * totread);
        if (*kmersarr == NULL) { cerr << "Thread " << my_rank << ": Could not allocate memory for " << totread << " reads" << endl; return -1; }
        for (int64_t i = 0; i < totread; i++) {
            (*kmersarr)[i] = (char*) malloc((kmer_length+1) * sizeof(char)); // extra character for NULL termination
            if ((*kmersarr)[i] == NULL) { cerr << "Thread " << my_rank << ": Could not allocate memory for read " << i << endl; return -1; }
        }
    
        *readlists = (READIDS*) malloc(sizeof(READIDS) * totread);
        *positionlists = (POSITIONS*) malloc(sizeof(POSITIONS) * totread);
        *counts = (int*) malloc(sizeof(int) * totread);
        if (*counts == NULL || *positionlists == NULL || *readlists == NULL) { cerr << "Thread " << my_rank << ": Could not allocate memory for 3 * " << totread << " reads" << endl; return -1; }
    }
   
   
	for(int64_t i=0; i< totread; ++i)
	{
		Kmer kmer(upacked[i].arr);
        
        // from C++11 standard:
        // 21.4.7.1 says that the pointer returned by c_str() must point to a buffer of length size()+1.
        // 21.4.5 says that the last element of this buffer must have a value of charT() -- in other words, the null character
		std::strcpy ((*kmersarr)[i], kmer.toString().c_str()); // (*kmersarr)[i] is of type char*
		(*counts)[i] = upacked[i].count;
		(*readlists)[i] = upacked[i].reads;
		(*positionlists)[i] = upacked[i].positions;
		//fprintf(stderr, "K-mer is named (%c) %s (%c) with count %d\n", (*lefts)[i], (*kmersarr)[i] , (*counts)[i]);
		//printf("Rank %d: K-mer is named %s (%c) with count %d\n", MYTHREAD, (*kmersarr)[i], (*counts)[i]);
	}
    
    delete [] upacked;

	return totread;
}

inline void setPosInCounts(char base, array<int,4> & array, int count)
{
   array[0] = 0;
   array[1] = 0;
   array[2] = 0;
   array[3] = 0;
   
   switch(base)
   {
      case 'A' :
         array[0] = count;
         break;
      case 'C' :
         array[1] = count;
         break;
      case 'G' :
         array[2] = count;
         break;
      case 'T' :
         array[3] = count;
         break;
   }
}

inline void copyCounts(array<int,4> & array, int *counts)
{
   array[0] = counts[0];
   array[1] = counts[1];
   array[2] = counts[2];
   array[3] = counts[3];
}

//! If cached_io is true, then offsetInUFX is not used by this function
/* commented-out since it's not used in the extracted standalone ufx code
int writeNewUFX(char *filename, int *merDepths, char *kmersAndExts, int64_t nNewEntries, int my_rank, int64_t offsetInUFX, int *leftCounts, int *rightCounts, int kmer_length, UFX_FILE *UFX_f)
{
   Kmer::set_k(kmer_length);
   ufxpack * packed = new ufxpack[nNewEntries];
   if (packed == NULL) { fprintf(stderr, "Thread %d: Could not allocate %lld new packed entries %lld bytes\n", my_rank, (long long int) nNewEntries, (long long int) nNewEntries * sizeof(ufxpack)); fflush(stderr); return -1; }
   array<int,4> leftCount;
   array<int,4> rightCount;
   int count;
   char base;
   int64_t i, posEntries = 0, posInCounts = 0;
   
   for (i = 0; i < nNewEntries; i++) {
      count = merDepths[i];
      Kmer kmer(&kmersAndExts[posEntries]);
      ufxpack upack;
      upack.arr = kmer.getArray();
      //base = kmersAndExts[posEntries+kmer_length];
      //setPosInCounts(base, leftCount, count);
      //base = kmersAndExts[posEntries+kmer_length+1];
      //setPosInCounts(base, rightCount, count);
      
      copyCounts(leftCount, &leftCounts[posInCounts]);
      copyCounts(rightCount, &rightCounts[posInCounts]);
      
      posInCounts += 4;
      posEntries += kmer_length + 2;
      
      PackIntoUFX(leftCount, rightCount, count, upack);
      packed[i] = upack;
   }

   char my_fname[MAX_FILE_PATH] = "";
   if (UFX_f->isCachedIO) {
       snprintf(my_fname, MAX_FILE_PATH, "/dev/shm/%s%d" GZIP_EXT, filename, my_rank);
       get_rank_path(my_fname, my_rank);
       UFXClose(UFX_f);
#ifndef NO_GZIP
       gzFile ufxFD = NULL;
       ufxFD = gzopen(my_fname, "a1");
       if (ufxFD == Z_NULL) { fprintf(stderr, "Thread %d: Could not gzopen %s for appending! %s!\n", my_rank, my_fname, strerror(errno)); fflush(stderr); return -1; }
       size_t bytesToWrite = nNewEntries * sizeof(ufxpack);
       size_t bytesWritten = gzwrite(ufxFD, packed, bytesToWrite);
       if (bytesToWrite != bytesWritten) { fprintf(stderr, "Thread %d: Could not append %lld new kmers to ufx file %s: (%lld bytes of %lld)! %s\n", my_rank, (long long int) nNewEntries, my_fname, (long long int) bytesWritten, (long long int) bytesToWrite, strerror(errno)); fflush(stderr); return -1; }
       int status = gzclose(ufxFD);  	
       if (status != Z_OK) { fprintf(stderr, "Thread %d: Could not close ufx file: %s! %s\n", my_rank, my_fname, strerror(errno)); fflush(stderr); return -1; }
#else
       FILE * ufxFD = NULL;
       ufxFD = ufx_fopen_track(my_fname, "a");
       size_t entriesWritten = fwrite(packed, sizeof(ufxpack), nNewEntries, ufxFD);
       if (entriesWritten != nNewEntries) { fprintf(stderr, "Thread %d: Could not append %lld new kmers to ufx file %s: (%lld entries written)! %s\n", my_rank, (long long int) nNewEntries, my_fname, (long long int) entriesWritten, strerror(errno)); fflush(stderr); return -1; }
       ufx_fclose_track(ufxFD);
#endif
   } else {
       FILE *ufxFD = UFX_f->f;
       strcpy(my_fname, filename);
       if (ufxFD == NULL) {
           ufxFD = ufx_fopen_track(my_fname, "r+");
       }
       int status = fseek(ufxFD, offsetInUFX, SEEK_SET);
       if (status != 0) { fprintf(stderr, "Thread %d: Could not fseek to %lld in %s.  %s!\n", my_rank, (long long int) offsetInUFX, my_fname, strerror(errno)); fflush(stderr); return -1; }
       size_t entriesWritten = fwrite(packed, sizeof(ufxpack), nNewEntries, ufxFD);
       if (entriesWritten != nNewEntries) { fprintf(stderr, "Thread %d: Could not write %lld new kmers (%llu bytes) to ufx file %s: (%lld entries written after %lld position) %s!\n", my_rank, (long long int) nNewEntries, (unsigned long long int) nNewEntries*sizeof(ufxpack), my_fname, (long long int) entriesWritten, (long long int) offsetInUFX, strerror(errno)); fflush(stderr); return -1; }
       if (UFX_f->f == NULL) ufx_fclose_track(ufxFD);
   }

   //fprintf(stderr, "Thread %d: writing new UFX to %s, number of records %ld\n", my_rank, my_fname, nNewEntries);

   delete [] packed;
   return 0;
}
*/

void DeAllocateAll(char *** kmersarr, int ** counts, char ** lefts, char ** rights, int64_t initialread)
{
    int64_t i;
    for (i = 0; i < initialread; i++)
        free((*kmersarr)[i]);
    if (*kmersarr)
        free(*kmersarr);
    if (*counts)
        free(*counts);
    if (*lefts)
        free(*lefts);
    if (*rights)
        free(*rights);
}

