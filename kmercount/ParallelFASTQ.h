#ifndef _PARALLEL_FASTQ_H_
#define _PARALLEL_FASTQ_H_

#include <mpi.h>
#include <string>
#include <vector>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <string.h>
#include <sstream>

#include "../common/Buffer.h"
#include "../common/mpi_common.h"
#include "../common/MPIUtils.h"
#include "../fqreader/fq_reader.h"

using namespace std;

class ParallelFASTQ
{
public:
    
    ParallelFASTQ(void)
        :nrecords(0), elapsed_t(0)
    {
        fqr = create_fq_reader();
        if (!fqr)
            DIE("Couldn't create fastq reader\n");
    }

    ~ParallelFASTQ()
    {
        close_fq(fqr);
        destroy_fq_reader(fqr);
    }

    void open(const char *filename, bool cached_io, const char* base_dir, long knownSize = -1)
    {
    	assert(base_dir != NULL);
        if(fqr->f) close_fq(fqr);
        open_fq(fqr, filename, cached_io, base_dir);
    }

    int get_max_read_len() {
        return fqr->max_read_len;
    }
    
    size_t fill_block(vector<string> & names, vector<string> & seqs, vector<string> & quals, size_t maxMemoryUsed)
    {
        double t = MPI_Wtime();
        size_t memUsed = 0, avgMemPerRead = 0;
        Buffer id = initBuffer(MAX_READ_NAME_LEN);
        Buffer seq = initBuffer(DEFAULT_READ_LEN);
        Buffer qual = initBuffer(DEFAULT_READ_LEN);;
        size_t records_read = 0;
        names.clear();
        seqs.clear();
        quals.clear();
        while ( memUsed + avgMemPerRead < maxMemoryUsed ) {
            if (!get_next_fq_record(fqr, id, seq, qual)) {
                elapsed_t += (MPI_Wtime() - t);
                freeBuffer(id);
                freeBuffer(seq);
                freeBuffer(qual);
                return records_read;
            }
            names.push_back(getStartBuffer(id));
            seqs.push_back(getStartBuffer(seq));
            quals.push_back(getStartBuffer(qual));
            records_read++;
            nrecords++;
            memUsed += getLengthBuffer(id) + getLengthBuffer(seq) + getLengthBuffer(qual) + 3 * sizeof(char*);
            avgMemPerRead = memUsed / records_read;
        }
        elapsed_t += (MPI_Wtime() - t);
        freeBuffer(id);
        freeBuffer(seq);
        freeBuffer(qual);
        return records_read;
    }

    int64_t getTotalRecordsRead() { return nrecords; }
    double get_elapsed_time() { return elapsed_t; }

private:
    fq_reader_t fqr;
    int64_t nrecords;
    double elapsed_t;
};

#endif
