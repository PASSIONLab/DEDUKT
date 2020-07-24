#ifndef _MPI_COMMON_H
#define _MPI_COMMON_H

#include <stdio.h>
#include <unistd.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

inline int get_rank(void) {int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank;}
inline int get_num_ranks(void) {int num_ranks; MPI_Comm_size(MPI_COMM_WORLD, &num_ranks); return num_ranks;}
inline void common_exit(int x) { 
	if (x) fprintf(stderr, "Thread %d, calling MPI_Abort(%d)\n", get_rank(), x);
	MPI_Abort(MPI_COMM_WORLD, x);
	_exit(x);
}

#ifdef __cplusplus
}
#endif

#include "common.h"

#endif
