#ifndef MPI_FILE_BUFFER_
#define MPI_FILE_BUFFER_

#include <mpi.h>

#include "../common/mpi_common.h"
#include "../common/Buffer.h"
#include "../common/MPIUtils.h"

#if defined (__cplusplus)
extern "C" {
#endif

typedef struct {
	MPI_File fh;
	MPI_Request request[2];
	MPI_Status status[2];
	Buffer buffer[2];
	int active, mode;
} _MPIFileBuffer;
typedef _MPIFileBuffer *MPIFileBuffer;

MPIFileBuffer initMPIFileBuffer(MPI_Comm comm, const char *filename, int mode, size_t bufferSize);
void freeMPIFileBuffer(MPIFileBuffer mfb);
void waitMPIFileBuffer(MPIFileBuffer mfb, int which);
void flushMPIFileBuffer(MPIFileBuffer mfb);
void setSizeMPIFileBuffer(MPIFileBuffer mfb, size_t newsize);
size_t getSizeMPIFileBuffer(MPIFileBuffer mfb);

void seekMPIFileBuffer(MPIFileBuffer mfb, size_t offset, int mode);
size_t writeMPIFileBuffer(MPIFileBuffer mfb, const void *data, size_t size);
size_t readMPIFileBuffer(MPIFileBuffer mfb, void *data, size_t size); 

#if defined (__cplusplus)
} // extern "C"
#endif


#endif // MPI_FILE_BUFFER_
