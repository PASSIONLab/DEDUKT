#include <assert.h>
#include "MPIFileBuffer.h"

/*
struct {
	MPI_File fh;
	MPI_Request request[2];
	MPI_Status status[2];
	Buffer buffer[2];
	int active, mode;
} _MPIFileBuffer;
typedef _MPIFileBuffer *MPIFileBuffer;
*/

MPIFileBuffer initMPIFileBuffer(MPI_Comm comm, const char *filename, int mode, size_t bufferSize) {
	MPIFileBuffer mfb = (MPIFileBuffer) malloc_chk(sizeof(_MPIFileBuffer));
	CHECK_MPI( MPI_File_open(comm, filename, mode, MPI_INFO_NULL, &(mfb->fh)) );
	mfb->request[0] = MPI_REQUEST_NULL;
	mfb->request[1] = MPI_REQUEST_NULL;
	mfb->buffer[0] = initBuffer(bufferSize);
	mfb->buffer[1] = initBuffer(bufferSize);
	mfb->active = 0;
	mfb->mode = mode;
	return mfb;
}
void freeMPIFileBuffer(MPIFileBuffer mfb) {
	flushMPIFileBuffer(mfb);
	flushMPIFileBuffer(mfb);
	assert( mfb->request[0] == MPI_REQUEST_NULL && mfb->request[1] == MPI_REQUEST_NULL ); // nothing pending
	CHECK_MPI( MPI_File_close( &(mfb->fh) ) );
	freeBuffer(mfb->buffer[0]); mfb->buffer[0] = NULL;
	freeBuffer(mfb->buffer[1]); mfb->buffer[1] = NULL;
	free_chk(mfb);
}

void waitMPIFileBuffer(MPIFileBuffer mfb, int which) {
	assert(which == 0 || which == 1);
	if (mfb->request[which] != MPI_REQUEST_NULL) {
		CHECK_MPI( MPI_Wait( &(mfb->request[which]), &(mfb->status[which]) ) );
		mfb->request[which] = MPI_REQUEST_NULL;
		int len;
		Buffer whichBuffer = mfb->buffer[which];
		CHECK_MPI( MPI_Get_count( &(mfb->status[which]), MPI_BYTE, &len) );
		if (getLengthBuffer(whichBuffer) == 0) {
			// is 0 length, so this is for reading, set length
			resetRawBuffer(whichBuffer, len);
		} else {
			// for writing, reset buffer to 0
			assert(len == getLengthBuffer(whichBuffer));	
			resetBuffer(whichBuffer);
		}
	}
}

void flushMPIFileBuffer(MPIFileBuffer mfb) {
	int mode = mfb->mode;
	assert ( (mode & MPI_MODE_WRONLY) == MPI_MODE_WRONLY || (mode & MPI_MODE_RDWR) == MPI_MODE_RDWR );

	int active = mfb->active % 2;
	Buffer activeBuffer = mfb->buffer[active];
	int len =  getLengthBuffer( activeBuffer );
	if ( len > 0 ) {
		assert( mfb->request[active] == MPI_REQUEST_NULL ); // nothing pending
		assert(isValidBuffer(activeBuffer));
		CHECK_MPI( MPI_File_iwrite(mfb->fh, getStartBuffer(activeBuffer), len, MPI_BYTE, &(mfb->request[active])) );
	}
	active = ++(mfb->active) % 2;
	waitMPIFileBuffer(mfb, active);
}

void setSizeMPIFileBuffer(MPIFileBuffer mfb, size_t newsize) {
	CHECK_MPI( MPI_File_set_size(mfb->fh, newsize) );
}

size_t getSizeMPIFileBuffer(MPIFileBuffer mfb) {
	MPI_Offset size;
	CHECK_MPI( MPI_File_get_size(mfb->fh, &size) );
	return size;
}

void seekMPIFileBuffer(MPIFileBuffer mfb, size_t offset, int whence) {
	for (int active = 0; active < 2 ; active++) {
		waitMPIFileBuffer(mfb, active);
	}
	CHECK_MPI( MPI_File_sync( mfb->fh ) );
	CHECK_MPI( MPI_File_seek( mfb->fh, offset, whence) );
	if (offset == 0 && whence == MPI_SEEK_SET) mfb->active = 0;
}

size_t writeMPIFileBuffer(MPIFileBuffer mfb, const void *data, size_t size) {
	int active = mfb->active % 2;
	Buffer activeBuffer = mfb->buffer[active];
	if (getLengthBuffer(activeBuffer) + size > getSizeBuffer(activeBuffer)) {
		flushMPIFileBuffer(mfb);
	}
	active = mfb->active % 2;
	assert(mfb->request[active] == MPI_REQUEST_NULL); // nothing pending
	return writeBuffer(mfb->buffer[active], data, size);
}

size_t _readMPIFileBuffer(MPIFileBuffer mfb, void *data, size_t size) {
	size_t readBytes = 0;
	int active = mfb->active % 2;
	Buffer activeBuffer = mfb->buffer[active];
	size_t len = getLengthBuffer(activeBuffer);
	size_t pos = getPosBuffer(activeBuffer);
	size_t bufsize = getSizeBuffer(activeBuffer);
	assert (len >= pos);
	if (len == pos) {
		if (len == 0 && mfb->request[active] == MPI_REQUEST_NULL && mfb->request[ (active+1) % 2] == MPI_REQUEST_NULL) {
			// first read, start an extra non-blocking read
			assert(len == 0 && pos == 0 && bufsize > 0);
			assert(isValidBuffer(activeBuffer));
			CHECK_MPI( MPI_File_iread(mfb->fh, getStartBuffer(activeBuffer), bufsize, MPI_BYTE, &(mfb->request[active]) ) );
			mfb->active++;
			active = mfb->active % 2;
			activeBuffer = mfb->buffer[active];
			waitMPIFileBuffer(mfb, active); // should be noop
		}
		// start next non-blocking read
		resetBuffer(activeBuffer);
		bufsize = getSizeBuffer(activeBuffer);
		assert( mfb->request[active] == MPI_REQUEST_NULL ); // nothing pending
		CHECK_MPI( MPI_File_iread(mfb->fh, getStartBuffer(activeBuffer), bufsize, MPI_BYTE, &(mfb->request[active]) ) );
		mfb->active++;
		active = mfb->active % 2;
		waitMPIFileBuffer(mfb, active);
		activeBuffer = mfb->buffer[active];
		len = getLengthBuffer(activeBuffer);
		if (len == 0) {
			// wait for the last read
			mfb->active++;
			active = mfb->active % 2;
			waitMPIFileBuffer(mfb, active);
			activeBuffer = mfb->buffer[active];
			len = getLengthBuffer(activeBuffer);
		}
		pos = getPosBuffer(activeBuffer);
		bufsize = getSizeBuffer(activeBuffer);
	}
	if (len == 0) return 0;
	int bufferLen = len - pos;
	int readLen = bufferLen > size ? size : bufferLen;
	char *d = (char*) data;
	if (readLen > 0) {
		readBytes += readBuffer(activeBuffer, d, readLen);
		d += readBytes;
	}
	if (bufferLen != readLen && bufferLen > 0) {
		// call recursively exactly once
		readBytes += _readMPIFileBuffer(mfb, d, size - readLen);
	}
	return readBytes;
}

size_t readMPIFileBuffer(MPIFileBuffer mfb, void *data, size_t size) {
	int readBytes = _readMPIFileBuffer(mfb, data, size);
	return readBytes;
}

