/*
 * FileOpenerMPI.cpp
 *
 */
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <cstdint>

#include "../common/mpi_common.h"
#include "../common/defines.h"
#include "FileOpenerMPI.h"

vector< filenamesize >  loadAllFiles(char * filename)
{
	uint64_t totalsize = 0;
	int numfiles = 0;
	vector< filenamesize > filesview;
	if(MYTHREAD == 0)
	{
		filenamesize fdata;
		ifstream allfiles(filename);
		if(!allfiles.is_open()) {
				cerr << "Could not open " << filename << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
		}
		allfiles.getline(fdata.filename, MAX_FILE_PATH);
		while(!allfiles.eof())
		{
				struct stat st;
				stat(fdata.filename, &st);
				fdata.filesize = st.st_size;

				filesview.push_back(fdata);
				cout << filesview.back().filename << " : " << filesview.back().filesize / (1024*1024) << " MB" << endl;
				allfiles.getline(fdata.filename,MAX_FILE_PATH);
				totalsize += fdata.filesize;
				numfiles++;
		}
	}
	MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&totalsize, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	if(MYTHREAD != 0) { filesview.resize(numfiles); }

	MPI_Datatype MPI_filedata;
	MPI_Type_contiguous(sizeof(filenamesize), MPI_CHAR, &MPI_filedata);
	MPI_Type_commit(&MPI_filedata);
	MPI_Bcast(filesview.data(), numfiles, MPI_filedata, 0, MPI_COMM_WORLD);
    return filesview;
}
