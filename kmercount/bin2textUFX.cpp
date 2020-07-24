#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <vector>
#include <sstream> 
#include <limits>
#include <array>
#include <map>
#include <tuple>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../common/common.h"
#include "../common/mpi_common.h"

StaticVars _sv = NULL;

#include "../common/MPIUtils.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "Deleter.h"
#include "Friends.h"

//#define CHECKONLY

using namespace std;
template <class T>
bool from_string(T & t, const string& s, std::ios_base& (*f)(std::ios_base&))
{
   istringstream iss(s);
   return !(iss >> f >> t).fail();
}

std::string itos(int n)
{
   const int max_size = std::numeric_limits<int>::digits10 + 1 /*sign*/ + 1 /*0-terminator*/;
   char buffer[max_size] = {0};
   sprintf(buffer, "%d", n);
   return std::string(buffer);
}
		
int myrank;

static const char *USAGE = "Usage: b2tufx -k KMER_LENGTH ufx.bin [entriesToDump]";

int main(int argc, char ** argv)
{	
	int nprocs; 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    double start_time = MPI_Wtime();

	char option;

	while ((option = getopt(argc, argv, "k:")) != -1) {
		switch(option) {
			case 'k':
				KMER_LENGTH = atoi(optarg);
				if (KMER_LENGTH >= MAX_KMER_SIZE) { SDIE("-k %d is too large for this binary, please choose one compiled with at least MAX_KMER_SIZE=%d", KMER_LENGTH, 32*(KMER_LENGTH+1+31/32)); }
				break;
			default:
				SDIE("Unknown option '%c'\n\t%s\n", option, USAGE);
		}
	}

	if (optind >= argc || KMER_LENGTH <= 0) {
		cerr << USAGE << "\n\tThis mpi code will create ufx.bin.full.txt" << endl;
		exit(1);
	}
	int mpiret = 0;
	
	MPI_Datatype datatype;
        MPI_Type_contiguous(sizeof(ufxpack), MPI_CHAR, &datatype );
        MPI_Type_commit(&datatype);
        int dsize;
        MPI_Type_size(datatype, &dsize);

	int64_t entries2dump;

	struct stat filestatus;
  	stat( argv[ optind ], &filestatus );
  	if (myrank == 0) {
		cout << filestatus.st_size << " bytes\n";
	}
	int64_t fileentries = static_cast<int64_t>(filestatus.st_size) / static_cast<int64_t>(dsize);
	if(argc - optind <= 1)
	{
		entries2dump = fileentries;

	}
	else
	{
		int64_t n;
		from_string(n,string(argv[optind+1]),std::dec);

		if(fileentries > n)
		{
			entries2dump = n;
		}
		else
		{
			entries2dump = fileentries;	// all we have
		}
	}
	Kmer::set_k(KMER_LENGTH);

	int64_t len = 0, bufSize = 16*1024*1024;
	int switcher=0;
	char *buf[2];
	buf[0] = (char*) malloc(bufSize);
	buf[1] = (char*) malloc(bufSize);
	if (buf[0] == NULL || buf[1] == NULL) {
		cerr << "Could not allocate buffers!" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Request req[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
	MPI_Status mpistatus;
	MPI_File out;
	{
		char outName[MAX_FILE_PATH];
		assert(strlen(argv[optind]) < MAX_FILE_PATH + 10);
		strcpy(outName, argv[optind]);
		strcpy(outName + strlen(outName), ".txt.full");
                MPI_File_delete(outName, MPI_INFO_NULL); 
                MPI_Barrier(MPI_COMM_WORLD);
		mpiret = MPI_File_open(MPI_COMM_WORLD, outName, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &out);
		if (mpiret != MPI_SUCCESS) { cerr << "Could not open " << outName << "!" << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
		if (myrank == 0) {
			cout << "Reading: " << argv[optind] << "\nCreating: " << outName << endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	int64_t perproc = entries2dump / static_cast<int64_t>(nprocs);
	int64_t perme;
	if(myrank == nprocs-1)	perme = entries2dump - perproc*(nprocs-1);
	else	perme = perproc;

  	cout << "Rank " << myrank <<  " dumping " << perme << " entries\n";

	int64_t mybegin = perproc * static_cast<int64_t>(myrank) ;
	FILE * f = fopen(argv[optind], "r");
	if(!f)
	{
		cerr << "Problem reading binary input file\n";
		return 1;
	}
		
	ufxpack upack;
	fseek (f, mybegin * static_cast<int64_t>(dsize), SEEK_SET );

#ifdef DEBUG
	FILE *pFile, *pFile2;
	{
		char fpath[MAX_FILE_PATH];
		sprintf(fpath, "%s.txt", argv[optind]);
		pFile = fopen_rank_path(fpath, "w", myrank);
		sprintf(fpath, "%s.odd.txt", argv[optind]);
		pFile2 = fopen_rank_path(fpath, "w", myrank);
	}
#endif

	string nonACTG("ACTG");
	string nonEXT("ACTGFX");
	int maxcount = 1;
	int mincount = 1;
	
	MPI_Barrier(MPI_COMM_WORLD);
	for(int64_t i=0; i< perme; ++i)
	{
		if (!fread(&upack, dsize, 1, f))
            continue;
		Kmer kmer(upack.arr);
		if (len + KMER_LENGTH + 40 >= bufSize) {
#ifdef DEBUG
			fprintf(stderr, "Thread %d: Startig write of %lld bytes\n", myrank, (lld) len);
#endif
			mpiret = MPI_File_iwrite_shared(out, buf[switcher%2], len, MPI_BYTE, &req[switcher%2]);
			if (mpiret != MPI_SUCCESS) { cerr << "Could not write " << len << " bytes" << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
			len = 0;
			switcher++;

#ifdef DEBUG
			fprintf(stderr, "Thread %d: waiting for previous write to finish\n", myrank);
#endif
			mpiret = MPI_Wait(&req[switcher%2], &mpistatus);
			if (mpiret != MPI_SUCCESS) { cerr << "Could not MPI_Wait" << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
#ifdef DEBUG
			int c;
			MPI_Get_count(&mpistatus, MPI_BYTE, &c);
			fprintf(stderr, "Thread %d: finished writing of %d bytes\n", myrank, c);
#endif

		}

		string seq = kmer.toString();
#ifdef DEBUG
/*
		if(upack.leftmin > upack.leftmax || upack.rightmin > upack.rightmax || upack.count < upack.leftmax || upack.count < upack.rightmax)
		{
			fprintf(pFile2, "%s\t\t%d: %c[%d,%d]%c[%d,%d]\n", seq, upack.count,
				upack.left, upack.leftmin, upack.leftmax,
				upack.right, upack.rightmin, upack.rightmax);
		}
*/
#endif
		if(maxcount < upack.count)	maxcount = upack.count;
		if(mincount > upack.count)	mincount = upack.count;

		if(seq.find_first_not_of(nonACTG) != string::npos)
		{
			cerr <<  "depth first search match non: " << nonACTG << " on " << seq << endl;
		}
		/*
		string leftext = string(1, upack.left);	//std::string(size_t , char ) constructor
		string rightext = string(1, upack.right);
		if(leftext.find_first_not_of(nonEXT) != string::npos)
		{
				cerr <<  "depth first search match non: " << nonEXT << " on " << seq << " with left ext " << leftext << endl;
		}
		if(rightext.find_first_not_of(nonEXT) != string::npos)
		{
				cerr <<  "depth first search match non: " << nonEXT << " on " << seq << " with right ext " << rightext << endl;
		}
		*/
		//int mylen = sprintf(buf[switcher%2] + len, "%s\t\t%c[%d,%d]%c[%d,%d]\n", seq.c_str(), upack.left, upack.leftmin, upack.leftmax, upack.right, upack.rightmin, upack.rightmax*/);
		int mylen = sprintf(buf[switcher%2] + len, "%s\n", seq.c_str());
		assert(mylen > 0);
		assert(mylen < KMER_LENGTH + 40);
		len += mylen;
#ifdef DEBUG
		//fprintf(pFile, "%s\t\t%c[%d,%d]%c[%d,%d]\n", seq.c_str(), upack.left, upack.leftmin, upack.leftmax, upack.right, upack.rightmin, upack.rightmax);
		fprintf(pFile, "%s\n", seq.c_str());
#endif
	}
	fclose(f);
	if (len > 0) {
#ifdef DEBUG
		fprintf(stderr, "Thread %d: last write of %lld bytes\n", myrank, (lld) len);
#endif
		mpiret = MPI_File_iwrite_shared(out, buf[switcher%2], len, MPI_BYTE, &req[switcher%2]);
		if (mpiret != MPI_SUCCESS) { cerr << "Could not write" << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
		len = 0;
		switcher++;

		mpiret = MPI_Wait(&req[switcher%2], &mpistatus);
		if (mpiret != MPI_SUCCESS) { cerr << "Could not MPI_Wait" << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
#ifdef DEBUG
		int c;
		MPI_Get_count(&mpistatus, MPI_BYTE, &c);
		fprintf(stderr, "Thread %d: finished writing of %d bytes\n", myrank, c);
#endif
	}
	switcher++;
	mpiret = MPI_Wait(&req[switcher%2], &mpistatus);
	if (mpiret != MPI_SUCCESS) { cerr << "Could not MPI_Wait" << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
#ifdef DEBUG
	{
		int c;
		MPI_Get_count(&mpistatus, MPI_BYTE, &c);
		fprintf(stderr, "Thread %d: finished writing of %d bytes\n", myrank, c);
	}
#endif

	MPI_Barrier(MPI_COMM_WORLD);
	mpiret = MPI_File_close(&out);
	if (mpiret != MPI_SUCCESS) { cerr << "Could not MPI_Close" << endl; MPI_Abort(MPI_COMM_WORLD, 1); }

#ifdef DEBUG
	fclose(pFile);
	fclose(pFile2);
#endif

	int glmaxcount, glmincount;
	MPI_Reduce(&maxcount, &glmaxcount, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&mincount, &glmincount, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	if(myrank == 0)
	{
		cout << "Max count is " << glmaxcount << endl;
		cout << "Min count is " << glmincount << endl;
	}
	
	free(buf[0]);
	free(buf[1]);

	double elapsed_time = MPI_Wtime() - start_time;
	double slowest_time;
	CHECK_MPI( MPI_Allreduce(&elapsed_time, &slowest_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) );
	if (!myrank) {
		char bname[MAX_FILE_PATH];
		printf("Overall time for %s is %.2f s\n", get_basename(bname, argv[0]), slowest_time);
	}

	MPI_Finalize();
	return 0;
}
