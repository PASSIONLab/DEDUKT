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

#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "Deleter.h"
#include "Friends.h"

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

static const char *USAGE = "b2tcount -k kmerlength ufx.bin";
int main(int argc, char ** argv)
{	
	char option;
	int nprocs; 

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
		SDIE("%s", USAGE);
	}

    	MPI_Init(&argc, &argv);
    	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	
	MPI_Datatype datatype;
        MPI_Type_contiguous(sizeof(kmerpack), MPI_CHAR, &datatype );
        MPI_Type_commit(&datatype);
        int dsize;
        MPI_Type_size(datatype, &dsize);

	int64_t entries2dump;

	struct stat filestatus;
  	stat( argv[ optind ], &filestatus );
  	cout << filestatus.st_size << " bytes\n";
	int64_t fileentries = static_cast<int64_t>(filestatus.st_size) / static_cast<int64_t>(dsize);
	int64_t threshold = 0;
	if(argc-optind <= 1)
	{
		entries2dump = fileentries;

	}
	else
	{
		int64_t n;
		from_string(n,string(argv[optind + 1]),std::dec);

		if(fileentries > n)
		{
			entries2dump = n;
		}
		else
		{
			entries2dump = fileentries;	// all we have
		}

		if(argc-optind > 2)
		{
			from_string(threshold,string(argv[optind + 2]),std::dec);
			cout << "Dumping only entries above " << threshold;
		}
	}
	Kmer::set_k(KMER_LENGTH);

	int64_t perproc = entries2dump / static_cast<int64_t>(nprocs);
	int64_t perme;
	if(myrank == nprocs-1)	perme = entries2dump - perproc*(nprocs-1);
	else	perme = perproc;

  	cout << "dumping " << perme << " entries\n";

	int64_t mybegin = perproc * static_cast<int64_t>(myrank) ;
	FILE * f = fopen(argv[optind], "r");
	if(!f)
	{
		cerr << "Problem reading binary input file\n";
		return 1;
	}
		
	kmerpack kpack;
	fseek (f, mybegin * static_cast<int64_t>(dsize), SEEK_SET );

	string outname = string(argv[optind])+ itos(myrank)+string(".txt");
	FILE * pFile = fopen (outname.c_str(),"w");
	for(int64_t i=0; i< perme; ++i)
	{
		if (!fread(&kpack, dsize, 1, f))
            continue;
		Kmer kmer(kpack.arr);
		if(kpack.count > threshold)
			fprintf(pFile, "%s\t\t%d\n", kmer.toString().c_str(), kpack.count);
	}
	fclose(f);
	fclose(pFile);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 0)
	{
		string command =  "cat";
        	for(int i=0; i<nprocs; ++i)
        	{
                	command += " ";
			string name = string(argv[optind])+ itos(i)+string(".txt");
                	command += name;
		}
        	command += " > ";
		string finalname = string(argv[optind])+ string(".txt.full");
        	command += finalname;
        	cout << command << endl;
        	if (!system(command.c_str()))
                cout << "Could not executed " << command.c_str() << endl;
	}

	MPI_Finalize();
	return 0;
}
