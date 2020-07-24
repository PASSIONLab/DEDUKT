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
#include <mpi.h>
#include <sys/stat.h>

#include "../common/hash_funcs.h"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "Deleter.h"
#include "Friends.h"
#include "SNPFunctions.h"

using namespace std;

#define FREQLOWER 40
#define FREQHIGHER 50

#ifdef MPI_IO
#undef MPI_IO
#endif

int nprocs;
int myrank;

size_t LoadFile(int ratio, int64_t begin, int64_t * bound, MPI_Datatype datatype, int dsize, char * filename, vector< pair<Kmer,int> > & merged)
{

#ifdef MPI_IO
	MPI_File file, P1file, P2file;
	int rc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);    
	if (rc != MPI_SUCCESS) {	MPI_Abort(MPI_COMM_WORLD, -95); }
#endif

	MPI_Status status;
	int readbyMPI;

	vector< vector<kmerpack> > packeds(ratio);
	size_t read = bound[ratio*(myrank+1)-1]-begin;
	kmerpack * allPacked = new kmerpack[read];

#ifdef MPI_IO
	MPI_File_set_view(file, begin * dsize, datatype, datatype, "external32", MPI_INFO_NULL); 
	rc = MPI_File_read_all(file, allPacked, read, datatype, &status);	
	MPI_Get_count(&status, datatype, &readbyMPI);
	if(rc != MPI_SUCCESS) {	MPI_Abort(MPI_COMM_WORLD, -101); }
	MPI_File_close(&file);
	MPI_Barrier(MPI_COMM_WORLD);
#else
	FILE * f = fopen(filename, "r");
	fseek (f, begin * static_cast<int64_t>(dsize), SEEK_SET );
	readbyMPI = fread(allPacked, dsize, read, f);
	fclose(f);
#endif
	if(myrank == 0) 
	{ 
		cout << "MPI read_all [at proc 0] got " << readbyMPI << " elements for " << filename << endl; 
		copy(allPacked[0].arr.begin(), allPacked[0].arr.end(), ostream_iterator<uint64_t>(cout," ")); cout << allPacked[0].count << endl; 
	}
	read = 0;
	for(int i=0; i< ratio; ++i)
	{
		size_t count = bound[ratio*myrank+i]-begin;
		packeds[i].resize(count);
		copy(allPacked+read, allPacked+read+count, packeds[i].begin());
		read += count;	// this is local pointer
		begin += count;

		if(myrank == 0)
		{
			for(int j=1; j< packeds[i].size(); ++j)
			{
				if(packeds[i][j] < packeds[i][j-1])
				{
					cout << "Sorting issue: "<< endl;
					cout << "Loc[" << j-1 << "]:";
					copy(packeds[i][j-1].arr.begin(), packeds[i][j-1].arr.end(), ostream_iterator<uint64_t>(cout," "));
					cout << ", Loc[" << j << "]:";
					copy(packeds[i][j].arr.begin(), packeds[i][j].arr.end(), ostream_iterator<uint64_t>(cout," "));
					cout << endl;
				}
			}
		}
	}
	delete [] allPacked;
	MPI_Barrier(MPI_COMM_WORLD);
	MergeAll( packeds, merged);	
	MPI_Barrier(MPI_COMM_WORLD);
	return (2*read);
}

class pair_first_compare
{
public:
   bool operator () (const pair<Kmer,int> & l, const pair<Kmer,int> & r)
   {
	return (l.first < r.first);
   }
};

int main(int argc, char ** argv)
{
    	MPI_Init(&argc, &argv);
    	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    
    	if(argc < 8)
    	{
        	cout << "Usage: ./SNPidentify offspringboun offspringkmers p1boun p1kmers p2boun p2kmers kmerLength" << endl;
		cout << "__boun files describe the processor boundaries in the __kmers files" << endl;
		cout << "p1__ and p2__ denote the parent files while offspring is self explanatory" << endl;
		cout << "number of lines on the __boun files should be an exact multiple of the number of processors used" << endl;
        	return 0;
    	}
	KMER_LENGTH = atoi(argv[7]);
	if (KMER_LENGTH < 1 || KMER_LENGTH >= MAX_KMER_SIZE) { cerr << "Invalid kmer length " << KMER_LENGTH << " for MAX_KMER_SIZE=" << MAX_KMER_SIZE << endl; MPI_Abort(MPI_COMM_WORLD, 1); }
	Kmer::set_k(KMER_LENGTH); 	// set static global k-value

	int64_t * OSbound; 
	int64_t * P1bound; 
	int64_t * P2bound; 
	int OSratio, P1ratio, P2ratio;
	if(myrank == 0)
	{
		int64_t point;
		int count = 0;
		ifstream inputOS(argv[1]);	
		while(inputOS >> point) count++;
		if(count % nprocs != 0)
			cout << "Number of lines in OS boundary file is not an exact multiple of #processors" << endl;
		OSratio = count / nprocs;
		cout << "OS ratio is " << OSratio << endl;
		inputOS.clear();	
		inputOS.seekg (0, inputOS.beg);

		count = 0;
		ifstream inputP1(argv[3]);	
		while(inputP1 >> point) count++;
		if(count % nprocs != 0)
			cout << "Number of lines in P1 boundary file is not an exact multiple of #processors" << endl;
		P1ratio = count / nprocs;
		cout << "P1 ratio is " << P1ratio << endl;
		inputP1.clear();	
		inputP1.seekg (0, inputP1.beg);

		count = 0;
		ifstream inputP2(argv[5]);	
		while(inputP2 >> point) count++;
		if(count % nprocs != 0)
			cout << "Number of lines in P2 boundary file is not an exact multiple of #processors" << endl;
		P2ratio = count / nprocs;
		cout << "P2 ratio is " << P2ratio << endl;
		inputP2.clear();	
		inputP2.seekg (0, inputP2.beg);

		OSbound = new int64_t[OSratio*nprocs];
		P1bound = new int64_t[P1ratio*nprocs];
		P2bound = new int64_t[P2ratio*nprocs];

		for(int i=0; i< OSratio*nprocs; i++)	inputOS >> OSbound[i];
		for(int i=0; i< P1ratio*nprocs; i++)	inputP1 >> P1bound[i];
		for(int i=0; i< P2ratio*nprocs; i++)	inputP2 >> P2bound[i];
	}
	MPI_Bcast(&OSratio, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&P1ratio, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&P2ratio, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(myrank != 0)
	{
		OSbound = new int64_t[OSratio*nprocs];
		P1bound = new int64_t[P1ratio*nprocs];
		P2bound = new int64_t[P2ratio*nprocs];
	}

	MPI_Bcast(OSbound, OSratio*nprocs, MPI_LONG_LONG, 0, MPI_COMM_WORLD); 
	MPI_Bcast(P1bound, P1ratio*nprocs, MPI_LONG_LONG, 0, MPI_COMM_WORLD); 
	MPI_Bcast(P2bound, P2ratio*nprocs, MPI_LONG_LONG, 0, MPI_COMM_WORLD); 


	MPI_Datatype datatype;
  	MPI_Type_contiguous(sizeof(kmerpack), MPI_CHAR, &datatype );
  	MPI_Type_commit(&datatype);
  	int dsize;
  	MPI_Type_size(datatype, &dsize);

	int64_t OSbegin, P1begin, P2begin;
	OSbegin = (myrank == 0)? 0: OSbound[OSratio * myrank-1];
	P1begin = (myrank == 0)? 0: P1bound[P1ratio * myrank-1];
	P2begin = (myrank == 0)? 0: P2bound[P2ratio * myrank-1];

	vector< pair<Kmer,int> > OSmerged;
	size_t OSread = LoadFile(OSratio, OSbegin,  OSbound, datatype, dsize, argv[2], OSmerged);
	if(myrank == 0)	cout << "Offspring k-mers merged\n";

	vector< pair<Kmer,int> > P1merged;
	size_t P1read = LoadFile(P1ratio, P1begin,  P1bound, datatype, dsize, argv[4], P1merged);
	if(myrank == 0)	cout << "Parent 1 k-mers merged\n";

	vector< pair<Kmer,int> > P2merged;
	size_t P2read = LoadFile(P2ratio, P2begin,  P2bound, datatype, dsize, argv[6], P2merged);
	if(myrank == 0)	cout << "Parent 2 k-mers merged\n";
	
	ofstream out;
        OpenDebugFile("SNPfile",  out);
	MPI_Barrier(MPI_COMM_WORLD);

	if(!is_sorted(OSmerged.begin(), OSmerged.end())) MPI_Abort(MPI_COMM_WORLD, -91); 
	if(!is_sorted(P1merged.begin(), P1merged.end())) MPI_Abort(MPI_COMM_WORLD, -92);
	if(!is_sorted(P2merged.begin(), P2merged.end())) MPI_Abort(MPI_COMM_WORLD, -93);

	pair_first_compare pfc;
	for(size_t i=1; i< OSread; ++i)	// at this point __begin has the total number of k-mers read
	{
		if((OSmerged[i].first).equalUpToLastBase(OSmerged[i-1].first))	
		{
			if((OSmerged[i].second >= FREQLOWER) && (OSmerged[i].second <= FREQHIGHER) && (OSmerged[i-1].second >= FREQLOWER) && (OSmerged[i-1].second <= FREQHIGHER))
			{
				// Don't accept tri/quad-allelic markers (those that are not biallelic)
				if(i > 1 && i < OSread-1)
				{
					if((OSmerged[i].first).equalUpToLastBase(OSmerged[i-2].first) || (OSmerged[i].first).equalUpToLastBase(OSmerged[i+1].first))
						continue;	// non biallelic
				}
			
				// benefit of nested ifs is to avoid unnecessary map::find() calls
				// lower_bound returns an iterator pointing to the first element in the range [first,last) which does not compare less than val

				pair<Kmer,int> alleleA = OSmerged[i-1];
				pair<Kmer,int> alleleB = OSmerged[i];

				auto p1found = lower_bound(P1merged.begin(), P1merged.end(), alleleA, pfc);
				if (p1found != P1merged.end() && !pfc(alleleA, *p1found))	// alleleA exists in P1
				{
					if(!binary_search(P1merged.begin(), P1merged.end(), alleleB, pfc)) // alleleB does not exist in P1
					{
						auto p2found = lower_bound(P2merged.begin(), P2merged.end(), alleleB, pfc);
						if(p2found != P2merged.end() && !pfc(alleleB, *p2found))	// alleleB exists in P2
						{
							if(!binary_search(P2merged.begin(), P2merged.end(), alleleA, pfc))	// alleleA  doesn't exist in P2
							{
								string strA = alleleA.first.toString();
								string strB = alleleB.first.toString();
								out << strA.substr(0,KMER_LENGTH-1) << "\t";
								out << strA.at(KMER_LENGTH-1) << "\t" <<  p1found->second << "\t" << alleleA.second << "\t";
								out << strB.at(KMER_LENGTH-1) << "\t" <<  p2found->second << "\t" << alleleB.second << endl;
							}
						}	
					}
				}
				else	// both cases can't happen at the same time
				{			
					p1found = lower_bound(P1merged.begin(), P1merged.end(), alleleB, pfc);
					if(p1found != P1merged.end() && !pfc(alleleB, *p1found))	// alleleB exists in P1
					{
						if(!binary_search(P1merged.begin(), P1merged.end(), alleleA, pfc))	// alleleA doesn't exist in P1
						{
							auto p2found = lower_bound(P2merged.begin(), P2merged.end(), alleleA, pfc);
							if(p2found != P2merged.end() && !pfc(alleleA, *p2found))	// alleleA exists in P2
							{
								if(!binary_search(P2merged.begin(), P2merged.end(), alleleB, pfc))	// alleleB  doesn't exist in P2
								{
									string strA = alleleA.first.toString();
									string strB = alleleB.first.toString();
									out << strB.substr(0,KMER_LENGTH-1) << "\t";
									out << strB.at(KMER_LENGTH-1) << "\t" <<  p1found->second << "\t" << alleleB.second << "\t";
									out << strA.at(KMER_LENGTH-1) << "\t" <<  p2found->second << "\t" << alleleA.second << endl;
								}
							}	
						}
					}
				}
			}
		}
	}
	out.close();
	MPI_Finalize();
	return 1;
}

	

