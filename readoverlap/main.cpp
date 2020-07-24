/*
 * main.cpp
 *
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <utility>
#include <tuple>
#include <cerrno>
#include <iostream>
#include <fstream>

#include "../version.h"
#include "../common/MPIUtils.h"
#include "../common/common.h"
#include "../kmercount/kmerhashtable.h"
#include "../kmercount/readufx.h"

namespace overlap {
	int reliable_max = MAX_NUM_READS;
}

#include "ReadOverlapper.h"

#ifdef HISTOGRAM
#define COUNT_THRESHOLD 100000
#define COUNT_THRESHOLD_HIGH 10000000
#define HIGH_BIN 100
#define HIGH_NUM_BINS ((COUNT_THRESHOLD_HIGH-COUNT_THRESHOLD)/HIGH_BIN)
#endif

static int cores_per_node = 4;

KmerCountsType* myUfxMapPartition = NULL;

int printOverlapsToFile(string outfilename, ReadOverlapMap* overlapResultMap) {
	ofstream myoutputfile;
	myoutputfile.open (outfilename);
	if (!myoutputfile.is_open()) {
		DIE("Could not open %s: %s\n", outfilename, strerror(errno));
		return -1;
	}
	DBG("Successfully opened file %s", outfilename);

	for(auto map_itr = overlapResultMap->begin(); map_itr != overlapResultMap->end(); map_itr++) {
		myoutputfile << get<0>(map_itr->first) << ' ' 	// read1
					 << get<1>(map_itr->first) << ' '; 	// read2
		for(auto set_itr = map_itr->second.begin(); set_itr != map_itr->second.end(); set_itr++) {
			myoutputfile << get<0>(*set_itr) << ' '		// position in read 1
						 << get<1>(*set_itr) << ' ';	// position in read 2
		}
		myoutputfile << '\n';
	}
	myoutputfile.close();

	return 0;
}

/*
 * Assumes: k-mers are unique in the input file
 */
void buildMapFromUfx(char* input_UFX_name, int kmer_len, int cached_io=0) {
	//performance statistics
	double start_time=0.0, ufx_load_time=0.0, store_time=0.0;
	double avg_glob_load_time=0.0, avg_glob_store_time=0.0;
	//
	// read-in input file
	//
	int64_t size;
	int64_t myshare;
	int dsize;
	UFX_FILE UFX_f;
	start_time = MPI_Wtime();
	int status = UFXInitOpen(input_UFX_name, &dsize, &myshare, THREADS, MYTHREAD, &size, cached_io, kmer_len, &UFX_f);
	if (status != 0) { DIE("Could not open %s%s: %s\n", input_UFX_name, cached_io ? " (shm)" : "", strerror(-status)); }
	ufx_load_time += MPI_Wtime() - start_time;
	// TODO: replace minMemory calculation with sizes of types and number of entries used for hashtable about to be built
	/* {
		if (MYTHREAD == 0) {
		  int64_t minMemory = 12*(sizeof(list_t) + sizeof(int64_t)) * size / 10 / 1024/ 1024/ THREADS;
		  serial_printf("Minimum required shared memory: %lld MB. (%lld ufx kmers) If memory runs out re-run with more -shared-heap\n", (lld) minMemory, (lld) size);
		}
	} */

	//
	// read kmers etc.
	//
	//init_LookupTable(); // necessary for packing routines TODO check if necessary for unpacking routines
	char **kmersarr = NULL;
	int *counts = NULL;
	READIDS *readlists = NULL;
	POSITIONS *positionlists = NULL;
	start_time = MPI_Wtime();
	int64_t kmers_read = UFXRead(dsize, &kmersarr, &counts, &readlists, &positionlists, myshare, MYTHREAD, kmer_len, &UFX_f, 0);
	ufx_load_time += MPI_Wtime() - start_time;

	if (kmers_read < 0) DIE("There was a problem reading from ufx");
	DBG("Thread %d: Loaded %lld kmers\n", MYTHREAD, (lld) kmers_read);
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	if (MYTHREAD == 0)
	{	printf("Threads done with I/O\n"); }

	// TODO calculate required memory?

	// TODO consider adding kmer validity check, see hipmer/contigs/buildUFXhashBin.c lines ~96-103
    /* requires refactoring kmer_handling code
    char rc_kmer[MAX_KMER_SIZE];
	assert(kmer_len < MAX_KMER_SIZE);
	memset(rc_kmer, 0, MAX_KMER_SIZE);
	rc_kmer[kmer_len] = '\0';
	bool is_least;
	*/
	start_time = MPI_Wtime();
	ASSERT(myUfxMapPartition != NULL,"");
	for (int64_t i = 0; i < kmers_read; i++) {
		//TODO consider assigning and redistributing kmers to processors based on their hashval
		// for now, just insert them directly into the hash table
		/*
		reverseComplementKmer(kmersarr[i], rc_kmer, kmer_len);
		is_least = 0;
		is_least = (strcmp(kmersarr[i], rc_kmer) > 0) ? 0  : 1 ;
		if (!is_least) {
		 WARN("Found invalid kmer at %lld of %lld: %.*s vs %.*s\n",
			  (lld) ptr, (lld) kmers_read, kmer_len, kmersarr[i], kmer_len, rc_kmer);
		}
		assert(is_least);
		*/

#ifdef KHASH
		myUfxMapPartition->insert( Kmer(kmersarr[i]).getArray(), make_tuple(readlists[i], positionlists[i], counts[i]));
#else
		myUfxMapPartition->insert(make_pair(Kmer(kmersarr[i]).getArray(), make_tuple(readlists[i], positionlists[i], counts[i])));
#endif

	}
//	assert( ((int64_t) myUfxMapPartition->size()) == myshare); // if this fails, input kmers are not unique or something else went wrong...
	store_time = MPI_Wtime() - start_time;

	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

	CHECK_MPI( MPI_Reduce(&ufx_load_time, &avg_glob_load_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_load_time /= THREADS;

	CHECK_MPI( MPI_Reduce(&store_time, &avg_glob_store_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_store_time /= THREADS;

	serial_printf("Average kmermatch loading time: %.3f s \n", avg_glob_load_time);
	serial_printf("Average hash map insertion (building) time: %.3f s \n", avg_glob_store_time);
}

int overlap_main(int argc, char** argv) {
	if(MYTHREAD==0) cout << "starting overlap_main" << endl; // debugging
	// performance statistics
	double start_time=0.0, io_time=0.0;
	// parse arguments
	int KMER_LENGTH = MAX_KMER_SIZE-1;
	const char *base_dir = ".";
	char *input_ufx_name;
	string output_name = "";
	bool opt_err = false;
	option_t *optList, *thisOpt;
	optList = NULL;
	optList = GetOptList(argc, argv, "i:o:k:h:B");
	while (optList != NULL) {
		thisOpt = optList;
		optList = optList->next;
		switch (thisOpt->option) {
		 case 'i':
			 if(MYTHREAD==0) cout << "parsing i " << endl; // debugging
			 if (thisOpt->argument == NULL || thisOpt->argument[0] == '\0') opt_err = true;
			input_ufx_name = thisOpt->argument;
			if(MYTHREAD==0) cout << "finished parsing i " << endl; // debugging
			break;
		 case 'o':
			output_name = thisOpt->argument;
			break;
		 case 'k':
			 if(MYTHREAD==0) cout << "parsing k " << endl; // debugging
			 ASSERT(&(thisOpt->argument) != NULL,"");
			 if(MYTHREAD==0) cout << "k is " << thisOpt->argument << endl; // debugging
			 KMER_LENGTH = atoi(thisOpt->argument);
			if (KMER_LENGTH >= MAX_KMER_SIZE)
				SDIE("Please compile with a higher MAX_KMER_SIZE (%d) for kmer length %d", MAX_KMER_SIZE, KMER_LENGTH);
			if (GET_KMER_PACKED_LEN(KMER_LENGTH) > MAX_KMER_PACKED_LENGTH)
				SDIE("Please compile with a higher MAX_KMER_PACKED_LENGTH (%d) for kmer length %d", MAX_KMER_PACKED_LENGTH, KMER_LENGTH);
			break;
		case 'h':
			overlap::reliable_max = atoi(thisOpt->argument);
			if (overlap::reliable_max < 2) { SDIE("Upper limit set by -h must be greater than 1, -h= %d", overlap::reliable_max); }
			if (overlap::reliable_max > MAX_NUM_READS) { SDIE("Upper limit set by -h is larger than MAX_NUM_READS, -h= %d. Use compilation with MAX_NUM_READS > %d", overlap::reliable_max, MAX_NUM_READS); }
			break;
		 case 'B':
			 if(MYTHREAD==0) cout << "parsing B " << endl; // debugging
			base_dir = thisOpt->argument;
			 if(MYTHREAD==0) cout << "finished parsing B " << endl; // debugging
			break;
			/* TODO debug, and add back to GetOptList() call
		 case 'N':
			 if(MYTHREAD==0) cout << "parsing N " << endl; // debugging
			cores_per_node = atoi(thisOpt->argument);
			if(MYTHREAD==0) cout << "cores_per_node is " << cores_per_node << endl; // debugging
			SET_CORES_PER_NODE(cores_per_node);
			 if(MYTHREAD==0) cout << "finished parsing N " << endl; // debugging
			break;
			*/
		 default:
			opt_err = true;
			break;
		}
   }
	if (opt_err || !input_ufx_name)
	{
		if(MYTHREAD  == 0)
		{
			cout << "Usage: ./overlap -i input_ufx_file -k kmer_length -h reliable_kmer_max [-B base_directory -N cores_per_node -o output_file_name]" << endl;
			cout << "'input_ufx_file' is output of the kmer-analysis (kmermatch) stage " << endl;
			cout << "'kmer_length' is the kmer length (k) used for the original kmer-analysis (kmermatch) stage (should match the input_ufx_file)" << endl;
			cout << "'reliable_kmer_max' is the reliable-kmer-upper-limit used for the original kmer-analysis" << endl;
			cout << "'base_directory' is the directory in which to execute, and is optional" << endl;
			cout << "'cores_per_node' is the hardware defined number of cores per node and is optional" << endl;
		}
		return 0;
	}

	if(MYTHREAD  == 0)
	{
		cout << "diBELLA read-to-read overlap step" << endl;
		cout << "You are running with the following settings:" << endl;
		cout << "Input kmer-to-reads file = " << input_ufx_name << endl;
		cout << "K-mer length (minimum contiguous overlap for any read pairing) = " << KMER_LENGTH << endl;
		cout << "Max k-mer length (internal) = " << MAX_KMER_SIZE << endl;
		cout << "Reliable k-mer max = " << overlap::reliable_max << endl;
	}

	int cached_io = 0;
	/*
	if (strstr(base_dir, "/dev/shm")) {
	   cached_io = 1;
	   if(MYTHREAD==0) { cout << "cached_io flipping " << cached_io << endl; } // debugging
	   LOGF("Opening UFX with cached_io\n");
	   if(MYTHREAD==0) { cout << "finished logging " << cached_io << endl; } // debugging
	}
	*/
	if(MYTHREAD==0) { cout << "cached_io is " << cached_io << endl; } // debugging

	//
	// read-in ReadId partitions
	//
	// TODO is this code still relevant?
	ReadId* endReadIds = new ReadId[THREADS](); ASSERT(endReadIds != NULL,""); // todo
	if(MYTHREAD==0) { cout << "successfully initialized endReadId[] " << endl; } // debugging
	start_time = MPI_Wtime();
	if (MYTHREAD == 0) {
		ifstream file(INIT_READIDS_FNAME);
		if(!file || !file.is_open()) { DIE("Could not open %s%s: %s\n", INIT_READIDS_FNAME, cached_io ? " (shm)" : "", strerror(errno)); }
		string str;
		int i = 0;
		while(getline(file, str, ' ') || i < THREADS) {
			istringstream iss(str);
			iss >> endReadIds[i];
			i++;
		}
		if (i < THREADS || !file.eof()) { DIE("Unexpected in %s. Rerun overlap from kmermatch output %s with %d threads. \n", INIT_READIDS_FNAME, THREADS); }
		file.close();

#ifdef DEBUG
		cout << "Rank 0 read the following per-thread read range limits from " << INIT_READIDS_FNAME << ": ";
		for (int i = 0; i < THREADS; i++) { cout << endReadIds[i] << " "; }
		cout << endl;
#endif
	}
	CHECK_MPI(MPI_Bcast(endReadIds, THREADS, MPI_UINT64_T, 0, MPI_COMM_WORLD)); // included in io_time because the alternative is for all threads to read from file
	DBG("my end ReadId is %lld", endReadIds[MYTHREAD]);
	io_time += MPI_Wtime() - start_time;

	//
	// build kmer->(reads and positions) hash map
	//
	myUfxMapPartition = new KmerCountsType(); // initialize before building map
	buildMapFromUfx(input_ufx_name, KMER_LENGTH, cached_io); // timing embedded

	//
	// initialize ReadOverlapper and compute abstract read by read matrix
	//
	overlap::ReadOverlapper* myoverlapper = new overlap::ReadOverlapper(overlap::reliable_max, endReadIds, myUfxMapPartition);
	myoverlapper->extractAndExchangeOverlappingReads();

	//
	// output results
	//
	start_time = MPI_Wtime();
	if (output_name == "") { output_name = "overlaps-" + to_string(KMER_LENGTH); }
	output_name = output_name + "_" + to_string(MYTHREAD);
	//printOverlapsToFile(output_name);
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD)); // timing
	io_time += MPI_Wtime() - start_time;
	serial_printf("Additional I/O time: %.3f s \n", io_time);

/*
#ifdef HISTOGRAM
	//
	// count
	//
	vector<uint64_t> hist(COUNT_THRESHOLD,0);
	vector<uint64_t> hist_high(HIGH_NUM_BINS,0);
	for (auto map_itr = overlapResultMap->begin(); map_itr != overlapResultMap->end(); map_itr++)
	{
		uint64_t cnt = (map_itr->second).size(); // the set contains the position pairs at which the respective reads (keys in the map) align
		ASSERT(cnt > 0, "overlap position set size is <= 0 !!!");
		if (cnt <= COUNT_THRESHOLD) { ++(hist[cnt]); }
		else if (cnt <= COUNT_THRESHOLD_HIGH) { ++(hist_high[(cnt - COUNT_THRESHOLD)/HIGH_BIN]); }
		else { ++(hist_high[HIGH_NUM_BINS]); }
	}
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	//
	// reduce histogram data
	//
	if (MYTHREAD > 0) {
		CHECK_MPI( MPI_Reduce(hist.data(), NULL, COUNT_THRESHOLD, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );	// receive buffer not significant at non-root
		CHECK_MPI( MPI_Reduce(hist_high.data(), NULL, HIGH_NUM_BINS, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );
	}
	else {
		CHECK_MPI( MPI_Reduce(MPI_IN_PLACE, hist.data(), COUNT_THRESHOLD, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );
		CHECK_MPI( MPI_Reduce(MPI_IN_PLACE, hist_high.data(), HIGH_NUM_BINS, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );

		stringstream ss, ss1, ss2;
		//
		// init file names
		//
		ss << KMER_LENGTH;
		string hname = "histogram_k";
		hname += ss.str() + "_overlaps";
		string hnamehigh = hname + "_beyond";
		ss1 << COUNT_THRESHOLD;
		hnamehigh += ss1.str();
		hnamehigh += "_binsize";
		ss2 << HIGH_BIN;
		hnamehigh += ss2.str();
		hname += ".txt";
		hnamehigh += ".txt";
		hname = getRankPath(hname.c_str(), -1);
		hnamehigh = getRankPath(hnamehigh.c_str(), -1);
		//
		// output histograms to files
		//
		ofstream hout(hname.c_str());
		int bin = 0; // should be no less than 1 shared kmer per read pair
		for(auto it = hist.begin(); it != hist.end(); it++) {
			hout << ++bin << ' ' << *it << "\n";
		}
		ofstream hhigh(hnamehigh.c_str());
		bin = COUNT_THRESHOLD;
		for(auto it = hist_high.begin(); it != hist_high.end(); it++) {
				hout << bin << ' ' << *it << "\n";
				bin += HIGH_BIN;
		}
		cout << "Generated read-overlap histograms" << endl;
	}
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
#endif // HISTOGRAM
*/
	//
	// cleanup
	//
	if (myUfxMapPartition) delete myUfxMapPartition;
	if (endReadIds) delete endReadIds;
	if (myoverlapper) delete myoverlapper;
	return 0;
}

#ifndef SINGLE_EXEC
StaticVars _sv = NULL;
int main(int argc, char **argv)
{
    CHECK_MPI( MPI_Init(&argc, &argv) );
    OPEN_MY_LOG("overlap");
    serial_printf("Starting diBELLA version %s on %d threads\n", DIBELLA_VERSION, THREADS);
    int ret = overlap_main(argc, argv);
    MPI_Finalize();
    return ret;
}
#endif

