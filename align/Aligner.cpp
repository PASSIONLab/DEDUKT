/*
 * Aligner.cpp
 *
 *      Author: mellis
 */
#include "Aligner.h"
#include <math.h>
#include <string>

#define SEQAN // TODO

#ifdef SEQAN
#include "seqan_interface.h" // include last
#endif

#ifdef EDLIB
#include edlib_kernel.h
#endif // EDLIB

#ifndef MAX_ALLTOALL_MEM
#define MAX_ALLTOALL_MEM (128*1024*1024) /* 128 MB */
//#define MAX_ALLTOALL_MEM (12*1024*1024) /* 12 MB */
#endif

namespace align {

//const int OUT_BATCH_SIZE = 1000000; // arbitrary

void Aligner::checkLocality(const ReadOverlapPair& currPair, string name1, string name2, int64_t* overlap_stats) {
	int bothLocal = checkAndOrderOverlapPair(currPair);
	switch(bothLocal) {
	    case 3: {
	    	overlap_stats[0]++;
	    	break;
	    }
	   case 1: { // first is nonlocal, second local
		   bufferReadId( currPair.readId1);
		   checkNameMapping( currPair.readId2, name2 );
		   (*readNameMap)[currPair.readId1] = name1; // add to local set to check against remote later
		   overlap_stats[1]++;
		   break;
	   }
	   case 2: { // second nonlocal, first local
		   bufferReadId( currPair.readId2);
		   checkNameMapping( currPair.readId1, name1 );
		   (*readNameMap)[currPair.readId2] = name2; // add to local set to check against remote later
		   overlap_stats[1]++;
		   break;
	   }
	   case 0: { // both nonlocal
		   bufferReadId( currPair.readId1);
		   bufferReadId( currPair.readId2);
		   (*readNameMap)[currPair.readId1] = name1; // add to local set to check against remote later
		   (*readNameMap)[currPair.readId2] = name2; // add to local set to check against remote later
		   overlap_stats[2]++;
		   break;
	   }
	}
	overlap_stats[3]++;
}

/*
 * Stores <code>seqs</code> in [ReadID -> string (sequences)] map, and <code>names</code> [ReadID -> string (names)] map,
 * skipping those sequences of length <= KMER_LENGTH.
 */
void Aligner::storeSequences(const std::vector<std::string>& seqs, const std::vector<std::string>& names, ReadId startIndex) {
	ASSERT( seqs.size() == names.size() , "mismatching number of sequences and sequence names!");
	size_t nskipped = 0;
	for (int i = 0; i < seqs.size(); i++) {
		ASSERT( ((std::string) seqs[i]).length() > 0, "was given sequence with length 0 to store!");
		ASSERT( ((std::string) seqs[i]).length() <= USHRT_MAX, "read "+to_string(startIndex)+" has length "
				+to_string( ( (std::string) seqs[i]).length() )+", which is larger than max represented ("+to_string(USHRT_MAX)+")");
		if (seqs[i].length() > KMER_LENGTH) { // see kmermatch.cpp ParseNPack - skips sequences not meeting this condition
			(*readMap)[startIndex] = seqs[i];
			(*readNameMap)[startIndex] = names[i];
			startIndex++;
		}
		else { nskipped++; }
	}
	ASSERT( readMap->size() == readNameMap->size() , "mismatching number of sequences and sequence names - some sequences may have a duplicated name or something else went wrong!");
	LOGF("%s: skipped %lld sequences (too short) while storing sequences and names \n", __FUNCTION__, (lld) nskipped);
}

/*
 * Adapted from kmermatch.cpp::ProcessFiles()
 *
 * Since we're not calling ParseNPack, there may be sequences with unfiltered 'N's and
 *   sequences <= KMER_LENGTH
 *
 * We may need to optimize this to use less memory at a time
 */
ReadId Aligner::loadSequences(const std::vector<filenamesize> & allfiles, std::vector<std::string>& allreads,
		std::vector<std::string>& allnames, bool cachedIO, const char* base_dir)
{
	// Buffers were used in kmers exchange
	//Buffer scratch1 = initBuffer(initBufferSize); // program will exit with error if sufficient memory can't be allocated
    //Buffer scratch2 = initBuffer(initBufferSize);

	ReadId numReads = 0;
    double pfqTime = 0.0;
    auto files_itr = allfiles.begin();
    int trunc_ret;
    double t01 = MPI_Wtime(), t02;
    while(files_itr != allfiles.end())
    {
        ParallelFASTQ *pfq = new ParallelFASTQ();
        pfq->open(files_itr->filename, cachedIO, base_dir, files_itr->filesize);
        LOGF("Opened %s of %lld size\n", files_itr->filename, (lld) files_itr->filesize);
        files_itr++;
        // once again, arbitrarily chosen - see ProudlyParallelCardinalityEstimate
        std::size_t upperlimit = MAX_ALLTOALL_MEM / 16;

        std::vector<std::string> names;
        std::vector<std::string> seqs;
        std::vector<std::string> quals;

        int moreflags[2], allmore2go[2], anymore2go;
        int &moreSeqs = moreflags[0], &moreFiles = moreflags[1];
        int &allmoreSeqs = allmore2go[0], &allmoreFiles = allmore2go[1];
        moreSeqs = 1; // assume more as file is just open
        moreFiles = (files_itr != allfiles.end());
        std::size_t fill_status;
        do { // extract raw data into seqs and quals
            DBG("Starting new round: moreSeqs=%d, moreFiles=%d\n", moreSeqs, moreFiles);
            MPI_Pcontrol(1,"FastqIO");
            do {
                DBG2("Filling a block: moreSeqs=%d, moreFiles=%d\n", moreSeqs, moreFiles);
                // fill a block, from this or the next file
                if (pfq && !moreSeqs) {
                     ASSERT(pfq != NULL,"");
                     double t = pfq->get_elapsed_time();
                     pfqTime += t;
                     DBG2("Closed last file %.3f sec\n", t);
                     delete pfq;
                     pfq = NULL;
                     if (files_itr != allfiles.end()) {
                         pfq = new ParallelFASTQ();
                         pfq->open(files_itr->filename, cachedIO, base_dir, files_itr->filesize);
                         DBG2("Opened new file %s of %lld size\n", files_itr->filename, (lld) files_itr->filesize);
                         files_itr++;
                         moreSeqs = 1;
                     }
                }
                moreFiles = (files_itr != allfiles.end());
                fill_status = 0;
                if(moreSeqs) {
                    fill_status = pfq->fill_block(names, seqs, quals, upperlimit);  // file_status is 0 if fpos >= end_fpos
                    long long llf = fill_status;
                    DBG2("Filled block to %lld\n",(lld)  llf);
                    ASSERT(fill_status == seqs.size(),"");
                    ASSERT(fill_status == quals.size(),"");
                }
                moreSeqs = (fill_status > 0);
            } while (moreFiles && !moreSeqs);
            MPI_Pcontrol(-1,"FastqIO");

            allreads.insert(allreads.end(), seqs.begin(), seqs.end());
            allnames.insert(allnames.end(), names.begin(), names.end());
            numReads += seqs.size();

            anymore2go = moreSeqs+moreFiles;// allmoreSeqs + allmoreFiles;
        } while(anymore2go);

        if (pfq) {
            double t = pfq->get_elapsed_time();
            pfqTime += t;
            DBG( "Closing last file: %.3f sec\n", t);
            delete pfq;
            pfq = NULL;
        }

    }	// end_of_loop_over_all_files

    t02 = MPI_Wtime();
    double tots[2], gtots[2] = {0.0, 0.0};
    tots[0] = pfqTime;
    tots[1] = t02 - t01;
    LOGF("Process Total times: fastq: %0.3f elapsed: %0.3f\n", pfqTime, tots[1]);
    CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
    CHECK_MPI( MPI_Reduce(&tots, &gtots, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
    if (MYTHREAD == 0) {
        int num_ranks = THREADS;
        cout << "Aligner:"<<__FUNCTION__<<" Average time taken for FASTQ reads is " << (gtots[0] / num_ranks) << ", myelapsed " << tots[0] << endl;
        cout << "Aligner:"<<__FUNCTION__<<" Average time taken for elapsed     is " << (gtots[1] / num_ranks) << ", myelapsed " << tots[1] << endl;
        cout << "Aligner:"<<__FUNCTION__<<" Loaded reads of " << (files_itr == allfiles.end() ? " ALL files " : files_itr->filename) << " in " << t02-t01 << " seconds" << endl;
    }

    LOGF("Finished loading sequences, total of %lld reads \n", numReads);
    t02 = MPI_Wtime(); // redefine after prints

    //freeBuffer(scratch1);
    //freeBuffer(scratch2);

    ASSERT( allreads.size() == numReads,"" );
    ASSERT( allnames.size() == numReads,"" );
    return numReads;
}

/*
inline ReadId Aligner::getLocalIndex(ReadId globalIndex) {
	ReadId localId=(globalIndex - globalReadIdOffset); \
	ASSERT(localId > 0, "Rank "+to_string(MYTHREAD)+": global index translation failure: local="+to_string(localId)
		+", global="+to_string(globalIndex)+", my global offset="+to_string(globalReadIdOffset)+", max ReadId="+to_string(readRangeEnds[MYTHREAD])); \
	return localId;
}
inline ReadId Aligner::getGlobalIndex(ReadId localIndex) {
	ReadId globalId=(localIndex + globalReadIdOffset); \
	ASSERT(globalId > 0 && globalId <= readRangeEnds[MYTHREAD], "Rank "+to_string(MYTHREAD)+": global index translation failure: globalId="
			+to_string(globalId)+", local="+to_string(localIndex)+", my global offset="+to_string(globalReadIdOffset)+", max ReadId="+to_string(readRangeEnds[MYTHREAD])); \
	return globalId;
}
*/

void Aligner::initReadRangesAndSequences(const std::vector<filenamesize> & allfiles, bool cachedIO, const char* base_dir) {
	// load reads and read names
	std::vector<std::string>* allreads = new std::vector<std::string>();
	std::vector<std::string>* allnames = new std::vector<std::string>();
	ReadId numReads = loadSequences(allfiles, *allreads, *allnames, cachedIO, base_dir);

	// exchange number of reads-per-processor to calculate read indices
	readRangeEnds = new ReadId[THREADS]();
	uint64_t sndReadCounts[THREADS];
	for (int i = 0; i < THREADS; i++) { sndReadCounts[i] = numReads; }
	uint64_t recvReadCounts[THREADS];
	CHECK_MPI( MPI_Alltoall(sndReadCounts, 1, MPI_UINT64_T, recvReadCounts, 1, MPI_UINT64_T, MPI_COMM_WORLD) );
	ReadId readEndIndex = 0;
	for (int i = 0; i < THREADS; i++) { // would use scan/exscan but need all other processor's ranges too
		//if(MYTHREAD == THREADS-1) DBG("recvReadCounts[ %lld ]\n", recvReadCounts[i]);
		readEndIndex += recvReadCounts[i];
		readRangeEnds[i] = readEndIndex;
	}
	globalReadIdOffset = (MYTHREAD==0? 0 : readRangeEnds[MYTHREAD-1]);
	storeSequences(*allreads, *allnames, globalReadIdOffset+1); // read ID's start at 1...
	delete allreads; // all stored in map now
	delete allnames; // all stored in map now

	DBG("last Read Id is %lld \n", readRangeEnds[MYTHREAD] );
	DBG("my globalReadIdOffset= %lld \n", globalReadIdOffset );

	if (!MYTHREAD) {
		uint64_t total_reads = 0;
		for (int i = 0; i < THREADS; i++) {
			total_reads += recvReadCounts[i];
		}
		serial_printf("Alignment: Total number of reads is %lld and average per rank is %lld\n", total_reads, ( (total_reads + THREADS - 1) / THREADS) );
	}
}

/*
 * only call on reads extracted from fastq partition
 */
inline void Aligner::checkNameMapping(const ReadId read, string name) {
	ASSERT ( (*readNameMap)[read] == name, "Rank " << MYTHREAD << ": found mismatching names in overlap vs. fastq file, ID=" << read << ", fastq_name: "<< (*readNameMap)[read] << ", overlap_file_name: "<< name);
}

/*
 * Returns:
 * 			3 if both are local
 * 			2 if the second is nonlocal
 * 			1 if the first is nonlocal
 * 			0 if both are nonlocal
 */
int Aligner::checkAndOrderOverlapPair(const ReadOverlapPair& currPair) {
	int firstInRange = inThreadRange(Aligner::readRangeEnds, MYTHREAD, currPair.readId1);
	int secondInRange = inThreadRange(Aligner::readRangeEnds, MYTHREAD, currPair.readId2);
	ASSERT ( (firstInRange > -1 && secondInRange > -1), "Rank " << MYTHREAD << " encountered invalid ReadId(s) "
			<< currPair.readId1 << " and " << currPair.readId2 << " in input file. The min ReadId is 1, and the max ReadId is " << Aligner::readRangeEnds[THREADS-1]);
	if (firstInRange == 1 && secondInRange == 1) return 3;
	else if (firstInRange == 1) return 2; // second is nonlocal
	else if (secondInRange == 1) return 1; // first is nonlocal
	return 0; // neither is local
}

/*
bool Aligner::loadPrevAlignmentSet(std::string prevalignmentsname) {
	if ( (&prevalignmentsname != NULL) && (!prevalignmentsname.empty()) ) {
		std::cout << "Rank " << MYTHREAD << ": loadPrevAlignmentSet: attempting to open file by name: " << prevalignmentsname << std::endl;
		std::ifstream prevalignments(prevalignmentsname.c_str(), std::ifstream::in);
		if (prevalignments.good()) {
			std::cout << "Rank " << MYTHREAD << ": file is good" << std::endl;
			prevAlignmentsSet = new std::set< std::pair<ReadId,ReadId> >();
			ASSERT (prevAlignmentSet != NULL,"");
			//std::istream prevstream(prevalignments);
			std::string line;
			ReadId r1, r2;
			while (std::getline(prevalignments, line)) {
				std::stringstream linestream(line);
				linestream >> r1;
				linestream >> r2;
				prevAlignmentsSet->insert( std::make_pair(r1, r2) );
			}
			std::cout << "Rank " << MYTHREAD << ": found " << prevAlignmentsSet->size() << " previous alignments" << endl;
		}
		else {
			std::cout << "Rank " << MYTHREAD << ": file of previous alignments, " << prevalignmentsname << ", failed to open " << endl;
		}
		prevalignments.close();
	}
	return ((prevAlignmentsSet != NULL) && (prevAlignmentsSet->size() > 0) );
}
*/

inline void Aligner::bufferReadId(const ReadId readid) {
	int owner1 = getOwnerProcId(Aligner::readRangeEnds, 0, THREADS, readid);
	Aligner::perProcRequestIds[owner1].insert( readid );
}

void Aligner::loadOverlapsPerThread(const char* overlap_filename) {
	//
	// check whether the file can be accessed
	//
	struct stat buffer;
	if (stat(overlap_filename, &buffer) < 0) {
	    int err = errno;
	    SWARN("Could not stat %s file. Error is: %s \n", overlap_filename, strerror(err));
	    if (!MYTHREAD) cout << "["<< __FUNCTION__ << " " << __LINE__ << "] cout overlap_filename=" << overlap_filename << endl;
	}
	//
	// initialization
	//
	ifstream infile(overlap_filename);
	string line;
	ReadOverlapPair currPair;
	string name1, name2;
	vector< pair<PosInRead, PosInRead> > shared_positions;
	int64_t overlap_stats[4]; // both local, one local, neither local, total
	memset(&overlap_stats, 0, 4*sizeof(int64_t));

	//
	// make sure the file is open, abort if not
	//
	if (!infile.is_open()) {
		SDIE("Rank %d: [%s %s] could not open file %s\n",
				MYTHREAD, __FILE__, __LINE__, overlap_filename);
	}
	//
	// read/store the file contents
	//
	while ( getline(infile, line) ) {
		stringstream linestream(line);
		linestream >> currPair.readId1 >> currPair.readId2;
		linestream >> name1 >> name2;
		std::pair<ReadId,ReadId> key = std::make_pair(currPair.readId1, currPair.readId2);
		checkLocality(currPair, name1, name2, overlap_stats);
		while (linestream) {
			linestream >> currPair.posInRead1 >> currPair.posInRead2;
			std::pair<PosInRead,PosInRead> val = std::make_pair(currPair.posInRead1, currPair.posInRead2);
			shared_positions.push_back(val);
		}
		(*readPairs)[key].insert( (*readPairs)[key].end(), shared_positions.begin(), shared_positions.end());
		shared_positions.clear();
	}

	//
	// clean-up
	//
	infile.close();

	int64_t rcv_stats[4];
	memset(rcv_stats, 0, sizeof(int64_t)*4);
	int64_t total_overlaps=0;
	CHECK_MPI( MPI_Reduce(&overlap_stats, &rcv_stats, 4, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );
	total_overlaps = rcv_stats[3];
	for (int i = 0; i < THREADS; i++) {
		rcv_stats[i] = (rcv_stats[i]+THREADS-1)/THREADS;
	}
	serial_printf("%s:%s: AVG per rank overlap statistics: %lld both local, %lld one nonlocal, %lld both nonlocal, %lld total\n",
			__FILE__,__FUNCTION__,rcv_stats[0],rcv_stats[1],rcv_stats[2],rcv_stats[3]);

	memset(rcv_stats, 0, sizeof(int64_t)*4);
	CHECK_MPI( MPI_Reduce(&overlap_stats, &rcv_stats, 4, MPI_INT64_T, MPI_MIN, 0, MPI_COMM_WORLD) );
	serial_printf("%s:%s: MIN per rank overlap statistics: %lld both local, %lld one nonlocal, %lld both nonlocal, %lld total\n",
				__FILE__,__FUNCTION__,rcv_stats[0],rcv_stats[1],rcv_stats[2],rcv_stats[3]);

	memset(rcv_stats, 0, sizeof(int64_t)*4);
	CHECK_MPI( MPI_Reduce(&overlap_stats, &rcv_stats, 4, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD) );
	serial_printf("%s:%s: MAX per rank overlap statistics: %lld both local, %lld one nonlocal, %lld both nonlocal, %lld total\n",
				__FILE__,__FUNCTION__,rcv_stats[0],rcv_stats[1],rcv_stats[2],rcv_stats[3]);

	serial_printf("%s:%s: TOTAL number of overlaps globally is %lld\n", __FILE__, __FUNCTION__, total_overlaps);
}

/*
 *  not in-use and probably broken
 */
void Aligner::loadOverlaps(const char* overlapfile) {
	cout << "hello new function! input overlap file is: " << overlapfile << endl;
	int num_nonlocal_readids = 0;
	struct stat file_stat;
	size_t total_number_of_bytes;
	size_t partitioned_bytes;
	size_t max_num_bytes;
	size_t buffer_overlap;
	size_t my_num_bytes;
	size_t my_offset;
	string line;
	ReadOverlapPair currPair;
	string name1, name2;
	vector< pair<PosInRead, PosInRead> > shared_positions;

	// open the input file
	ifstream infile(overlapfile);
	if (!infile.is_open()) {
		DIE("Could not open the input file: %s...\n", overlapfile);
	}

	// calculate my partition and offset
	stat(overlapfile, &file_stat);
	total_number_of_bytes = file_stat.st_size;
	partitioned_bytes = total_number_of_bytes / THREADS;
	buffer_overlap = partitioned_bytes*0.2; // arbitrary multiplier
	max_num_bytes = partitioned_bytes + total_number_of_bytes % THREADS;
	my_num_bytes = (MYTHREAD == THREADS-1)? max_num_bytes
			: partitioned_bytes + buffer_overlap;
	my_offset = MYTHREAD * partitioned_bytes;

	// seek to my position in the file
	if (MYTHREAD > 0) {
		infile.seekg(my_offset-1);
		if (infile.get() != '\n') { getline(infile, line); }
	}

	//cout << "my_offset=" << my_offset << ", my_num_bytes=" << my_num_bytes << ", infile.tellg()=" << infile.tellg()<< endl;  //debugging
	//cout << "infile.tellg() < (my_offset + my_num_bytes) ? " << (infile.tellg() < (my_offset + my_num_bytes)? "true" : "false") << endl;  //debugging
	while ( infile.tellg() < (my_offset + my_num_bytes) && getline(infile, line) ) {
		//cout << line << endl; //debugging
		stringstream linestream(line);
		linestream >> currPair.readId1 >> currPair.readId2 >> name1 >> name2;
		std::pair<ReadId,ReadId> key = std::make_pair(currPair.readId1, currPair.readId2);
		// cout << "Parsed pair: " << currPair.readId1 << " " << currPair.readId2 << ", name1="<<name1<<", name2="<<name2<<endl; //debugging
		//checkLocality(currPair, name1, name2, num_nonlocal_readids); // commented because signature is outdated and this function isn't being used
		while (linestream) {
			linestream >> currPair.posInRead1 >> currPair.posInRead2;
			std::pair<PosInRead,PosInRead> val = std::make_pair(currPair.posInRead1, currPair.posInRead2);
			shared_positions.push_back(val);
		}
		(*readPairs)[key].insert((*readPairs)[key].end(), shared_positions.begin(), shared_positions.end());
		shared_positions.clear();
	}

	// clean-up
	infile.close();
}

/*
 * Assumes
 * - at least one of the two reads in a read pair is local
 * - lines in overlap input correspond 1:1 to lines in (previous or current) alignment output
 *   (one overlap list or alignment, respectively, per read pair, per line)
 *
 *   TODO: this code is not robust, the partitioning for various inputs needs work
 *   TODO: use collective reading
 *   TODO: use Buffer instances
 */
int Aligner::loadOverlapsMpi(char* infilename, bool checkPrev) {
	// TODO revisit types and assertions
	ASSERT(readRangeEnds != NULL,"");
	DBG( "Input file name is: %s \n", infilename);

	MPI_File filehandle;
	MPI_Status status;
	int rcv_count;
	int64_t total_number_of_bytes;
	int64_t max_number_of_bytes;
	int64_t my_offset;
	size_t number_of_bytes_partitioned;
	size_t my_number_of_bytes;
	MPI_Offset buffer_overlap;
	char *read_buffer;
	char *start_ptr;
	char *end_ptr;
	const char* delim = " \0";
	char *name1, *name2;
	char *temp;
	ReadOverlapPair currPair;
	vector< pair<PosInRead, PosInRead> > shared_positions;
	int64_t num_nonlocal_readids[4];
	memset(num_nonlocal_readids, 0, sizeof(int64_t)*4);

	CHECK_MPI(MPI_File_open(MPI_COMM_WORLD, infilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &filehandle));

	CHECK_MPI(MPI_File_get_size(filehandle, (MPI_Offset*) &total_number_of_bytes));
	DBG("total_number_of_bytes = %lld in file %s\n", total_number_of_bytes, infilename);
	assert(total_number_of_bytes > 0); // TODO change the handling of empty-file case

	/* calculate the number of bytes-to-read per thread */
	number_of_bytes_partitioned = total_number_of_bytes / THREADS;
	buffer_overlap = number_of_bytes_partitioned*0.2; // arbitrary multiplier
	max_number_of_bytes = number_of_bytes_partitioned + total_number_of_bytes % THREADS;
	my_number_of_bytes = (MYTHREAD == THREADS-1)? max_number_of_bytes
			: number_of_bytes_partitioned + buffer_overlap;
	//cout << "total_number_of_bytes=" << total_number_of_bytes << ", number_of_bytes_partitioned=" << number_of_bytes_partitioned
	//		<< ", buffer_overlap=" << buffer_overlap << ", max_number_of_bytes=" << max_number_of_bytes << ", my_number_of_bytes=" << my_number_of_bytes << endl;

	/* allocate read buffers */
	if (my_number_of_bytes >= INT_MAX) {
		DIE("Not enough per-thread memory to read the input file: %s. Consider running on more nodes.\n", infilename);
	}
	read_buffer = (char*) malloc(my_number_of_bytes);
	ASSERT(read_buffer != NULL, "read_buffer is null!");
	//cout << "after malloc, read_buffer=" << (void*) read_buffer << endl; //debugging
	memset(read_buffer, 0, sizeof(char)*my_number_of_bytes);
	DBG("allocated %d bytes\n", my_number_of_bytes);

	/* calculate offset and begin reading */
    my_offset = MYTHREAD * number_of_bytes_partitioned;
    DBG("my offset = %lld\n", my_offset);
	CHECK_MPI(MPI_File_read_at_all(filehandle, my_offset, read_buffer, (int) my_number_of_bytes, MPI_CHAR, &status));
	MPI_Get_count(&status, MPI_CHAR, &rcv_count);
	//cout << "MPI received " << rcv_count << " chars" << endl; //debugging
	ASSERT(rcv_count > 0, "MPI read 0 characters from input file!");

	/* advance start_ptr to the character following the next newline */
	start_ptr = read_buffer;
	if (MYTHREAD > 0 && *(start_ptr-1) != '\n') {
		while(*start_ptr != '\n') start_ptr++;
		start_ptr++;
	}
	DBG("advanced start_ptr %d charcters \n", (start_ptr-read_buffer) );

	/* advance end_ptr to the next newline or EOF */
	if (MYTHREAD == THREADS-1) { end_ptr = read_buffer + my_number_of_bytes; }
	else {
		end_ptr = read_buffer + number_of_bytes_partitioned;
		while(*end_ptr != '\n' && end_ptr < (read_buffer + my_number_of_bytes)) { end_ptr++; }
		if (*end_ptr != '\n') {
			DIE("Cannot find end-of-line within %lld characters of original end-point. " \
					"Consider increasing buffer overlap or whether file may be corrupt.\n", buffer_overlap);
		}
	}
	//*end_ptr='\0'; //TODO this might be buggy
	DBG("advanced end_ptr %d charcters \n", ( end_ptr - (read_buffer + number_of_bytes_partitioned) ) );

	// for each character in char* read_buffer between start_ptr and end_ptr (not \0)
	// seek the first space, the characters before that are the first readId
	// seek the second space, the characters before that are the second readId
	// seek the third space, characters before that are the first read name
	// seek the fourth space, characters before that are the second read name
	// after that, each space separated series of characters are positions
	// the newline is the last character for a readid pair list
	// don't advance read_buffer, want to free it later

	// debugging
	int num_null = 0;
	while (!isdigit(*start_ptr)) {
		start_ptr++;
		num_null++;
	}
	DBG("Rank %d: number of non-digits preceding the input overlaps: %d", MYTHREAD, num_null);

	/*
	 * BUG
	 * for some reason, read_buffer doesn't contain anything, start_ptr and read_buffer are not null but...
	 */
	//cout << "start_ptr=" << (void*) start_ptr << ", end_ptr-start_ptr="<< (end_ptr-start_ptr) <<endl;
	//cout << ", read_buffer is null? "<< (read_buffer == NULL? "true" : "false") << endl;
	//cout << "loop condition is " << ((start_ptr && (end_ptr-start_ptr) >= 2)? "true" : "false") << endl;
	//cout << "start_ptr ? " << (start_ptr? "true" : "false") << endl;
	while(start_ptr && (end_ptr-start_ptr) >= 2) { // no more read pairs possible within 2 (actually 7) characters of the end
		// parse the ReadId pair
		currPair.readId1 = (ReadId) strtoull(start_ptr, &start_ptr, 10);
		DBG("extracted ReadId=%lld, start_ptr-read_buffer=%d, end_ptr-start_ptr=%d \n",
				currPair.readId1, (start_ptr-read_buffer), (end_ptr-start_ptr) );
		assert(errno != EINVAL);
		assert(errno != ERANGE);
		currPair.readId2 = (ReadId) strtoull(start_ptr, &start_ptr, 10);
		DBG("extracted ReadId=%lld, start_ptr-read_buffer=%d, end_ptr-start_ptr=%d \n",
				currPair.readId2, (start_ptr-read_buffer), (end_ptr-start_ptr) );
		assert(errno != EINVAL);
		assert(errno != ERANGE);

		// parse the corresponding name pair
		name1 = strtok_r(start_ptr, delim, &temp);
		start_ptr = temp;
		name2 = strtok_r(start_ptr, delim, &temp);
		start_ptr = temp;

		DBG("Rank %d: name1=%s, name2=%s\n", MYTHREAD, name1, name2);

		// add the ReadId pair to the appropriate buffers
		checkLocality(currPair, name1, name2, num_nonlocal_readids);

		// parse-out all the k-mer positions shared between the ReadId pair (at least 1)
		do {
			currPair.posInRead1 = (PosInRead) strtoul(start_ptr, &start_ptr, 10);
			DBG("extracted position=%hu, start_ptr-read_buffer=%d, end_ptr-start_ptr=%d \n",
					currPair.posInRead1, (start_ptr-read_buffer), (end_ptr-start_ptr) );
			assert(errno != EINVAL);
			assert(errno != ERANGE);
			currPair.posInRead2 = (PosInRead) strtoul(start_ptr, &start_ptr, 10);
			DBG("extracted position=%hu, start_ptr-read_buffer=%d, end_ptr-start_ptr=%d \n",
					currPair.posInRead2, (start_ptr-read_buffer), (end_ptr-start_ptr) );
			assert(errno != EINVAL);
			assert(errno != ERANGE);
			std::pair<PosInRead,PosInRead> val = std::make_pair(currPair.posInRead1, currPair.posInRead2);
			shared_positions.push_back(val);
		} while (start_ptr < end_ptr && *start_ptr != '\n' && *start_ptr != EOF);
		std::pair<ReadId,ReadId> key = std::make_pair(currPair.readId1, currPair.readId2);
		(*readPairs)[key].insert((*readPairs)[key].end(), shared_positions.begin(), shared_positions.end());
		shared_positions.clear();
		DBG("about to start newline or exit: (*start_ptr)=%c, (*start_ptr+1)=%c \n", (*start_ptr), (*start_ptr+1));
		start_ptr++;
	}

	/* clean-up */
	free(read_buffer);
	CHECK_MPI(MPI_File_close(&filehandle));

	LOGF("Number of non-local read ids is %d, in %d total read pairs \n", num_nonlocal_readids, readPairs->size());
	return 0;
}

void Aligner::exchangeSequences() {
	// performance statistics
	double start_time=0.0, tot_packing_time=0.0, tot_exchange_time=0.0, tot_store_time=0.0;
	double tot_glob_packing_time=0.0, tot_glob_store_time=0.0;
	double avg_glob_packing_time=0.0, avg_glob_store_time=0.0;
	double tot_glob_exchange_time=0.0, max_glob_exchange_time=0.0, min_glob_exchange_time=0.0, avg_glob_exchange_time=0.0;
	int exch_itr = -1;
	int nonlocal_read_cnt = 0;
	/*
	 * The maximum number of characters that can be exchanged per thread per round,
	 * due to MPI's limitation on counts and displacements (int representation).
	 * Assumes there are no more than 2^16 characters per sequence (uncompressed)
	 * -- current long read sequencing limitation. This may need to be updated soon
	 * as some of the latest data sets report sequences over 100,000 (with ~10,000 on average).
	 */
	int MAX_READS_PER_EXCHANGE = INT32_MAX/UINT16_MAX/THREADS;
	ASSERT(MAX_READS_PER_EXCHANGE > 0, "MAX_READS_PER_EXCHANGE is <= 0. Need to refactor sequence exchange...");

	int tot_readids_wanted, // restricts displacement indices to max handleable by MPI Alltoallv
		tot_cnt_owned,
		tot_chars_owned,
		tot_chars_rcvd;
	char* chars_index;
	// initialize per thread buffers
	// naming convention:
	// 		"wanted" == this processor does not own and wants
	//		"owned"  == this processor owns and another processors wants
	// (I know, it's sad when spelling out a naming convention is necessary...)
	int* cnts_for_wanted;
	int* displs_for_wanted_ids;
	int* cnts_for_owned;
	int* cnts_tobe_owned;
	int* displs_for_owned_ids;
	int* displs_for_owned_seqs;
	int* cnts_seqs_owned;
	int* dspls_seqs_owned;

	cnts_for_wanted = new int[THREADS]();		ASSERT(cnts_for_wanted != NULL,""); // TODO handle
	displs_for_wanted_ids = new int[THREADS]();	ASSERT(displs_for_wanted_ids != NULL,""); // TODO handle
	cnts_for_owned = new int[THREADS]();		ASSERT(cnts_for_owned != NULL,""); // TODO handle
	cnts_tobe_owned = new int[THREADS]();		ASSERT(cnts_for_owned != NULL,""); // TODO handle
	displs_for_owned_ids = new int[THREADS]();	ASSERT(displs_for_owned_ids != NULL,""); // TODO handle
	displs_for_owned_seqs = new int[THREADS]();	ASSERT(displs_for_owned_ids != NULL,""); // TODO handle
	cnts_seqs_owned = new int[THREADS]();		ASSERT(cnts_seqs_owned != NULL,""); // TODO handle
	dspls_seqs_owned = new int[THREADS]();		ASSERT(dspls_seqs_owned != NULL,""); // TODO handle
	size_t bytes_per_buffer = THREADS*sizeof(int);

	Buffer buffRidsWanted = initBuffer(MAX_ALLTOALL_MEM);
	Buffer buffRidsOwned = initBuffer(MAX_ALLTOALL_MEM);
	Buffer buffCharsSnd = initBuffer(MAX_ALLTOALL_MEM);
	Buffer buffCharsRcv = initBuffer(MAX_ALLTOALL_MEM);
	Buffer buffLengthsSnd = initBuffer(MAX_ALLTOALL_MEM);
	Buffer buffLengthsRcv = initBuffer(MAX_ALLTOALL_MEM);

	ReadId* wantedReadIds;
	ReadId* ownedReadIds;
	char* seqs_owned;
	char* seqs_names_owned;
	char* rcvseqs;
	char* rcvnames;
	PosInRead* sndseqlengths;
	PosInRead* rcvseqlengths;

	bool moreToSnd=false,
		 moreToRcv=false;

	/*
	 * Sends a maximum of INT8_MAX Read ID's (with associated sequences and name tags) per thread, per loop.
	 * Technically, we could send up to UINT16_MAX without an overflow error -- since individual sequence lengths
	 * are up to UINT16_MAX, this would ensure that we send less than (2^15 * 2^16 = 2^31) characters each loop.
	 * However, in practice, the communication buffers can handle much less than that maximum in a single Alltoallv.
	 * We could also send more ReadId's upfront, with an additional inner loop for the actual sequences,
	 * but simplicity is nice sometimes...
	 */
	do { // while moreToSnd || moreToRcv
		exch_itr++;
		moreToSnd=false,
		moreToRcv=false;

		memset(cnts_for_wanted, 0, bytes_per_buffer);
		memset(displs_for_wanted_ids, 0, bytes_per_buffer);
		memset(cnts_for_owned, 0, bytes_per_buffer);
		memset(cnts_tobe_owned, 0, bytes_per_buffer);
		memset(displs_for_owned_ids, 0, bytes_per_buffer);
		memset(displs_for_owned_seqs, 0, bytes_per_buffer);
		memset(cnts_seqs_owned, 0, bytes_per_buffer);
		memset(dspls_seqs_owned, 0, bytes_per_buffer);

		resetBuffer(buffRidsWanted);
		resetBuffer(buffRidsOwned);
		resetBuffer(buffCharsSnd);
		resetBuffer(buffCharsRcv);
		resetBuffer(buffLengthsSnd);
		resetBuffer(buffLengthsRcv);

		//////////////////////////// Begin ReadId exchange ///////////////////////////////

		//
		// initialize buffer for requested number of ReadIds from each other processor
		//
		start_time = MPI_Wtime();
		tot_readids_wanted = 0;
		for (int i = 0; i < THREADS; i++) {
			displs_for_wanted_ids[i] = tot_readids_wanted;
			cnts_for_wanted[i] = min(perProcRequestIds[i].size(), (size_t) MAX_READS_PER_EXCHANGE );
			tot_readids_wanted += cnts_for_wanted[i];
			if (cnts_for_wanted[i] < 0 || tot_readids_wanted < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			DBG("Asking for %d reads from rank %d \n", cnts_for_wanted[i], i);
		}
		tot_packing_time += MPI_Wtime() - start_time;

		// (not timed)
		// calculate and report min, max, avg, and total number of global reads requested
		//
		CHECK_MPI( MPI_Reduce(&tot_readids_wanted, &nonlocal_read_cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) );
		if (nonlocal_read_cnt < 0) {
			serial_printf("Aligner:%s:iteration %d: (overflow error) %d total reads requested\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
		}
		else {
			serial_printf("Aligner:%s:iteration %d: %d total reads requested\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
			serial_printf("Aligner:%s:iteration %d: %d average reads requested per rank from all other ranks\n", __FUNCTION__, exch_itr, nonlocal_read_cnt/THREADS);
		}
		nonlocal_read_cnt=0;
		CHECK_MPI( MPI_Reduce(&tot_readids_wanted, &nonlocal_read_cnt, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD) );
		serial_printf("Aligner:%s:iteration %d: %d min reads requested\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
		nonlocal_read_cnt=0;
		CHECK_MPI( MPI_Reduce(&tot_readids_wanted, &nonlocal_read_cnt, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD) );
		serial_printf("Aligner:%s:iteration %d: %d max reads requested\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
		nonlocal_read_cnt=0;

		//
		// exchange number of ReadIds requested from/for each processor
		//
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoall(cnts_for_wanted, 1, MPI_INT, cnts_for_owned, 1, MPI_INT, MPI_COMM_WORLD) );
		tot_exchange_time += MPI_Wtime() - start_time;

		//
		// initialize buffers of ReadIds to send (requesting respective reads from other processors)
		//
		start_time = MPI_Wtime();
		growBuffer(buffRidsWanted, tot_readids_wanted * sizeof(ReadId));
		wantedReadIds = (ReadId*) getStartBuffer(buffRidsWanted);
		DBG("successfully allocated send buffer of %d ReadIds, buffer size=%lld, buffer len=%lld \n", tot_readids_wanted, buffRidsWanted->size, buffRidsWanted->len);
		for (int i = 0; i < THREADS; i++) {
			int index = 0;
			auto itr = perProcRequestIds[i].begin();
			for (itr ; itr != perProcRequestIds[i].end() && index < cnts_for_wanted[i]; itr++)
			{
				wantedReadIds[displs_for_wanted_ids[i] + index] = (ReadId) *itr;
				index++;
			}
			ASSERT(index == cnts_for_wanted[i],"cnts_for_wanted["+to_string(i)+"]="+to_string(cnts_for_wanted[i])+", but index="+to_string(index));
			perProcRequestIds[i].erase(perProcRequestIds[i].begin(), itr); // erases to the current itr location
			moreToSnd = (!perProcRequestIds[i].empty()? true : moreToSnd);
			DBG("erased perProcRequestIds[%d].size() is %zu \n", i, perProcRequestIds[i].size());
		}

		//
		// initialize receive and receive displacement buffers from previous alltoall
		//
		tot_cnt_owned = 0;
		for (int i = 0; i < THREADS; i++) {
			displs_for_owned_ids[i] = tot_cnt_owned;
			tot_cnt_owned += cnts_for_owned[i];
			if (cnts_for_owned[i] < 0 || tot_cnt_owned < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			DBG("Expecting %d requested ReadIds from rank %d at position rdispls[%d]=%d \n",cnts_for_owned[i], i, i, displs_for_owned_ids[i]);
		}
		growBuffer(buffRidsOwned, tot_cnt_owned * sizeof(ReadId));
		ownedReadIds = (ReadId*) getStartBuffer(buffRidsOwned);
		DBG("Successfully initialized receive buffer of %d ReadIds, buffer size=%lld, buffer len=%lld \n",tot_cnt_owned, buffRidsOwned->size, buffRidsOwned->len);
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange ReadIds of requested sequences
		//
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoallv(wantedReadIds, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT64_T,
				ownedReadIds, cnts_for_owned, displs_for_owned_ids, MPI_UINT64_T, MPI_COMM_WORLD) );
		DBG("successfully exchanged ReadIds \n");
		tot_exchange_time += MPI_Wtime() - start_time;

		//////////////////////////// end ReadId exchange ///////////////////////////////

		//
		// for each processor's requested reads by ReadId, determine the lengths of the respective reads
		//
		start_time = MPI_Wtime();
		tot_chars_owned = 0;
		size_t seqlen;
		memset(cnts_seqs_owned, 0, sizeof(int) * THREADS);
		memset(dspls_seqs_owned, 0, sizeof(int) * THREADS);
		growBuffer(buffLengthsSnd, tot_cnt_owned * sizeof(PosInRead) );
		sndseqlengths = (PosInRead*) getStartBuffer(buffLengthsSnd);
		DBG("Successfully initialized sndseqlengths of %d PosInRead, buffer size=%lld, buffer len=%lld \n", tot_cnt_owned, buffLengthsSnd->size, buffLengthsSnd->len);
		for (int i = 0; i < THREADS; i++) {
			dspls_seqs_owned[i] += tot_chars_owned;
			for (int j = displs_for_owned_ids[i]; j < (i==THREADS-1? tot_cnt_owned : displs_for_owned_ids[i+1]); j++) { // j should index rcvbuffer directly
				ReadId currId = ownedReadIds[j];
				ASSERT(readMap->count(currId) > 0,"Rank "+to_string(MYTHREAD)+": was asked for global ReadId="+to_string(currId)
						+", but it does not exist in readMap. global offset is "+to_string(globalReadIdOffset)
						+". last ReadId is "+to_string(readRangeEnds[MYTHREAD])); // if the ID was requested at all, it should be on this processor
				seqlen = ( (*readMap)[currId] ).length();
				//DBG("sequence-to-send %lld has length %d \n", currId, seqlen);
				sndseqlengths[j] = seqlen;
				cnts_seqs_owned[i] += seqlen;
			}
			tot_chars_owned += cnts_seqs_owned[i];
			if (cnts_seqs_owned[i] < 0 || tot_chars_owned < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		DBG("finished calculating read lengths, total %d chars \n", tot_chars_owned);
		growBuffer(buffLengthsRcv, tot_readids_wanted*sizeof(PosInRead)); // one length for each ReadId I requested
		rcvseqlengths = (PosInRead*) getStartBuffer(buffLengthsRcv);
		DBG("Successfully initialized rcvseqlengths of %d PosInRead, buffer size=%lld, buffer len=%lld \n", tot_readids_wanted, buffLengthsRcv->size, buffLengthsRcv->len);
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange sequence lengths
		//
		start_time = MPI_Wtime();
		// cnts and displacements stay the same since the last exchange because we're exchanging 1 length per ReadId
#ifdef TIGHT_READS
		CHECK_MPI( MPI_Alltoallv(sndseqlengths, cnts_for_owned, displs_for_owned_ids, MPI_UINT16_T,
				rcvseqlengths, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT16_T, MPI_COMM_WORLD) );
#else
		CHECK_MPI( MPI_Alltoallv(sndseqlengths, cnts_for_owned, displs_for_owned_ids, MPI_UINT32_T,
				rcvseqlengths, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT32_T, MPI_COMM_WORLD) );
#endif
		tot_exchange_time += MPI_Wtime() - start_time;
		LOGF("finished exchanging sequence lengths in %.2f s", (MPI_Wtime() - start_time) );

		// (not timed)
		// calculate and report min and max number of read bytes to be sent
		//
		CHECK_MPI( MPI_Reduce(&tot_chars_owned, &nonlocal_read_cnt, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD) );
		serial_printf("Aligner:%s:iteration %d: sending %d min read bytes per rank\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
		nonlocal_read_cnt=0;
		CHECK_MPI( MPI_Reduce(&tot_chars_owned, &nonlocal_read_cnt, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD) );
		serial_printf("Aligner:%s:iteration %d: sending %d max read bytes per rank\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
		nonlocal_read_cnt=0;

		//
		// copy actual sequences into send buffer
		//
		start_time = MPI_Wtime();
		growBuffer(buffCharsSnd, tot_chars_owned * sizeof(char));
		seqs_owned = (char*) getStartBuffer(buffCharsSnd);
		DBG("allocated seqs_owned buffer in %.2f s\n", (MPI_Wtime() - start_time));
		chars_index = seqs_owned;
		std::string val;
		for (int i = 0; i < tot_cnt_owned; i++) {
			val = (*readMap)[ ownedReadIds[i] ];
			val.copy( chars_index, sndseqlengths[i]);
			chars_index += sndseqlengths[i];
		}
		DBG("finished filling send buffers with sequences in %.2f s\n", (MPI_Wtime() - start_time));

		//
		// init rcv buffers based on number of chars expected from each processor
		//
		tot_chars_rcvd = 0;
		for (int i = 0; i < THREADS; i++) {
			displs_for_owned_seqs[i] = tot_chars_rcvd;
			for (int j = displs_for_wanted_ids[i]; j < (displs_for_wanted_ids[i] + cnts_for_wanted[i]); j++) {
				cnts_tobe_owned[i] += rcvseqlengths[j];
				tot_chars_rcvd += rcvseqlengths[j];
			}
			if (tot_chars_rcvd < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		growBuffer(buffCharsRcv, tot_chars_rcvd * sizeof(char));
		rcvseqs = (char*) getStartBuffer(buffCharsRcv);
		tot_packing_time += MPI_Wtime() - start_time;
		DBG("successfully initialized receive buffers in %.2f s \n", (MPI_Wtime() - start_time));

		// (not timed)
		// calculate and report min, max, avg, and total number of read bytes to be exchanged
		//
		CHECK_MPI( MPI_Reduce(&tot_chars_rcvd, &nonlocal_read_cnt, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD) );
		serial_printf("Aligner:%s:iteration %d: receiving %d min read bytes per rank\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
		nonlocal_read_cnt=0;
		CHECK_MPI( MPI_Reduce(&tot_chars_rcvd, &nonlocal_read_cnt, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD) );
		serial_printf("Aligner:%s:iteration %d: receiving %d max read bytes per rank\n", __FUNCTION__, exch_itr, nonlocal_read_cnt);
		nonlocal_read_cnt=0;

		//
		// exchange sequence data
		//
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoallv(seqs_owned, cnts_seqs_owned, dspls_seqs_owned, MPI_CHAR,
				rcvseqs, cnts_tobe_owned, displs_for_owned_seqs, MPI_CHAR, MPI_COMM_WORLD) );
		tot_exchange_time += MPI_Wtime() - start_time;
		DBG("finished exchanging sequence data in %.2f s, cleaning up \n", (MPI_Wtime() - start_time));

		//
		// store received sequence data in readMap
		//
		start_time = MPI_Wtime();
		std::string* tostore;
		chars_index = rcvseqs; // chars_index reuse
		DBG("storing received sequences \n");
		for (int i = 0; i < tot_readids_wanted; i++) { // for each sequence received
			DBG("tot_wanted_readids=%d, chars_index=%p, i=%d, rcvseqlengths[i]=%d for ReadId=%lld \n", tot_readids_wanted,
					chars_index, i, rcvseqlengths[i], wantedReadIds[i]);
			tostore = new std::string(chars_index, rcvseqlengths[i]);
			chars_index += rcvseqlengths[i];
			(*readMap)[ wantedReadIds[i] ] = *tostore;
		}
		tot_store_time += MPI_Wtime() - start_time;
		DBG("Finished storing local sequences in %.2f s, doing intermediate clean-up \n", (MPI_Wtime() - start_time));
		CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

		// reset for reuse within same iteration
		resetBuffer(buffCharsSnd);
		resetBuffer(buffCharsRcv);

		//////////////////////////// BEGIN NAME TAG EXCHANGE ////////////////////////////
		start_time = MPI_Wtime();
		std::memset(sndseqlengths, 0, sizeof(PosInRead) * tot_cnt_owned);
		std::memset(cnts_seqs_owned, 0, sizeof(int) * THREADS);
		std::memset(dspls_seqs_owned, 0, sizeof(int) * THREADS);
		std::memset(rcvseqlengths, 0, sizeof(PosInRead) * tot_readids_wanted);

		//
		// for each name requested by ReadID, calculate the sum and displacements per thread
		//
		tot_chars_owned = 0;
		for (int i = 0; i < THREADS; i++) {
			dspls_seqs_owned[i] += tot_chars_owned;
			for (int j = displs_for_owned_ids[i]; j < (i==THREADS-1? tot_cnt_owned : displs_for_owned_ids[i+1]); j++) { // j should index rcvbuffer directly
				ReadId currId = ownedReadIds[j];
				ASSERT(readNameMap->count(currId) > 0,"Rank "+to_string(MYTHREAD)+": was asked for global ReadId="+to_string(currId)
						+", but it does not exist in readMap. global offset is "+to_string(globalReadIdOffset)
						+". last ReadId is "+to_string(readRangeEnds[MYTHREAD])); // if the ID was requested at all, it should be on this processor
				seqlen = ( (*readNameMap)[currId] ).length();
				DBG("name-to-send %lld has length %d \n", currId, seqlen);
				sndseqlengths[j] = seqlen;
				cnts_seqs_owned[i] += seqlen;
			}
			tot_chars_owned += cnts_seqs_owned[i];
			if (tot_chars_owned < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		DBG("finished calculating name lengths, total %d chars \n", tot_chars_owned);
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange name lengths and displacements
		//
		start_time = MPI_Wtime();
		// cnts and displacements stay the same since the last exchange because we're exchanging 1 length per ReadId
		//verified sndseqlengths, cnts_seqs_owned, displs_for_owned_ids, rcvseqlengths, cnts_for_wanted, displs_for_wanted
#ifdef TIGHT_READS
		CHECK_MPI( MPI_Alltoallv(sndseqlengths, cnts_for_owned, displs_for_owned_ids, MPI_UINT16_T,
				rcvseqlengths, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT16_T, MPI_COMM_WORLD) );
#else
		CHECK_MPI( MPI_Alltoallv(sndseqlengths, cnts_for_owned, displs_for_owned_ids, MPI_UINT32_T,
				rcvseqlengths, cnts_for_wanted, displs_for_wanted_ids, MPI_UINT32_T, MPI_COMM_WORLD) );
#endif
		tot_exchange_time += MPI_Wtime() - start_time;

		//
		// initialize receive buffers and package names for exchange
		//
		start_time = MPI_Wtime();
		growBuffer(buffCharsSnd, tot_chars_owned * sizeof(char));
		seqs_names_owned = (char*) getStartBuffer(buffCharsSnd);
		chars_index = seqs_names_owned;
		val = "";
		for (int i = 0; i < tot_cnt_owned; i++) {
			val = (*readNameMap)[ ownedReadIds[i] ];
			val.copy( chars_index, sndseqlengths[i]);
			chars_index += sndseqlengths[i];
		}
		DBG("finished filling send buffers with sequences \n");

		//
		// prep for reuse,
		// owned cnts and displs to be used as "to be owned" buffers
		//
		std::memset(cnts_for_owned, 0, sizeof(int) * THREADS);
		std::memset(displs_for_owned_seqs, 0, sizeof(int) * THREADS);

		//
		// init rcv buffers based on number of chars expected from each processor
		//
		tot_chars_rcvd = 0;
		for (int i = 0; i < THREADS; i++) {
			displs_for_owned_seqs[i] = tot_chars_rcvd;
			for (int j = displs_for_wanted_ids[i]; j < (displs_for_wanted_ids[i] + cnts_for_wanted[i]); j++) {
				cnts_for_owned[i] += rcvseqlengths[j];
				tot_chars_rcvd += rcvseqlengths[j];
			}
			if (tot_chars_rcvd < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}

		growBuffer(buffCharsRcv, tot_chars_rcvd * sizeof(char));
		rcvnames = (char*) getStartBuffer(buffCharsRcv);
		DBG("successfully initialized receive buffers \n");
		tot_packing_time += MPI_Wtime() - start_time;

		//
		// exchange actual names
		//
		start_time = MPI_Wtime();
		CHECK_MPI( MPI_Alltoallv(seqs_names_owned, cnts_seqs_owned, dspls_seqs_owned, MPI_CHAR,
				rcvnames, cnts_for_owned, displs_for_owned_seqs, MPI_CHAR, MPI_COMM_WORLD) );
		tot_exchange_time += MPI_Wtime() - start_time;
		DBG("exchanged sequence data, cleaning up \n");

		//
		// check received names match names loaded from overlap file (no inconsistent numbering/tagging)
		//
		// *originally, the tags weren't provided in the overlap file and were retrieved along with the sequences
		//
		start_time = MPI_Wtime();
		tostore = NULL;
		chars_index = rcvnames; // chars_index reuse
		DBG("storing received names \n");
		for (int i = 0; i < tot_readids_wanted; i++) { // for each sequence received
			if(rcvseqlengths[i] < 0) {
				cerr << "Rank " << MYTHREAD << ": ["<<__FILE__ << " " << __LINE__ << "] detected overflow in sequenceExchange " << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			DBG("tot_readids_wanted=%d, chars_index=%p, i=%d, rcvseqlengths[i]=%d for ReadId=%lld \n", tot_readids_wanted,
					chars_index, i, rcvseqlengths[i], wantedReadIds[i]);
			tostore = new std::string(chars_index, rcvseqlengths[i]);
			chars_index += rcvseqlengths[i];
			ASSERT ( (*readNameMap)[ wantedReadIds[i] ] == *tostore, "Rank " << MYTHREAD << ": locally stored name for read ID=" << wantedReadIds[i] << ", does not match remote name, " << (*readNameMap)[ wantedReadIds[i] ] << " != " << *tostore);
			//(*readNameMap)[ wantedReadIds[i] ] = *tostore; // originally, the tags weren't provided in the overlap file and were retrieved along with the sequences, now we just check they match
		}
		DBG("Finished checking local sequence names, doing final clean-up \n");
		tot_store_time += MPI_Wtime() - start_time;
		CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

		//////////////////////////// END NAME TAG EXCHANGE ////////////////////////////
		DBG("has more to send? %s \n", (moreToSnd? "yes" : "no"));
		CHECK_MPI( MPI_Allreduce(&moreToSnd, &moreToRcv, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD));
		DBG("has more to receive? %s \n", (moreToRcv? "yes" : "no"));
	} while (moreToSnd || moreToRcv);

	//
	// do final cleanup
	//
	delete[] cnts_for_owned;
	delete[] cnts_tobe_owned;
	delete[] displs_for_owned_ids;
	delete[] displs_for_owned_seqs;
	delete[] cnts_for_wanted;
	delete[] displs_for_wanted_ids;
	delete[] cnts_seqs_owned;
	delete[] dspls_seqs_owned;
	freeBuffer(buffRidsWanted);
	freeBuffer(buffRidsOwned);
	freeBuffer(buffCharsSnd);
	freeBuffer(buffCharsRcv);
	freeBuffer(buffLengthsSnd);
	freeBuffer(buffLengthsRcv);

	//
	// output performance statistics
	//
	CHECK_MPI( MPI_Reduce(&tot_store_time, &tot_glob_store_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_store_time = tot_glob_store_time/THREADS;

	CHECK_MPI( MPI_Reduce(&tot_packing_time, &tot_glob_packing_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_packing_time = tot_glob_packing_time/THREADS;

	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &min_glob_exchange_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &max_glob_exchange_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&tot_exchange_time, &tot_glob_exchange_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	avg_glob_exchange_time = tot_glob_exchange_time/THREADS;

	serial_printf("Aligner:%s Finished sequence exchange \n", __FUNCTION__);
	serial_printf("Aligner:%s Sequence Exchange: total (sum) and average (packing) times were %.2f s and %.2f s \n",
		__FUNCTION__, tot_glob_packing_time, avg_glob_packing_time);
	serial_printf("Aligner:%s Sequence Exchange: total (sum) and average (store) times were %.2f s and %.2f s \n",
		__FUNCTION__, tot_glob_store_time, avg_glob_store_time);
	serial_printf("Aligner:%s Sequence Exchange: total (sum), average, min, and max (exchange) times (s) were: %.2f, %.2f, %.2f, %.2f \n",
		__FUNCTION__, tot_glob_exchange_time, avg_glob_exchange_time, min_glob_exchange_time, max_glob_exchange_time);
}

/*
 * Assumes: sequences indexable using readPairMap are stored in Aligner::readMap
 */
void Aligner::computeLocalAlignments
	(ReadOverlapMap* overlaplist, std::vector<alignment_data>& allAlignments, std::string outfilename)
{
	// performance statistics
	unsigned long long num_overlaps = overlaplist->size();
	unsigned long long global_num_overlaps, alignments_attempted = 0;
	long glob_min_algnmnts=0, glob_max_algnmnts=0, glob_avg_algnmnts=0, glob_tot_algnmnts=0;
	double agg_seqan_time=0.0, agg_io_time=0.0, algn_elapsed_time=0.0;
	double start_time=0.0, intrm_start_time=0.0;
	double glob_avg_time=0.0, glob_min_time=0.0, glob_max_time=0.0, glob_tot_time=0.0;
	LOGF("About to compute alignments for %d read pairs. \n", num_overlaps);

	std::string readA, readB;
	ReadId readIdA, readIdB;
	PosInRead posInA, posInB;
	bool do_align;
	int align_len, editDist, score;
	int num_failed=0;
	int num_readpairs_aligned = 0;
	int64_t file_offset = 0;
	start_time = MPI_Wtime();
	// for each candidate overlap (overlapping pair of reads)
	for (auto itr = overlaplist->begin(); itr != overlaplist->end(); itr++) {
		readIdA = (itr->first).first;
		readIdB = (itr->first).second;
		ASSERT( readMap->count(readIdA) > 0, "globalId="+to_string(readIdA));
		ASSERT( readMap->count(readIdB) > 0, "globalId="+to_string(readIdB));
		readA = (*readMap)[ readIdA ];	ASSERT(!readA.empty() && readA.length()>=kmerLen,"");
		readB = (*readMap)[ readIdB ];	ASSERT(!readB.empty() && readB.length()>=kmerLen,"");
		int end = std::min((itr->second).size(), (size_t) maxPerPairAlignments);
		alignment_data* keepData = NULL;
		int best_score;
		auto pos_itr = (itr->second).begin();
		ASSERT(pos_itr != (itr->second).end(),"");
		// for each seed pair
		for (int i = 0; i < end && pos_itr != (itr->second).end(); pos_itr++) {
			posInA = pos_itr->first;
			ASSERT(posInA+kmerLen <= readA.length() && posInA >= 0,"Rank "+to_string(MYTHREAD)+": readIdA="
					+to_string(readIdA)+", posInA="+to_string(posInA)+", kmerLen="+to_string(kmerLen)
					+", readA.length()="+to_string(readA.length())+" readIdB=" +to_string(readIdB)
					+"posInB="+to_string(posInB)+", readB.length()="+to_string(readB.length()));
			posInB = pos_itr->second;
			ASSERT(posInB+kmerLen <= readB.length() && posInB >= 0,"Rank "+to_string(MYTHREAD)+": readIdB="
					+to_string(readIdB)+"posInB="+to_string(posInB)+", kmerLen="+to_string(kmerLen)
					+", readB.length()="+to_string(readB.length()));
			DBG("about to compute alignment for [%lld, %lld, %d, %d], length(readA)= %d, length(readB)=%d, k=%d %s \n", readIdA, readIdB,
				posInA, posInB, readA.length(), readB.length(), kmerLen, (_skip_algnmnt_krnl? ", but _skip_alignment set" : ""));
			if (!_skip_algnmnt_krnl) {
	#ifdef SEQAN
				alignment_data currData(readIdA, readIdB);
				agg_seqan_time += seqan_seedExtend(readA, readB, posInA, posInB, kmerLen, currData);
				bool keep = false;
				if (myscoring.adaptive) {
					// directly adapted from BELLA code
					int diffCol = currData.r1_al_end - currData.r1_al_start;
					int diffRow = currData.r2_al_end - currData.r2_al_start;
					int minLeft = min(currData.r1_al_start, currData.r2_al_start);
					int minRight = min(currData.length_r1 - currData.r1_al_end, currData.length_r2 - currData.r2_al_end);
					int overlap = minLeft+minRight+(diffCol+diffRow)/2;

					double threshold = (1 - myscoring.chernoff_delta) * (myscoring.adaptive_phi * overlap);
					keep = myscoring.aligntoend ?
							(currData.score >= threshold) && aligned_to_end(currData)
							: (currData.score >= threshold);
				}
				else { keep = currData.score >= myscoring.static_threshold; }
				if (keep && (keepData == NULL || currData.score > best_score)) {
					keepData = &currData;
					best_score = currData.score;
				}
	#endif
	#ifdef EDLIB // TODO get rid of edlib stuff
				do_align = edlibOp((*itr), align_len, editDist);
				if (do_align) { //store alignments
					allAlignments.push_back( alignment_data (  *itr, align_len, editDist) );
				}
	#endif
			}
			i++;
			alignments_attempted++;
		}
		if(keepData != NULL) {
			num_readpairs_aligned++;
			keepData->similarity = (itr->second).size();
			allAlignments.push_back( *keepData );
// TODO outputAlignmentsToFileWithMpi is collective, so if we write the output in batches with it, we need to synchronize with other processes, but as it is, they have different numbers of alignments to do
//			if (num_readpairs_aligned % OUT_BATCH_SIZE == 0) {
//#ifndef BENCHMARKONLY
//				intrm_start_time = MPI_Wtime();
//				file_offset = outputAlignmentsToFileWithMpi(&allAlignments, outfilename, file_offset);
//				agg_io_time += (MPI_Wtime() - intrm_start_time);
//#endif
//				allAlignments.clear();
//			}
		}
	}
	algn_elapsed_time = MPI_Wtime() - start_time;
#ifndef BENCHMARKONLY
	intrm_start_time = MPI_Wtime();
	file_offset = outputAlignmentsToFileWithMpi(&allAlignments, outfilename, file_offset);
	agg_io_time += (MPI_Wtime() - intrm_start_time);
#endif
	allAlignments.clear();

	// log individual process time for later analysis
	if (num_failed > 0) { cout << "Rank " << MYTHREAD << ": failed to compute a total of " << num_failed << " alignments. See thread log for details. " << endl; }
	LOGF("Process total times: local alignment elapsed time for %ld alignment attempts is approx. %0.3f\n", alignments_attempted, agg_seqan_time);

	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	CHECK_MPI( MPI_Reduce( &num_overlaps, &global_num_overlaps, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
	serial_printf("Total global number of overlaps: %lld\n", global_num_overlaps);

	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	CHECK_MPI( MPI_Reduce(&alignments_attempted, &glob_min_algnmnts, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&alignments_attempted, &glob_max_algnmnts, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&alignments_attempted, &glob_tot_algnmnts, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_avg_algnmnts = glob_tot_algnmnts/THREADS;

	CHECK_MPI( MPI_Reduce(&algn_elapsed_time, &glob_min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&algn_elapsed_time, &glob_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&algn_elapsed_time, &glob_tot_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_avg_time = glob_tot_time/THREADS;
	serial_printf("Aligner:%s Local alignment maximum (global elapsed) time for %ld total attempted alignments is %0.3f\n",
				__FUNCTION__, glob_tot_algnmnts, glob_max_time);

	serial_printf("Aligner:%s Min, max, and avg number of local alignments computed: %ld, %ld, %ld \n",
				 __FUNCTION__, glob_min_algnmnts, glob_max_algnmnts, glob_avg_algnmnts);

	CHECK_MPI( MPI_Reduce(&agg_seqan_time, &glob_min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&agg_seqan_time, &glob_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&agg_seqan_time, &glob_tot_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_avg_time = glob_tot_time/THREADS;
	serial_printf("Aligner:%s Min, max, avg and total (sum) time (s) for SeqAn seed and extend: %.3f, %.3f, %.3f, %.3f \n",
				 __FUNCTION__, glob_min_time, glob_max_time, glob_avg_time, glob_tot_time);

	CHECK_MPI( MPI_Reduce(&agg_io_time, &glob_min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&agg_io_time, &glob_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&agg_io_time, &glob_tot_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
	glob_avg_time = glob_tot_time/THREADS;
	serial_printf("Aligner:%s Min, max, avg and total (sum) time (s) for alignment IO: %.3f, %.3f, %.3f, %.3f \n",
				 __FUNCTION__, glob_min_time, glob_max_time, glob_avg_time, glob_tot_time);
}

/* function imported from pipeline-refactor version */
/*
void Aligner::computeAlignmentsFromOverlapsInMemory(bool cachedIO, char* inputFqFiles,
       std::vector<alignment_data>& allAlignments, string outfilename)
{
       // load sequences and sequence names
       std::vector<filenamesize> allfiles = loadAllFiles(inputFqFiles);
       std::vector<std::string>* allreads = new std::vector<std::string>();
       std::vector<std::string>* allnames = new std::vector<std::string>();
       ReadId numReads = loadSequences(allfiles, *allreads, *allnames, cachedIO);
       storeSequences(*allreads, *allnames, globalReadIdOffset+1); // read ID's start at 1...
       delete allreads; // all stored in map now
       delete allnames; // all stored in map now

       // find and buffer any non-local reads by ID
       int num_nonlocal_readids = 0;
       ReadOverlapPair currPair;
       for(auto itr = readPairs->begin(); itr != readPairs->end(); itr++) {
               currPair.readId1 = get<0>(itr->first);
               currPair.readId2 = get<1>(itr->first);
               checkLocality(currPair, num_nonlocal_readids); //TODO checkLocality should be rewritten generically
       }

       // exchange necessary read sequence requests before computing alignments
       exchangeSequences();

       // compute and output all alignments
       computeLocalAlignments(readPairs, allAlignments, outfilename); // compute alignments, timing code embedded
 }
 */

/*
 * main/primary function
 */
void Aligner::computeAlignmentsFromOverlaps(char* inputOverlaps, bool loadperthread, std::string prevAlignments, std::string outfilename,
		char* inputFqFiles, bool cachedFastq, bool cachedOverlaps, const char* base_dir, std::vector<alignment_data>& allAlignments, bool okClearData)
{
#if defined(DEBUG) && defined(NDEBUG)
	if (MYTHREAD==0) { cout << "NDEBUG is defined!" << endl; }
#endif
	double start_time=0.0, elapsed_time=0.0;

	// load sequences from FASTQ file - read ranges must be initialized before loading overlaps
	std::vector<filenamesize> allfiles = loadAllFiles(inputFqFiles);
	initReadRangesAndSequences(allfiles, cachedFastq, base_dir); // timing code embedded
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

	// load previous alignment pair Ids to avoid loading corresponding overlap records and recomputing alignments
	/* there's some weird bug in this optimization TBD
	start_time = MPI_Wtime();
	std::cout << "Rank " << MYTHREAD << ": invoking loadPrevAlignmentSet with string: " << prevAlignments << std::endl;
	bool checkPrev = loadPrevAlignmentSet(prevAlignments);
	std::cout << "Rank " << MYTHREAD << " loadPrevAlignmentSet returned: " << checkPrev << std::endl;
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD)); // debugging
	if ( (!checkPrev) && (prevAlignmentsSet != NULL) ) {
		delete prevAlignmentsSet;
		std::cout << "Rank " << MYTHREAD << " deleted  prevAlignmentsSet " << std::endl;
	} // save the memory ...
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD)); // timing
	elapsed_time = MPI_Wtime() - start_time;
	if (checkPrev) { serial_printf("Loaded previous alignment IDs in %.2f s\n", elapsed_time); }
	*/

	// read-in overlaps from inFile
	start_time = MPI_Wtime();
	if(cachedOverlaps || loadperthread) {
		char overlapfilepath[MAX_FILE_PATH];
		memset(overlapfilepath, '\0', sizeof(char)*MAX_FILE_PATH);
		if (cachedOverlaps) { sprintf(overlapfilepath, VMFILE_PATH); }
		strcat(overlapfilepath, inputOverlaps);
		get_rank_path(overlapfilepath, MYTHREAD);
		DBG("loading overlaps from per thread files, %s, cachedOverlaps=%d.", overlapfilepath, cachedOverlaps);
		loadOverlapsPerThread(overlapfilepath);
	} else {
		DBG("loading overlaps from shared file %s.", inputOverlaps);
		loadOverlapsMpi(inputOverlaps, false);
	}
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
	elapsed_time = MPI_Wtime() - start_time;
	serial_printf("Loaded overlaps in %.2f s\n", elapsed_time);
	LOGF("Loaded %lld overlaps\n", readPairs->size());

	// exchange necessary read sequence requests before computing alignments
	exchangeSequences();

	computeLocalAlignments(readPairs, allAlignments, outfilename); // compute alignments, timing code embedded

	if(okClearData) readPairs->clear();
}

// TODO replace with MPI_I/O
// TODO save this implementation just in case comparison later becomes interesting
int Aligner::outputAlignmentsToFile
	(std::vector<alignment_data>* alignmentsToOutput, std::string outfilename)
{
	std::ofstream myoutputfile;
	myoutputfile.open (outfilename, ios::out | ios::ate | ios::app ); // add alignment data to end of file
	if (!myoutputfile.is_open()) {
		DIE("Could not open %s: %s\n", outfilename, strerror(errno));
		return -1;
	}
	LOGF("successfully opened file %s, number of alignments to output is %d %s \n", outfilename,
			alignmentsToOutput->size(), (_skip_algnmnt_krnl? "(_skip_alignment is set)":"") );

	for(auto itr = alignmentsToOutput->begin(); itr != alignmentsToOutput->end(); itr++) {
		myoutputfile << (*readNameMap)[ (*itr).r1_id ] << " ";
		myoutputfile << (*readNameMap)[ (*itr).r2_id ] << " ";
		myoutputfile << (*itr) << std::endl;
	}
	myoutputfile.close();

	return 0;
}

int64_t Aligner::outputAlignmentsToFileWithMpi(std::vector<alignment_data>* alignmentsToOutput, std::string outfilename, int64_t file_offset)
{
	// initialization
	int buffer_size = MAX_READ_NAME_LEN*1024*1024; // arbitrary
	Buffer outBuffer = initBuffer(buffer_size);
	char* outChars;
	char* bufferIndex;
	std::stringstream ss;
	std::string tmp_string;
	size_t tmp_length_bytes;
	int64_t my_num_bytes;
	int64_t myoffset;
	int myrank, nprocs;
	CHECK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &nprocs) );
	CHECK_MPI( MPI_Comm_rank(MPI_COMM_WORLD, &myrank) );
	MPI_File outfile;
	MPI_Status status;
	int byte_count;
	bool more_to_write = (alignmentsToOutput->size() > 0);
	bool glob_more_to_write;
	auto itr = alignmentsToOutput->begin();
	int64_t end_offset = 0;

	CHECK_MPI( MPI_Allreduce(&more_to_write, &glob_more_to_write, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD));
	if (!glob_more_to_write) return file_offset;
	// TODO update file deletion when checkpoint-restart code in place
	//if(myrank == 0) MPI_File_delete((char*) outfilename.c_str(), MPI_INFO_NULL); // delete first (if it exists do not check error status)
	//CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );

	// open the output file collectively
	CHECK_MPI( MPI_File_open(MPI_COMM_WORLD, (char*) outfilename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile) );
	//begin write loop
	while (glob_more_to_write) {
		my_num_bytes = 0;
		myoffset = 0;
		more_to_write = 0;
		resetBuffer(outBuffer);
		// growBuffer?
		outChars = getStartBuffer(outBuffer);
		for(; itr != alignmentsToOutput->end(); ) { // don't advance itr until the alignment is copied to the write buffer
				ss << (*readNameMap)[ (*itr).r1_id ] << " ";
				ss << (*readNameMap)[ (*itr).r2_id ] << " ";
				ss << (*itr) << std::endl;
				tmp_string = ss.str();
				tmp_length_bytes = tmp_string.length();
				if (buffer_size >= my_num_bytes + tmp_length_bytes) {
					//strcat(outChars + my_num_bytes, tmp_string.c_str());
					memcpy(outChars + my_num_bytes, (const char*) tmp_string.c_str(), tmp_length_bytes);
					my_num_bytes += tmp_length_bytes;
					itr++;
					ss.str(""); // reset the stringstream, clear() just resets the status flags
					ss.clear();
				}
				else {
					more_to_write = 1;
					break;
				}
		}
		// MPI all write
		CHECK_MPI( MPI_Exscan(&my_num_bytes, &myoffset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD) );
		CHECK_MPI( MPI_File_write_at_all(outfile, file_offset+myoffset, outChars, my_num_bytes, MPI_CHAR, &status) );
		CHECK_MPI( MPI_Get_count(&status, MPI_CHAR, &byte_count) );
		assert( byte_count == my_num_bytes );
		CHECK_MPI( MPI_Allreduce(&more_to_write, &glob_more_to_write, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD));
	};
	//cleanup
	CHECK_MPI( MPI_File_close(&outfile) );
	freeBuffer(outBuffer);
	// calculate new offset
	file_offset += myoffset;
	CHECK_MPI(MPI_Allreduce(&file_offset, &end_offset, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD));
	return end_offset;
}

Aligner::Aligner(int kmerlength, int maxNumAlignPerPair, bool skipKernel, ReadId* readRanges, ReadOverlapMap* overlapMap)
: kmerLen(kmerlength), maxPerPairAlignments(maxNumAlignPerPair), skip_kernel(skipKernel), readRangeEnds(readRanges), readPairs(overlapMap)
{
       perProcRequestIds = new std::set< ReadId >[THREADS];
       globalReadIdOffset = (MYTHREAD==0? 0 : readRangeEnds[MYTHREAD-1]);
       readMap = new std::unordered_map<ReadId, std::string>();
       readNameMap = new std::unordered_map<ReadId, std::string>();
}

Aligner::Aligner(int kmerlength, int maxNumAlignPerPair, bool skipKernel, char *filename)
	: kmerLen(kmerlength), maxPerPairAlignments(maxNumAlignPerPair), skip_kernel(skipKernel)
{
	// delayed readRangeEnds initialization intentional
	globalReadIdOffset = 0;
	perProcRequestIds = new std::set< ReadId >[THREADS];
	readPairs = new ReadOverlapMap(); // TODO look at alternate constructors reserving space given readRanges up-front
	readMap = new std::unordered_map<ReadId, std::string>();
	readNameMap = new std::unordered_map<ReadId, std::string>();
	if (filename == NULL) {
		SWARN("Failed to open configuration file ", filename, " using default settings.");
	} else {
		ifstream infile(filename);
		if (infile.is_open()) {
			do {
				string line;
				// Skip comments
				while (!infile.eof() && infile.peek()=='#') {
					std::getline(infile, line);
				}
				if (infile.eof()) {
					SWARN("No uncommented lines in ", filename, " using default settings.");
					break;
				}
				for (auto field : {&myscoring.match, &myscoring.mismatch, &myscoring.gap, &myscoring.xdrop, &myscoring.static_threshold} ) {
					infile >> *field;
				}
				for (auto field : {&myscoring.adaptive, &myscoring.aligntoend}) {
					infile >> *field;
				}
				for (auto field : {&myscoring.error_rate, &myscoring.chernoff_delta, &myscoring.relaxation_margin}) {
					infile >> *field;
				}
				if (!MYTHREAD) {
					std::cout << "Successfully initialized pairwise alignment settings from input file "<< filename << std::endl;
				}
			} while (0);
			infile.close();
		}
	}
	seqan::Score<int, seqan::Simple>(myscoring.match, myscoring.mismatch, myscoring.gap);
	myscoring.set_adaptive_threshold();
	if (!MYTHREAD) { std::cout << "Pairwise alignment settings: " << myscoring << std::endl; }
}

Aligner::~Aligner() {
	if (readRangeEnds != NULL) { delete readRangeEnds; }
	if (perProcRequestIds != NULL) { delete[] perProcRequestIds; }
	if (readPairs != NULL) { delete readPairs; }
	if (readMap != NULL) { delete readMap; }
	if (readNameMap != NULL) { delete readNameMap; }
	//if (prevAlignmentsSet != NULL) { delete prevAlignmentsSet; }
}


} /* namespace align */
