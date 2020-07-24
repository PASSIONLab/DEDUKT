/*
 * main.cpp
 *
 *      Author: mellis@lbl.gov
 */
#include <string.h>
#include <fstream>
#include <algorithm>
#include <climits>

#include "../version.h"
#include "../readoverlap/ReadOverlapper.h"
#include "../common/MPIUtils.h"
#include "../common/common.h"

#include <stdio.h>

#include "align.h"
#include "Aligner.h"

bool _skip_algnmnt_krnl = false; // internal option for development purposes

namespace align {
	const int MAX_PPALIGNMENTS_DFLT = USHRT_MAX; // reads are not longer than this
	const int MAX_PPALIGNMENTS_DFLT_MIN = 2; // empirically, most read pairs have 2 common k-mers
}

using namespace std;

//TODO add parameters for various alignment kernel options
int alignment_main(int argc, char** argv) {
	// performance statistics
	double start_time=0.0, io_time=0.0;

	const char *base_dir = ".";
	int kmer_len = min(17, MAX_KMER_SIZE-1); // minimum reasonable seed length
	char* inputFQ = NULL;
	char* inputOverlaps = NULL;
	char *bella_settings_file = NULL;
	bool cachedFastq = false;
	bool cachedOverlaps = false;
	bool loadperthread = false;
	string prevAlignments = ""; // TODO finish parameterizing
	string outputName = "";
	int maxPerPairAlignments = align::MAX_PPALIGNMENTS_DFLT;

	bool opt_err = false;
	option_t *optList, *thisOpt;
	optList = NULL;
	optList = GetOptList(argc, argv, "pBa:i:m:b:l:o:q:k:");
	while (optList != NULL) {
		thisOpt = optList;
		optList = optList->next;
		switch (thisOpt->option) {
			case 'i':
				inputFQ = thisOpt->argument;
				break;
			case 'm':
				inputOverlaps = thisOpt->argument;
				if (inputOverlaps == NULL) { SDIE("Please enter a valid output file name after -m. (none entered)"); }
				break;
			case 'b':
				bella_settings_file = thisOpt->argument;
				if (bella_settings_file == NULL) { SDIE("A valid file name is required after -b"); }
				break;
			case 'l':
				prevAlignments = ((string)thisOpt->argument);
				break;
			case 'o':
				outputName = thisOpt->argument;
				break;
			case 'q':
				maxPerPairAlignments = strtol(thisOpt->argument, NULL, 10);
				if (maxPerPairAlignments < 1) {
					serial_printf("Maximum alignments per overlapping read pair has been set to %d. Setting skip-alignment-kernel to TRUE.\n", maxPerPairAlignments);
					_skip_algnmnt_krnl = true;
				}
				break;
			case 'k':
				kmer_len = atoi(thisOpt->argument);
				if (kmer_len >= MAX_KMER_SIZE)
					SDIE("Please compile with a higher MAX_KMER_SIZE (%d) for kmer length %d\n", MAX_KMER_SIZE, kmer_len);
				if (GET_KMER_PACKED_LEN(kmer_len) > MAX_KMER_PACKED_LENGTH)
					SDIE("Please compile with a higher MAX_KMER_PACKED_LENGTH (%d) for kmer length %d\n", MAX_KMER_PACKED_LENGTH, kmer_len);
				break;
			case 'B':
				base_dir = "/dev/shm";
				cachedFastq = 1;
				break;
			case 'a':
				cachedOverlaps = true;
				serial_printf("%s:%s: set cached overlaps\n",__FILE__,__FUNCTION__);
                break;
			case 'p':
				loadperthread = true;
				break;
			/*
			case 'N':
				cores_per_node = atoi(thisOpt->argument);
				SET_CORES_PER_NODE(cores_per_node);
				break;
			*/
			default:
				opt_err = true;
				break;
		}
	}

	if (opt_err || !inputFQ || (inputOverlaps==NULL && !loadperthread))
	{
		if(MYTHREAD  == 0)
		{
			cout << "Usage: ./alignment -i input_fastq_list -m input_overlaps [ -k kmer_length -q max_per_pair -o output_file_name -b alignment_settings_file -B -a]" << endl;
			cout << "'input_fastq_list' is a fofn file, or a txt file containing the list of input fastq files (1 per line) " << endl;
			cout << "'input_overlaps' is output of the read-overlap (./overlap) stage " << endl;
			cout << "'alignment_settings_file' is the name of a file with settings for the pairwise aligner such as x-drop value, mismatch, match, and gap scoring, etc." << endl;
			cout << "'kmer_length' is the k-mer length (k) used as the minimum seed length for alignment, usually though not necessarily the same value used for the original kmer-analysis (kmermatch) stage" << endl;
			cout << "'max_per_pair is the maximum number of alignments to attempt per overlapping read pair and is optional. If the supplied value is <1, the alignment kernel will be skipped." << endl;
			cout << "Setting -B specifies that the input fastq files have been cached in /dev/shm/ " << endl;
			cout << "Setting -a specifies that the input overlap file(s) has been cached in /dev/shm/ " << endl;
		}
		return 0;
	}

	if (outputName == "") { outputName = "alignments-" + to_string(kmer_len); }

	if(MYTHREAD  == 0)
	{
		cout << "************* diBELLA read alignment step ************* " << endl;
		cout << "You are running with the following settings:" << endl;
		cout << "\tFastq caching is " << (cachedFastq? "" : "NOT ") << "enabled." << endl;
		cout << "\tOverlaps caching is " << (cachedOverlaps? "" : "NOT ") << "enabled." << endl;
		if (cachedFastq || cachedOverlaps) { cout << "\tCached I/O base directory is: " << base_dir << endl; }
		cout << "\tK-mer (minimum seed) length = " << kmer_len << endl;
		cout << "\tMax k-mer length (internal) = " << MAX_KMER_SIZE << endl;
		cout << "\tMax alignments to attempt per overlapping read pair = " << maxPerPairAlignments << endl;
		cout << "\tinput fastq file list: " << inputFQ << endl;
		cout << "\tinput overlap file: " << inputOverlaps << endl;
		cout << "\tinput pairwise alignment settings file: " << bella_settings_file << endl;
		cout << "\toutput file name: " << outputName << endl;
		cout << "\tskipping alignment kernel? " << (_skip_algnmnt_krnl? "True" : "False") << endl;
#ifdef BENCHMARKONLY
		cout << "\tBENCHMARKONLY set" << endl;
#else
		cout << "\tBENCHMARKONLY unset" << endl;
#endif
	}
	//inputOverlaps += "_" + to_string(MYTHREAD);

	/* for individual output files
	outputName += "_%d";
	char outputfilepath[MAX_FILE_PATH];
	snprintf(outputfilepath, MAX_FILE_PATH, outputName.c_str(), MYTHREAD);
	char *outputfile = get_rank_path(outputfilepath, MYTHREAD);
	*/

	//
	// initialize and run alignment
	//
	vector<align::alignment_data>* allAlignments = new vector<align::alignment_data>();
	align::Aligner* aligner = new align::Aligner(kmer_len, maxPerPairAlignments, _skip_algnmnt_krnl, bella_settings_file);
	aligner->computeAlignmentsFromOverlaps(inputOverlaps, loadperthread, prevAlignments, outputName.c_str(),
				inputFQ, cachedFastq, cachedOverlaps, base_dir, *allAlignments, true);
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

	// flush alignments to a file
	/* now does it in batches
	start_time = MPI_Wtime();
	aligner->outputAlignmentsToFile(allAlignments, outputName);
	CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD)); // timing
	io_time += MPI_Wtime() - start_time;
	serial_printf("Additional I/O time: %.2f s \n", io_time);
	*/

	//
	// clean-up
	//
	if (allAlignments) delete allAlignments;
	if (aligner) delete aligner;

	return 0;
}

#ifndef SINGLE_EXEC
StaticVars _sv = NULL;
int main(int argc, char **argv)
{
    CHECK_MPI( MPI_Init(&argc, &argv) );
    OPEN_MY_LOG("alignment");
    serial_printf("Starting diBELLA version %s on %d threads\n", DIBELLA_VERSION, THREADS);
    int ret = alignment_main(argc, argv);
    MPI_Finalize();
    return ret;
}
#endif


