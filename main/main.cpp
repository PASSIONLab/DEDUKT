/*
 * main.cpp
 *
 *  To be replaced with hipmer-like pipeline
 */
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <libgen.h>
#include <mpi.h>

#include "../version.h"
#include "../common/MPIUtils.h" // CHECK_MPI()
#include "../common/common.h"
//#include "../common/optlist.h"
// #include "../readoverlap/overlap.h" // overlap_main()
// #include "../align/align.h" // align_main()
#include "../kmercount/kmermatch.h"
#include "../kmercount/kmermatch.h" // kmermatch_main()
#include "../loadfq/loadfq.h" // loadfq_main()

#include "exec_utils.h"
#include "config.h"

//#define NO_GZIP 1 // unsupported option from hipmer
#define SINGLE_EXEC // TODO: this may need to be defined in the cmakelist instead
StaticVars _sv = NULL; // has to be initialized by the main executable - see note in StaticVars

#ifdef __linux__
  #define BASE_DIR(b) ((b) ? "/dev/shm" : ".")
#else
  #define BASE_DIR(b) "."
#endif
#define NOT_SCAFF 0
#define IS_SCAFF 1
#define MERGE_FLAG 0

static char _name_buff[MAX_FILE_PATH];

using namespace std;

static char *name_i(const char *name, int k)
{
    snprintf(_name_buff, MAX_FILE_PATH, "%s-%d", name, k);
    return _name_buff;
}

static void setup_fofns(const char *all_inputs_fofn, cfg_t *cfg)
{
    FILE *all_inputs_f = fopen_chk(all_inputs_fofn, "w");
    for (int i = 0; i < cfg->num_libs; i++) {
        lib_t *lib = &(cfg->libs[i]);
        char libname[MAX_FILE_PATH];
        snprintf(libname, MAX_FILE_PATH, "%s.fofn", lib->name);
        FILE *lib_f = fopen_chk(libname, "w");
        if (lib->files_per_pair == 2) {
            char *filename = strdup(lib->files);
            char *splitpt = strchr(filename, ',');
            splitpt[0] = 0;
            if (lib->for_contigging)
                fprintf(all_inputs_f, "%s\n%s\n", filename, splitpt + 1);
            fprintf(lib_f, "%s\n%s\n", filename, splitpt + 1);
            free(filename);
        } else {
            if (lib->for_contigging)
                fprintf(all_inputs_f, "%s\n", lib->files);
            fprintf(lib_f, "%s\n", lib->files);
        }
        fclose(lib_f);
    }
    fclose_track(all_inputs_f);
}

int main(int argc, char ** argv) {
	CHECK_MPI( MPI_Init(&argc, &argv) );
	int kmer_len;
	int reliable_max = MAX_NUM_READS;
	int err_threshold = 2;
	int seed_distance = 1000;
	int max_seeds = 2;
	bool cached_overlaps = false;
	bool use_perthread = false;
	bool purge_vm = 1;
	int type = 2;
	int window = 1;
	int mini_len = 3;

	//
	// setup memory usage tracking
	//
    double start_mem_free = 0;
    double start_mem_used = get_used_mem_gb();
    if (!MYTHREAD) {
    		start_mem_free = get_free_mem_gb();
    		cout << "RANK 0: start_mem_free = " << start_mem_free << endl;
    }

    // set up the per_thread directories
	OPEN_MY_LOG("diBELLA");
	init_diags();
	DBG("started DBG recording\n");
	LOGF("started LOGF recording\n");

	serial_printf("Starting diBELLA version %s on %d threads\n", DIBELLA_VERSION, THREADS);

	#ifdef DEBUG
		serial_printf("Debug build with log level %d\n", HIPMER_VERBOSITY);
	#endif

    char cwd[MAX_FILE_PATH];
    serial_printf(" Current working directory is %s \n", getcwd(cwd, MAX_FILE_PATH));
    //
    // set AUTORESTART_UFX option
    //
    // maybe useful in future, see hipmer/main/main.c ~526 for example usage
    /*
    char *autorestart_ufx_str = getenv("AUTORESTART_FIRST");
	int autorestart_ufx = 0;
	if (autorestart_ufx_str && atoi(autorestart_ufx_str) == 1) {
		serial_printf("# AUTORESTART_FIRST (restart after first stage to free memory \n");
		autorestart_ufx = 1;
	}
	*/

	//
	// setup the input files
	//
    char inufx[255];
	cfg_t cfg;
	cfg.num_mer_sizes = 1; // change from metahipmer
	cfg.cores_per_node = 1; //default
	cfg.cached_io=0;
	cfg.runHLL=0;
	cfg.save_intermediates = 0;
	//get_config(argc, argv, &cfg);
	option_t *this_opt;
	//char *myoptions = get_option_list();
	bool skip_algnmnt_krnl = false;
	char *bella_settings = NULL;
	option_t *opt_list = GetOptList(argc, argv, "i:e:h:DpBHESa:k:P:x:b:y:s:u:m:q:d:N:t:w:n:");
	char fail[10000];
	char buf[100];
	char *all_inputs_fofn = NULL;
	char *overlaps_fname = NULL;
	const char* default_stages = "all";
	cfg.stages=strdup(default_stages);
	while (opt_list) {
		this_opt = opt_list;
		opt_list = opt_list->next;
		switch (this_opt->option) {
			case 'i': {
				all_inputs_fofn = strdup(this_opt->argument);
				serial_printf("all_inputs_fofn set to %s\n", all_inputs_fofn);
				break;
			}
			case 'b': {
				bella_settings = strdup(this_opt->argument);
				break;
			}
			case 'N': cfg.cores_per_node = atoi(this_opt->argument); SET_CORES_PER_NODE(cfg.cores_per_node); break;
			case 'D': {
				purge_vm=0;
				serial_printf("won't purge VM\n");
				break;
			}
			case 'B': cfg.cached_io = 1; break; // not currently in optlist
			case 'a': {
				cached_overlaps = 1;
				serial_printf("%s:%s: cached overlaps set\n", __FILE__,__FUNCTION__);
				break;
			}
			case 'S': cfg.save_intermediates = 1; break; // not currently in optlist
			//case 'I': cfg.illumina_version = atoi(this_opt->argument); break;
			//case 'X': cfg.sw_cache_mb = atoi(this_opt->argument); break;
			//case 'Y': cfg.kmer_cache_mb = atoi(this_opt->argument); break;
			case 'k': // TODO if k is not set I think ufx will fallback on MAX_KMER_SIZE-1, but need to check/test...
			{
				cfg.mer_sizes[0] = strtol(this_opt->argument, NULL, 10);
				break;
			}
			case 'e': {
				err_threshold = strtol(this_opt->argument, NULL, 10);
				break;
			}
			case 't':
				type = strtol(this_opt->argument, NULL, 10);
				break;
			case 'w':
				window = strtol(this_opt->argument, NULL, 10);
				break;
			case 'n':
				mini_len = strtol(this_opt->argument, NULL, 10);
				break;

			case 's': cfg.stages = strdup(this_opt->argument); break;
			case 'P': cfg.pePath = strdup(this_opt->argument); break;
			case 'H': cfg.runHLL = 1; break;
			case 'h': CFG_ERR("Usage: %s %s\n", basename(argv[0]), _usage); break;
			case 'u':
				reliable_max = strtol(this_opt->argument, NULL, 10);
				if (reliable_max <= err_threshold) { SDIE("Upper limit set by -u must be greater than the error threshold set by -e.", reliable_max, err_threshold); }
				if (reliable_max < 2) { SDIE("Upper limit set by -u must be greater than 1, -u= %d", reliable_max); }
				if (reliable_max > MAX_NUM_READS) { SDIE("Upper limit set by -u is larger than MAX_NUM_READS, -u= %d. Use compilation with MAX_NUM_READS > %d", reliable_max, MAX_NUM_READS); }
				break;
			case 'm':
				overlaps_fname = this_opt->argument;
				break;
			case 'q': {
				max_seeds = strtol(this_opt->argument, NULL, 10);
				if (max_seeds < 1) { skip_algnmnt_krnl = true; }
				break;
			}
			case 'd': {
				seed_distance = strtol(this_opt->argument, NULL, 10);
				if (seed_distance < 1) {
					seed_distance=1000;
					SWARN("Minimum seed distance must be greater than or equal to 1 (-d %s). The seed distance has been reset to %d for this run.", this_opt->argument, seed_distance);
				}
				break;
			}
			case 'p':
				use_perthread = 1;
				break;
			default:
				sprintf(buf, "\tInvalid Option: %c\n", this_opt->option);
				strcat(fail, buf);
				break;
			}
	}
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
#ifndef __linux__
    if (cfg.cached_io) {
    	SWARN("Cached I/O is not yet supported on non-linux platforms. The pipeline will be run without this optimization.");
    	cfg.cached_io=0;
    }
#endif
	if (!all_inputs_fofn)
	{
		if(MYTHREAD  == 0)
		{
			cout << "Usage: ./main -i filename -e error_threshold -u reliable_upper_limit <-p prefix> <-m overlaps_file_name> <-s> <-H>" << endl;
			cout << "'filename' is a text file containing the paths to each file to be counted - one on each line -" << endl;
			cout << "'error_threshold' is the lowest number of occurrences of a k-mer for it to be considered a real (non-error) k-mer" << endl;
			cout << "reliable_upper_limit sets the maximum number of reads that will be recorded per k-mer" << endl;
			cout << "'prefix' is optional and if set, the code prints all k-mers with the prefix" << endl;
			cout << "'overlaps_file_name' (necessarily) specifies the input file name (if running the alignment stage alone), or (optionally) specifies a file name for the overlap output." << endl;
			cout << "-s is optional and if set, shared memory will be used to cache all files" << endl;
			cout << "-H is optional and if set, then HyperLogLog will be used to estimate the required bloomfilter" << endl;
			cout << "-N sets cores-per node (default 1)" << endl;
			cout << "-B is optional and if set, all fastq files will be loaded into virtual memory before k-mer analysis runs" << endl;
			cout << "-p overlaps will be output in per thread files between stages" << endl;
		}
		return 0;
	}

	if ( cfg.stages && strstr(cfg.stages, name_i("alignment", kmer_len)) == cfg.stages ) { // if the first stage specified is alignment, use overlaps file specified by -m
		if (overlaps_fname == NULL) {
			serial_printf("First stage specified is alignment, but no input overlaps specified. Exiting.\n");
			return 0;
		}
	}
	else if (overlaps_fname == NULL) {
		size_t name_length = strlen(OVERLAP_OUT_FNAME);
		overlaps_fname = (char*) malloc(name_length+1);
		sprintf(overlaps_fname, OVERLAP_OUT_FNAME);
		serial_printf("No output overlap file name provided. Overlaps will be written to %s.\n", overlaps_fname );
	}

	if (cfg.cached_io) {
		//char all_libs[4096] = "";
		//get_all_lib_names(&cfg, all_libs);
		serial_printf("all_inputs_fofn=%s\n", all_inputs_fofn);
		const char* stage = "loadfq";
		exec_stage( (char*) stage, NOT_SCAFF, loadfq_main, "loadfq",
				   "-N %d", cfg.cores_per_node,
				   "-f %s", all_inputs_fofn,
				   "-B %s", BASE_DIR(cfg.cached_io),
				   "-D %d", purge_vm,
					NULL);
	}

	//
	// generate k-mer -> read list, position list, #occurrences map (k-mer analysis)
	//
	//char inufx[255], all_inputs_fofn_kmer[255];
	kmer_len = cfg.mer_sizes[0]; // only use 1 k-mer size
	sprintf(inufx, "%s.ufx.bin", all_inputs_fofn);
	/*
	if (cfg.cached_io) {
		// if cached_io then don't have separate instances of ufx output for each kmer length since
		// we likely can't restart anyway
		strcpy(all_inputs_fofn_kmer, all_inputs_fofn);
		sprintf(inufx, "%s.ufx.bin", all_inputs_fofn);
	} else {
		// separate copies of ufx output to restart at any stage
		sprintf(all_inputs_fofn_kmer, "%s-%d", all_inputs_fofn, kmer_len);
		if (!MYTHREAD) {
			if (link(all_inputs_fofn, all_inputs_fofn_kmer) != 0) {
				if (errno != 17)
					DIE("Could not link %s to %s: %s\n", all_inputs_fofn, all_inputs_fofn_kmer,
						strerror(errno));
			}
		}
		sprintf(inufx, "%s-%d.ufx.bin", all_inputs_fofn, kmer_len);
		CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
	}
	*/
	//cfg.stages=NULL; // stages option not currently implemented

	//
	// computes read overlaps (read by read map) from k-mer map
	//
	int kmermatch_ran;
	/* -P not currently used in ufx
	char *pe_path_fmt = NULL;
	char *pe_path = NULL;
	if (cfg.pePath) {
		pe_path_fmt = strdup("-P %s");
		pe_path = strdup(cfg.pePath);
	}
	*/
	kmermatch_ran = exec_stage(cfg.stages, NOT_SCAFF, kmermatch_main, name_i("kmermatch", kmer_len),
						 "-k %d", kmer_len, "-t %d", type, "-w %d", window, "-n %d", mini_len, "-i %s", all_inputs_fofn, "-B %B", cfg.cached_io,
						 "-e %d", err_threshold, "-u %d", reliable_max, "-a %B", cached_overlaps,
						 "-P %B", use_perthread,
						 "-m %s", overlaps_fname, "-d %d", seed_distance, "-q %d", max_seeds,
						 "-H %B", cfg.runHLL, /*"-E %d", cfg.extension_cutoff,*/
						 "-S %B", cfg.save_intermediates,
						 /*pe_path_fmt, pe_path, NULL);*/
						 NULL);
	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
	/*
	if (cfg.pePath) {
		free(pe_path_fmt);
		free(pe_path);
	}
	*/
	serial_printf("%s %d : finished kmermatch stage with code %d\n", __FILE__, __LINE__, kmermatch_ran);

	// could use the histogram to calculate things at this point... see hipmer/main/main.c ~471

// #ifndef KMERA_ONLY
// 	//
// 	// compute local alignments for each read-read pair
// 	//
// 	int alignment_ran;
// 	char alignments_out[1024];
// 	sprintf(alignments_out, "alignments-%d", kmer_len);
// 	alignment_ran = exec_stage(cfg.stages, NOT_SCAFF, alignment_main, name_i("alignment", kmer_len),
// 						 "-k %d", kmer_len, "-q %d", max_seeds, "-b %s", bella_settings,
// 						 "-m %s", overlaps_fname, "-o %s", alignments_out, "-i %s", all_inputs_fofn,
// 						 "-B %B", cfg.cached_io, "-a %B", cached_overlaps, "-p %B", use_perthread, NULL);
// 	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
// #endif

	//
	// clean-up
	//
	free(all_inputs_fofn);
	CLOSE_MY_LOG;
	fini_diags();
	CHECK_MPI( MPI_Finalize() );
	return 0;
}

