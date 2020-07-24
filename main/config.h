/*
  The idea here is to get as many parameters as possible from the config file
 */

#ifndef __EXEC_CONFIG_H
#define __EXEC_CONFIG_H

#include <iostream>

#ifndef MAX_MER_SIZES
#define MAX_MER_SIZES 1 // dibella handles only single genome (vs. meta genome) right now
#endif

#ifndef MAX_LIBRARIES
#define MAX_LIBRARIES 128
#endif

#define CFG_ERR(fmt,...)                                        \
    do {                                                        \
        if (!MYTHREAD)                                          \
            fprintf(stderr, KLRED fmt KNORM, ##__VA_ARGS__);    \
        fflush(stderr);                                         \
        EXIT_FUNC(100);                                         \
    } while (0)

static const char *_usage = "-f config_file -N cores_per_node -B (cached_io) -S (save_intermediates) -I illumina_version "
    "-X sw_cache_mb -Y kmer_cache_mb -s stages -P firstTimeKmerDir -H (runHLL) -C (canonicalizeOutput) -h (help)";

typedef struct {
    char *files;
    char *name;
    int ins_avg;
    int std_dev;
    int read_len;
    int innie;
    int rev_comp;
    int for_contigging;
    int ono_setid;
    int for_gapclosing;
    int five_p_wiggle;
    int three_p_wiggle;
    int files_per_pair;
    int for_splinting;
    int num;
} lib_t;

typedef struct {
    char config_fname[MAX_FILE_PATH];
    int num_libs;
    lib_t libs[MAX_LIBRARIES];
    int is_metagenome;
    int is_diploid;
    int mer_size; // deprecated!
    int num_mer_sizes;
    int mer_sizes[MAX_MER_SIZES];
    int scaff_mer_size;
    int cores_per_node;
    int min_contig_len;
    int min_depth_cutoff;
    int extension_cutoff;
    double dynamic_min_depth;
    double alpha;
    double beta;
    double tau;
    double error_rate;
    int qual_offset;
    int bubble_min_depth_cutoff;
    int cached_io;
    int sw_cache_mb;
    int kmer_cache_mb;
    int nodes;
    int max_ono_setid;
    char ono_pair_thresholds[100];
    double gap_close_rpt_depth_ratio;
    int gap_close_aggressive;
    int colorize_output;
    int skip_la;
    double la_depth_dist_thres;
    double la_nvote_diff_thres;
    int la_mer_size_min;
    int la_mer_size_max;
    int la_hi_qual_thres;
    int la_qual_thres;
    double la_min_viable_depth;
    double la_min_expected_depth;
    char *stages;
    int assm_scaff_len_cutoff;
    int canonicalize;
    int num_scaff_loops;
    int illumina_version;
    int save_intermediates;
    char *pePath;
    int runHLL;
} cfg_t;

static char *get_option_list();
void get_list_lib_params(cfg_t *cfg, int *max_readlen, char *libs, char *list_insert_sizes, 
                         char *list_insert_sigmas, char *list_revcomps, char *list_5p_wiggles, 
                         char *list_3p_wiggles);
void get_all_lib_names(cfg_t *cfg, char *all_libs);
void get_config(int argc, char **argv, cfg_t *cfg);// { std::cout<< "hello from config.h" << std::endl; }
void print_lib(lib_t *lib);

#endif
