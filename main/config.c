#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <libgen.h>

#include "../common/MPIUtils.h"
#include "../common/common.h"
#include "../common/optlist.h"
//#include "../common/timers.h"
#include "../common/StaticVars.h" //SET_CORES_PER_NODE
#include "config.h"
#include "exec_utils.h"

#define MAX_CFG_LINE_LEN 4096 

#define MAX_LINE_LEN 4096

//extern char *_use_col;

/*
shared int _qual_offset = 0;

static int get_qual_offset(int illumina_version, lib_t *lib)
{
    if (!MYTHREAD) {
        _qual_offset = 64;
        if (!illumina_version) {
            char *libname = lib->files;
            if (lib->files_per_pair == 2)
                libname = strchr(lib->files, ',') + 1;
                
            double startEstimate = now();
            serial_printf("%s# Estimating quality offset from %s as ", _use_col, libname);
            FILE *f = fopen_chk(libname, "r");
            char buf[MAX_LINE_LEN];
            int i = 0;
            int use_next_line = 0;
            while (fgets(buf, MAX_LINE_LEN, f)) {
                i++;
                if (i > 50000)
                    break;
                if (buf[0] == '+') {
                    use_next_line = 1;
                } else if (use_next_line) {
                    use_next_line = 0;
                    if (strpbrk(buf, "\"#$%&'()*,./0123456789:;<=>?")) {
                        _qual_offset = 33;
                        break;
                    }
                }
            }
            fclose_track(f);
            serial_printf("%d in %0.3f s" KNORM "\n", _qual_offset, now() - startEstimate);
        } else if (illumina_version == 13 || illumina_version == 15) {
            _qual_offset = 64;
        } else if (illumina_version == 18) {
            _qual_offset = 33;
        } else {
            _qual_offset = 0;
        }
    }
    UPC_LOGGED_BARRIER;
    return _qual_offset;
}
*/

/*
static int get_kmer_cache_mb(int cores_per_node)
{
    int shared_heap_mb = get_shared_heap_mb();
    if (!shared_heap_mb) {
        serial_printf("%s# Cannot determine UPC shared heap size and hence cannot determine KMER cache size" KNORM "\n",
                      _use_col);
        return 0;
    }
    // set kmer cache as up to 20% of the shared memory of a node
    int cache_size = shared_heap_mb * cores_per_node * 0.2;
    int min_cache = 128;
    if (cache_size < min_cache)
        cache_size = min_cache;
    return cache_size;
}
*/

/*
static int get_sw_cache_mb(int cores_per_node, int kmer_cache_mb)
{
    int shared_heap_mb = get_shared_heap_mb();
    if (!shared_heap_mb) {
        serial_printf("%s# Cannot determine UPC shared heap size and hence cannot determine SW cache size" KNORM "\n", 
                      _use_col);
        return 0;
    }
    // set the SW cache to up to 40% of the remaining shared memory
    int cache_size = (shared_heap_mb * cores_per_node - kmer_cache_mb) * 0.4;
    // make sure it's no more than 40% of the kmer cache
    int max_cache = kmer_cache_mb * 0.4;
    if (cache_size > max_cache)
        cache_size = max_cache;
    return cache_size;
}
*/

static int get_bool(const char *val)
{
    if (strlen(val) >= 1) {
        if (val[0] == '0') 
            return 0;
        if (val[0] == '1')
            return 1;
    }
    return -1;
}

static char* get_next_lib_tok(char **saveptr, int line)
{
    char *tok = strtok_r(NULL, " \t", saveptr);
    if (!tok)
        CFG_ERR("Invalid or missing token at line %d\n", line);
    return tok;
}

static int get_bool_lib(const char *key, char **saveptr, int line)
{
    char *val = get_next_lib_tok(saveptr, line);
    int b = get_bool(val);
    if (b == -1)
        CFG_ERR("Error on line %d: %s must be either 1 or 0, not '%s'\n", line, key, val);
    return b;
}

static int get_str_param(const char *key, char *token, char *val, char *s) 
{
    if (strcmp(key, token) == 0) {
        strcpy(s, val);
        if (s[strlen(s) - 1] == '\n')
            s[strlen(s) - 1] = 0;
        return 1;
    }
    return 0;
}

static int get_int_param(const char *key, char *token, char *val, int *i) 
{
    if (strcmp(key, token) == 0) {
        *i = atoi(val);
        return 1;
    }
    return 0;
}

static int get_double_param(const char *key, char *token, char *val, double *d)
{
    if (strcmp(key, token) == 0) {
        *d = atof(val);
        return 1;
    }
    return 0;
}

static int get_bool_param(const char *key, char *token, char *val, int *b)
{
    if (strcmp(key, token) == 0) {
        *b = get_bool(val);
        if (*b == -1)
            CFG_ERR("Error: %s must be either 1 or 0, not '%s'\n", key, val);
        return 1;
    }
    return 0;
}

static int get_int_array_param(const char *key, char *token, char *val, int *arr, int *num) 
{
    if (strcmp(key, token) != 0)
        return 0;

    char *saveptr;
    char *tok = strtok_r(val, ",", &saveptr);
    *num = 0;
    while (tok) {
        arr[*num] = atoi(tok);
        (*num)++;
        tok = strtok_r(NULL, ",", &saveptr);
    }
    if (!(*num))
        return 0;
    return 1;
}

static int skip_line(const char *line)
{
    while (*line) {
        if (!isspace(*line)) {
            if (*line == '#')
                return 1;
            return 0;
        }
        line++;
    }
    return 1;        
}

static void read_config_file(cfg_t *cfg)
{
    FILE *f = fopen_chk(cfg->config_fname, "r");
    char buf[MAX_CFG_LINE_LEN];
    char *saveptr;
    int line = -1;
    cfg->max_ono_setid = 0;
    while (fgets(buf, MAX_CFG_LINE_LEN, f)) {
        line++;
        if (skip_line(buf))
            continue;
        char *key = strtok_r(buf, " \t", &saveptr);
        if (!key) 
            CFG_ERR("Error parsing config file at line %d: expected a key string\n", line);
        if (strcmp(key, "lib_seq") == 0) {
            lib_t *lib = &cfg->libs[cfg->num_libs];
            lib->files = strdup(get_next_lib_tok(&saveptr, line));
            lib->name = strdup(get_next_lib_tok(&saveptr, line));
            lib->ins_avg = atoi(get_next_lib_tok(&saveptr, line));
            lib->std_dev = atoi(get_next_lib_tok(&saveptr, line));
            lib->read_len = atoi(get_next_lib_tok(&saveptr, line));
            lib->innie = get_bool_lib("innie", &saveptr, line);
            lib->rev_comp = get_bool_lib("revcomp", &saveptr, line);
            lib->for_contigging = get_bool_lib("contigging", &saveptr, line);
            lib->ono_setid = atoi(get_next_lib_tok(&saveptr, line));
            lib->for_gapclosing = get_bool_lib("gapclosing", &saveptr, line);
            lib->five_p_wiggle = atoi(get_next_lib_tok(&saveptr, line));
            lib->three_p_wiggle = atoi(get_next_lib_tok(&saveptr, line));
            lib->files_per_pair = atoi(get_next_lib_tok(&saveptr, line));
            lib->for_splinting =  get_bool_lib("splinting", &saveptr, line);

            if (lib->ono_setid > cfg->max_ono_setid)
                cfg->max_ono_setid = lib->ono_setid;

            if (strchr(lib->files, ',')) {
                if (lib->files_per_pair != 2) 
                    CFG_ERR("Invalid config file at line %d: %s has a comma in the filename field but %d"
                            " in FilesPerPair\n",
                            line, lib->files, lib->files_per_pair);
                if (!lib->read_len) 
                    CFG_ERR("Invalid fastqs for %s on line %d.  Currently variable length paired read file in"
                            " two separate files is NOT supported.  Please convert them to interleaved files"
                            " prior to executing HipMer.  You may want to use the "
                            "${HIPMER_INSTALL}/bin/interleave_fastq executable\n", 
                            lib->files, line);
            }
            if (!lib->five_p_wiggle) {
                lib->five_p_wiggle = 5;
                serial_printf("# In config line %d: Setting 5' wiggle room in %s to the default of 5" KNORM "\n",
                              line, lib->files);
            }
            if (!lib->three_p_wiggle) {
                lib->three_p_wiggle = 5;
                serial_printf("# In config line %d: Setting 3' wiggle room in %s to the default of 5" KNORM "\n",
                              line, lib->files);
            }
            lib->num = cfg->num_libs++;
        } else {
            char *val = strtok_r(NULL, " \t", &saveptr);
            if (!val)
                CFG_ERR("Error parsing config file at line %d: expected a value\n", line);
            if (get_int_param("min_contig_len", key, val, &cfg->min_contig_len)) continue;
            if (get_int_param("min_depth_cutoff", key, val, &cfg->min_depth_cutoff)) continue;
            if (get_int_param("extension_cutoff", key, val, &cfg->extension_cutoff)) continue;
            if (get_int_param("mer_size", key, val, &cfg->mer_size)) continue;
            if (get_int_array_param("mer_sizes", key, val, cfg->mer_sizes, &cfg->num_mer_sizes)) continue;
            if (get_int_param("scaff_mer_size", key, val, &cfg->scaff_mer_size)) continue;
            if (get_double_param("dynamic_min_depth", key, val, &cfg->dynamic_min_depth)) continue;
            if (get_double_param("alpha", key, val, &cfg->alpha)) continue;
            if (get_double_param("beta", key, val, &cfg->beta)) continue;
            if (get_double_param("tau", key, val, &cfg->tau)) continue;
            if (get_double_param("error_rate", key, val, &cfg->error_rate)) continue;
            if (get_int_param("bubble_min_depth_cutoff", key, val, &cfg->bubble_min_depth_cutoff)) continue;
            if (get_bool_param("is_metagenome", key, val, &cfg->is_metagenome)) continue;
            if (get_bool_param("is_diploid", key, val, &cfg->is_diploid)) continue;
            if (get_str_param("ono_pair_thresholds", key, val, cfg->ono_pair_thresholds)) continue;
            if (get_double_param("gap_close_rpt_depth_ratio", key, val, &cfg->gap_close_rpt_depth_ratio)) continue;
            if (get_bool_param("gap_close_aggressive", key, val, &cfg->gap_close_aggressive)) continue;
            if (get_int_param("assm_scaff_len_cutoff", key, val, &cfg->assm_scaff_len_cutoff)) continue;
            if (get_int_param("num_scaff_loops", key, val, &cfg->num_scaff_loops)) continue;
            // all the local assembly params
            if (get_int_param("skip_la", key, val, &cfg->skip_la)) continue;
            if (get_double_param("la_depth_dist_thres", key, val, &cfg->la_depth_dist_thres)) continue;
            if (get_double_param("la_nvote_diff_thres", key, val, &cfg->la_nvote_diff_thres)) continue;
            if (get_int_param("la_mer_size_min", key, val, &cfg->la_mer_size_min)) continue;
            if (get_int_param("la_mer_size_max", key, val, &cfg->la_mer_size_max)) continue;
            if (get_int_param("la_qual_thres", key, val, &cfg->la_qual_thres)) continue;
            if (get_int_param("la_hi_qual_thres", key, val, &cfg->la_hi_qual_thres)) continue;
            if (get_double_param("la_min_viable_depth", key, val, &cfg->la_min_viable_depth)) continue;
            if (get_double_param("la_min_expected_depth", key, val, &cfg->la_min_expected_depth)) continue;
            if (get_int_param("illumina_version", key, val, &cfg->illumina_version)) continue;
            if (get_bool_param("save_intermediates", key, val, &cfg->save_intermediates)) continue;
            if (!MYTHREAD)
                WARN("Ignoring unrecognized key in config file at line %d: %s\n", line + 1, key);
        }
    }
    fclose_track(f);
}

static char *get_option_list()
{
    char *usage = strdup(_usage);
    // get option characters from usage
    char *options = (char*) malloc_chk(1000);
    int pos = 0;
    char *aux = NULL;
    char *token = strtok_r(usage, " ", &aux);
    while (token) {
        if (token[0] == '-') 
            options[pos++] = token[1];
        else if (token[0] != '(') 
            options[pos++] = ':';
        token = strtok_r(NULL, " ", &aux);
    }
    options[pos] = 0;
    return options;
}

static void clean_end(char *list)
{
    int len = strlen(list);
    if (list[len - 1] == ',')
        list[len - 1] = 0;
}

void add_list_elem(char *list, char *elem, int i,  int n)
{
    strcat(list, elem);
    if (i < n - 1)
        strcat(list, ",");
}

void get_list_lib_params(cfg_t *cfg, int *max_readlen, char *libs, char *list_insert_sizes, 
                         char *list_insert_sigmas, char *list_revcomps, char *list_5p_wiggle, 
                         char *list_3p_wiggle)
{
    char buf[256];
    *max_readlen = 0;
    for (int i = 0; i < cfg->num_libs; i++) {
        lib_t *lib = &cfg->libs[i];
        if (!lib->for_gapclosing)
            continue;
        if (lib->read_len > *max_readlen)
            *max_readlen = lib->read_len;
        add_list_elem(libs, lib->name, i, cfg->num_libs);
        
        sprintf(buf, "%d", lib->ins_avg);
        add_list_elem(list_insert_sizes, buf, i, cfg->num_libs);
        sprintf(buf, "%d", lib->std_dev);
        add_list_elem(list_insert_sigmas, buf, i, cfg->num_libs);
        
        sprintf(buf, "%d", lib->rev_comp);
        add_list_elem(list_revcomps, buf, i, cfg->num_libs);
        sprintf(buf, "%d", lib->five_p_wiggle);
        add_list_elem(list_5p_wiggle, buf, i, cfg->num_libs);
        sprintf(buf, "%d", lib->three_p_wiggle);
        add_list_elem(list_3p_wiggle, buf, i, cfg->num_libs);
    }
    clean_end(libs);
    clean_end(list_insert_sizes);
    clean_end(list_insert_sigmas);
    clean_end(list_revcomps);
    clean_end(list_5p_wiggle);
    clean_end(list_3p_wiggle);
}

void get_all_lib_names(cfg_t *cfg, char *all_libs)
{
    char buf[256];
    for (int i = 0; i < cfg->num_libs; i++) {
        lib_t *lib = &cfg->libs[i];
        add_list_elem(all_libs, lib->name, i, cfg->num_libs);
    }
    clean_end(all_libs);
}

void get_config(int argc, char **argv, cfg_t *cfg) 
{
    if (argc == 1) {
        CFG_ERR("Usage: %s %s\n", basename(argv[0]), _usage);
    }
    // defaults are all 0
    memset(cfg, 0, sizeof(cfg_t));

    // by default only loop once through scaffolding
    cfg->num_scaff_loops = 1;

    // defaults for progressiveRelativeDepth
    cfg->alpha = 0.1;
    cfg->beta = 0.2;
    cfg->tau = 2.0;

    // except defaults for localassm
    cfg->skip_la = 0;
    cfg->la_depth_dist_thres = 0.5;
    cfg->la_nvote_diff_thres = 1.5;
    cfg->la_mer_size_min = 13;
    cfg->la_mer_size_max = -1;
    cfg->la_qual_thres = 10;
    cfg->la_hi_qual_thres = 20;
    cfg->la_min_viable_depth = 0.2;
    cfg->la_min_expected_depth = 0.5;
    cfg->pePath = NULL;
    cfg->kmer_cache_mb = -1;
    cfg->sw_cache_mb = -1;
    cfg->assm_scaff_len_cutoff = 200;
    cfg->extension_cutoff = 20;

    // get command line params
    option_t *this_opt;
    char *myoptions = get_option_list();
    option_t *opt_list = GetOptList(argc, argv, myoptions);
    char fail[10000] = "";
    char buf[100];
    while (opt_list) {
        this_opt = opt_list;
        opt_list = opt_list->next;
        switch (this_opt->option) {
        case 'f': strcpy(cfg->config_fname, this_opt->argument); break; // TODO change to just feed ufx the fastq files in the normal way
        case 'N': cfg->cores_per_node = atoi(this_opt->argument); SET_CORES_PER_NODE(cfg->cores_per_node); break;
        case 'B': cfg->cached_io = 1; break;
        case 'S': cfg->save_intermediates = 1; break;
        case 'Q': {
        		cfg->qual_offset = atoi(this_opt->argument);
        		printf("cfg->qual_offset set to %d \n", cfg->qual_offset);
        		break;
        }
        //case 'I': cfg->illumina_version = atoi(this_opt->argument); break;
        //case 'X': cfg->sw_cache_mb = atoi(this_opt->argument); break;
        //case 'Y': cfg->kmer_cache_mb = atoi(this_opt->argument); break;
        case 'k': // if k is not set I think ufx will fallback on MAX_KMER_SIZE-1
        {
        		cfg->num_mer_sizes = 1;
        	    cfg->mer_sizes[0] = cfg->mer_size;
        	    break;
        }
        case 'e': cfg->min_depth_cutoff = atoi(this_opt->argument); break;
        case 's': cfg->stages = strdup(this_opt->argument); break;
        case 'P': cfg->pePath = strdup(this_opt->argument); break;
        case 'H': cfg->runHLL = 1; break;
        case 'C': cfg->canonicalize = 1; break;
        case 'h': CFG_ERR("Usage: %s %s\n", basename(argv[0]), _usage);
        default:
            sprintf(buf, "\tInvalid Option: %c\n", this_opt->option);
            strcat(fail, buf);
            break;
        }
    }

    // check command line params
    if (!cfg->config_fname[0])
        strcat(fail, "\tconfig_file required (-f)\n");
    if (!cfg->cores_per_node) {
        // estimate cores per node from procfs
        cfg->cores_per_node = get_cores_per_node();
        if (cfg->cores_per_node <= 0)
            strcat(fail, "Cannot determine cores per node from /proc/cpuinfo. Please specify with -N\n");
        else 
            serial_printf("# Determined cores per node as %d" KNORM "\n", cfg->cores_per_node);
    }
    if (cfg->cores_per_node > THREADS) {
        serial_printf("# Cores per node (%d) is greater than threads (%d), restricting to thread count" KNORM "\n",
                      cfg->cores_per_node, THREADS);
        cfg->cores_per_node = THREADS;
    }
    /*//TODO
    if (cfg->kmer_cache_mb < 0) {
        cfg->kmer_cache_mb = get_kmer_cache_mb(cfg->cores_per_node);
        if (cfg->kmer_cache_mb <= 0)
            strcat(fail, "Cannot determine KMER cache size. Please specify with -Y\n");
        else
            serial_printf("%s# Calculated KMER cache size as %d MB" KNORM "\n", _use_col, cfg->kmer_cache_mb);
    }
    if (cfg->sw_cache_mb < 0) {
        cfg->sw_cache_mb = get_sw_cache_mb(cfg->cores_per_node, cfg->kmer_cache_mb);
        if (cfg->sw_cache_mb <= 0)
            strcat(fail, "Cannot determine software cache size. Please specify with -X\n");
        else
            serial_printf("%s# Calculated software cache size as %d MB" KNORM "\n", _use_col, cfg->sw_cache_mb);
    }
    */
    if (fail[0]) {
        CFG_ERR("Error in command line parameters:\n%s\nUsage: %s %s\n", 
                fail, basename(argv[0]), _usage);
    }

    // some defaults
    /* artifacts from hipmer
    strcpy(cfg->ono_pair_thresholds, "1,10");
    cfg->gap_close_rpt_depth_ratio = 0;
    cfg->gap_close_aggressive = 1;
    */

    // now read in config file
    /*
     * code from hipmer. dibella doesn't yet use standard format config files. we just pass the fastq file list to ufx.
     */
    /*
    read_config_file(cfg);

    // check config file params
    if (!cfg->num_mer_sizes) {
        if (!cfg->mer_size) { // legacy parameter
            strcat(fail, "\tneed to set at least one entry in mer_sizes\n");
        } else {
            SWARN("the config option mer_size is deprecated, please use the mer_sizes parameter in the future\n");
            cfg->num_mer_sizes = 1;
            cfg->mer_sizes[0] = cfg->mer_size;
            cfg->mer_size = 0; // to detect if they are both set!
        }
    }
    if (cfg->num_mer_sizes > MAX_MER_SIZES) {
        strcat(fail, "\ttoo many mer sizes were specified, please reduce or recompile with a larger MAX_MER_SIZES configuration parameter\n");
    }
    if (cfg->mer_size) {
        strcat(fail, "\tthe config option mer_size is deprecated, and mer_sizes is also defined.  Please just specify mer_sizes\n");
    }
    for (int i = 0; i < cfg->num_mer_sizes - 1; i++) {
        if (cfg->mer_sizes[i] >= cfg->mer_sizes[i + 1]) {
            strcat(fail, "\tmer sizes need to be increasing\n");
            break;
        }
    }
    if (!cfg->num_libs)
        strcat(fail, "\tat least one lib_seq required\n");
    if (!cfg->min_depth_cutoff) {
        SWARN("min_depth_cutoff is not specified, choosing 2.  Please consider setting this important parameter considering your expected coverage of your data set.\n");
        cfg->min_depth_cutoff = 2;
    }
    if (!cfg->dynamic_min_depth) {
        SWARN("dynamic_min_depth is not specified, choosing 1.0 (i.e. disabled). Please consider setting this important parameter considering your data set kmer-error rates.\n");
        cfg->dynamic_min_depth = 1.0;
    }

    if (cfg->is_metagenome) {
        if (cfg->max_ono_setid > 1 && cfg->num_scaff_loops > 1)
            strcat(fail, "\tcannot have multiple oNo rounds for metagenomes\n");
        if (!cfg->alpha)
            strcat(fail, "\talpha is required with is_metagenome\n");
        if (!cfg->beta)
            strcat(fail, "\tbeta is required with is_metagenome\n");
        if (!cfg->tau)
            strcat(fail, "\ttau is required with is_metagenome\n");
        if (!cfg->error_rate) {
            if (cfg->dynamic_min_depth < 1.0) {
                cfg->error_rate = cfg->dynamic_min_depth;
            } else {
                cfg->error_rate = 0.9;
                SWARN("Setting error_rate to %f, Please consider setting this important metagenome parameter and/or dynamic_min_depth.\n", cfg->error_rate); 
            }
        }
        if (!cfg->skip_la) {
            if (cfg->la_mer_size_max <= 0) {
                // using the default, so set it to 2+ the max
                cfg->la_mer_size_max = 2 + cfg->mer_sizes[cfg->num_mer_sizes - 1];
            }
            if (cfg->la_mer_size_min >= cfg->la_mer_size_max)
                strcat(fail, "\tla_mer_size_max must be larger than la_mer_size_min\n");
            if (cfg->num_mer_sizes && cfg->mer_sizes[cfg->num_mer_sizes - 1] > cfg->la_mer_size_max) {
                SWARN("largest mer_size must be larger than la_mer_size_max (it was %d, the max is %d)\n", cfg->la_mer_size_max, cfg->mer_sizes[cfg->num_mer_sizes-1]);
                cfg->la_mer_size_max = 2 + cfg->mer_sizes[cfg->num_mer_sizes - 1];
            }
        }
    }
    if (fail[0]) {
        CFG_ERR("Error in config file:\n%s\nUsage: %s %s\n", fail, basename(argv[0]), _usage);
    }
    if (!cfg->scaff_mer_size)
        cfg->scaff_mer_size = cfg->mer_sizes[cfg->num_mer_sizes - 1];
    if (!cfg->min_contig_len)
        cfg->min_contig_len = cfg->mer_sizes[0];
    if (!cfg->is_diploid && !cfg->is_metagenome && cfg->min_contig_len < cfg->mer_sizes[0] * 2) {
        serial_printf("# Setting min_contig_len to %d because this is not a diploid or metagenome" KNORM "\n",
                      cfg->mer_sizes[0] * 2);
        cfg->min_contig_len = cfg->mer_sizes[0] * 2;
    }
    */

    cfg->nodes = THREADS / cfg->cores_per_node;

    /*// TODO
    cfg->qual_offset = get_qual_offset(cfg->illumina_version, &cfg->libs[0]);
    if (cfg->qual_offset <= 0) 
        CFG_ERR("Error in command line parameters:\n%s\nUsage: %s %s\n", 
                "Could not determine quality offset. Please specify with -I\n", basename(argv[0]), _usage);
	*/

    // print params
    if (!MYTHREAD) {
        serial_printf("\n%s" KNORM "\n", HASH_BAR);
        serial_printf("# Parameters for %s:" KNORM "\n", argv[0]);
        serial_printf("#" KNORM "\n");
        serial_printf("#  config_file:               %s" KNORM "\n", cfg->config_fname);
        serial_printf("#  Libs:                      ");
        for (int i = 0; i < cfg->num_libs; i++)
            serial_printf("%s ", cfg->libs[i].files);
        serial_printf("" KNORM "\n");
        serial_printf("#  cores_per_node:            %d" KNORM "\n", cfg->cores_per_node);
        	serial_printf("#  nodes:                     %d" KNORM "\n", cfg->nodes);
        serial_printf("#  min_depth_cutoff:          %d" KNORM "\n", cfg->min_depth_cutoff);
        serial_printf("#  extension_cutoff:          %d" KNORM "\n", cfg->extension_cutoff);
        serial_printf("#  dynamic_min_depth:         %.3f" KNORM "\n", cfg->dynamic_min_depth);
        serial_printf("#  illumina_version::         %d" KNORM "\n", cfg->illumina_version);
        serial_printf("#  qual_offset:               %d" KNORM "\n", cfg->qual_offset);
        serial_printf("#  mer_sizes:                 ");
        for (int i = 0; i < cfg->num_mer_sizes; i++)
            serial_printf("%d ", cfg->mer_sizes[i]);
        serial_printf("" KNORM "\n");
        if (cfg->is_metagenome) {
            serial_printf("#  alpha:                     %.3f" KNORM "\n", cfg->alpha);
            serial_printf("#  beta:                      %.3f" KNORM "\n", cfg->beta);
            serial_printf("#  tau:                       %.3f" KNORM "\n", cfg->tau);
            serial_printf("#  error_rate:                %.3f" KNORM "\n", cfg->error_rate);
            serial_printf("#  skip_la:                   %d" KNORM "\n", cfg->skip_la);
            serial_printf("#  la_depth_dist_thres:       %.2f" KNORM "\n", cfg->la_depth_dist_thres);
            serial_printf("#  la_nvote_diff_thres:       %.2f" KNORM "\n", cfg->la_nvote_diff_thres);
            serial_printf("#  la_mer_size_min:           %d" KNORM "\n", cfg->la_mer_size_min);
            serial_printf("#  la_mer_size_max:           %d" KNORM "\n", cfg->la_mer_size_max);
            serial_printf("#  la_qual_thres:             %d" KNORM "\n", cfg->la_qual_thres);
            serial_printf("#  la_hi_qual_thres:          %d" KNORM "\n", cfg->la_hi_qual_thres);
            serial_printf("#  la_min_viable_depth:       %.2f" KNORM "\n", cfg->la_min_viable_depth);
            serial_printf("#  la_min_expected_depth:     %.2f" KNORM "\n", cfg->la_min_expected_depth);
        } else {
            serial_printf("#  min_contig_len:            %d" KNORM "\n", cfg->min_contig_len);
        }
        serial_printf("#  scaff_mer_size:            %d" KNORM "\n", cfg->scaff_mer_size);

        if (cfg->is_metagenome || cfg->is_diploid)
            serial_printf("#  bubble_min_depth_cutoff:   %d" KNORM "\n", cfg->bubble_min_depth_cutoff);
        serial_printf("#  is_diploid:                %s" KNORM "\n", cfg->is_diploid ? "true" : "false");
        serial_printf("#  is_metagenome:             %s" KNORM "\n", cfg->is_metagenome ? "true" : "false");
        serial_printf("#  kmer_cache_mb:             %d" KNORM "\n", cfg->kmer_cache_mb);
        serial_printf("#  sw_cache_mb:               %d" KNORM "\n", cfg->sw_cache_mb);
        serial_printf("#  ono_pair_thresholds:       %s" KNORM "\n", cfg->ono_pair_thresholds);
        serial_printf("#  gap_close_rpt_depth_ratio: %.3f" KNORM "\n", cfg->gap_close_rpt_depth_ratio);
        serial_printf("#  gap_close_aggressive:      %s" KNORM "\n", cfg->gap_close_aggressive ? "true" : "false");
        serial_printf("#  use cached IO:             %s" KNORM "\n", cfg->cached_io ? "true" : "false");
        serial_printf("#  save_intermediates:        %s" KNORM "\n", cfg->save_intermediates ? "true" : "false");
        serial_printf("#  assm_scaff_len_cutoff:     %d" KNORM "\n", cfg->assm_scaff_len_cutoff);
        serial_printf("#  num_scaff_loops:           %d" KNORM "\n", cfg->num_scaff_loops);
        serial_printf("#  canonicalize_output:       %s" KNORM "\n", cfg->canonicalize ? "true" : "false");
        serial_printf("#  run_HLL:                   %s" KNORM "\n", cfg->runHLL ? "true" : "false");
        if (cfg->pePath)
            serial_printf("#  first_time_kmer_dir:       %s" KNORM "\n", cfg->pePath);
        if (cfg->stages)
            serial_printf("#  running stages:            %s" KNORM "\n", cfg->stages);
        serial_printf("%s" KNORM "\n", HASH_BAR);
    }
    serial_printf("hello!\n");
    //free_chk(myoptions); // TODO this segfaults on the first call from main in the refactored version, but it's not clear why
    serial_printf("hello?\n");
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    serial_printf("hello??\n");
}

void print_lib(lib_t *lib) 
{
    printf(KLGREEN);
    printf("LLLLib %s\n", lib->name);
    printf("LLL\tfiles          %s\n", lib->files);
    printf("LLL\tins_avg        %d\n", lib->ins_avg);
    printf("LLL\tstd_dev        %d\n", lib->std_dev);
    printf("LLL\tread_len       %d\n", lib->read_len);
    printf("LLL\tinnie          %d\n", lib->innie);
    printf("LLL\trev_comp       %d\n", lib->rev_comp);
    printf("LLL\tfive_p_wiggle  %d\n", lib->five_p_wiggle);
    printf("LLL\tthree_p_wiggle %d\n", lib->three_p_wiggle);
    printf("LLL\tfiles_per_pair %d\n", lib->files_per_pair);
    printf("LLL\tfor_contigging %d\n", lib->for_contigging);
    printf("LLL\tfor_splinting  %d\n", lib->for_splinting);
    printf("LLL\tfor_gapclosing %d\n", lib->for_gapclosing);
    printf("LLL\tono_setid      %d\n", lib->ono_setid);
    printf(KNORM);
}
