#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <time.h>
//#include <upc.h>
//#include <upc_tick.h>


#include "../common/MPIUtils.h"
//#include "../common/upc_common.h"
//#include "../common/utils.h" // TODO there may still be dependencies on this file
//#include "../common/timers.h"
#include "../common/common.h"
#include "../common/Buffer.h"

#include "exec_utils.h"

#define MAX_STAGES 12
#define MAX_ARGS 24

static double _stage_timings[MAX_STAGES];
static double _stage_mem[MAX_STAGES];
static double _stage_smaps_mem[MAX_STAGES];
static double _stage_peak_mem[MAX_STAGES];
static char *_stage_names[MAX_STAGES];
static int _nstages = 0;

//extern char *_use_col;

static int _dryrun = 0;

static char *get_current_time(char *tbuf)
{
    time_t t = time(NULL);
    strftime(tbuf, 255, "%D %T", localtime(&t));
    return tbuf;
}

int exec_stage(char *stages, int is_scaff, int (*stage_main)(int, char**), const char *stage_name, ...)
{
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
#ifdef TRACK_FILE_OPENS
    if (!MYTHREAD)
        fprintf(stderr, "%s\n", stage_name);
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
#endif
    if (stages) DBG("%s %d stages: %s \n", __FILE__, __LINE__, (stages? stages: "null") );
    //serial_printf("%s %d strstr(stages, \"-end\") %s", __FILE__, __LINE__,  strstr(stages, "-end"));
    if (stages && strstr(stages, "-end")) {
        // make sure we run all stages after this one
        if (strstr(stages, stage_name)) {
            serial_printf(">>>>>>>> Running stage %s to END <<<<<<<\n", stage_name);
            strcpy(stages, "all");
        }
    }
    if (stages && strcmp(stages, "all") != 0 && !strstr(stages, stage_name) &&
        !(strcmp(stages, "scaffolding") == 0 && is_scaff)) {
        LOGF("Skipping stage %s\n", stage_name);
        return 0;
    }
    // check for dryrun
    char *dryrun_str = getenv("DRYRUN");
    if (dryrun_str && atoi(dryrun_str) == 1)
        _dryrun = 1;
    char **args = (char**) calloc_chk(MAX_ARGS, sizeof(char*));
    Buffer *argBuffers = (Buffer*) malloc_chk(MAX_ARGS * sizeof(Buffer));
    for (int i = 0; i < MAX_ARGS; i++)
        argBuffers[i] = initBuffer(64);
    char *arg;
    int d;
    int64_t ld;
    double f;
    char *s;
    strcpyBuffer(argBuffers[0], stage_name);
    int nargs = 1;
    va_list ap;
    va_start(ap, stage_name);
    while ((arg = va_arg(ap, char*)) != NULL) {
        resetBuffer(argBuffers[nargs]); // in case it was previously a binary
        strcpyBuffer(argBuffers[nargs], arg);
        nargs++;
        char *fmt = strchr(arg, '%');
        if (fmt) {
            // strip out format characters
            truncateBuffer(argBuffers[nargs - 1], 2);
            switch (fmt[1]) {
            case 'B': // this is a binary flag
                d = va_arg(ap, int);
                // there is only a single param
                nargs--;
                // now strip out even this one if the flag is off
                if (!d)
                    nargs--;
                break;
            case 'd':
                d = va_arg(ap, int);
                printfBuffer(argBuffers[nargs], fmt, d);
                break;
            case 'f':
                f = va_arg(ap, double);
                printfBuffer(argBuffers[nargs], fmt, f);
                break;
            case 's':
                s = va_arg(ap, char*);
                printfBuffer(argBuffers[nargs], fmt, s);
                break;
            case 'l':
                ld = va_arg(ap, int64_t);
                printfBuffer(argBuffers[nargs], "%ld", ld);
                break;
            default:
                SDIE("Error in format string for exec stage %s: %s (only support d, f, s, B, l)\n",
                     stage_name, arg);
                return 0;
            }
            nargs++;
        }
    }
    va_end(ap);
    serial_printf("\n");
    serial_printf("%s" KNORM "\n", HASH_BAR);
    serial_printf("# Starting stage ");
    for (int i = 0; i < nargs; i++) {
        args[i] = getStartBuffer(argBuffers[i]);
        serial_printf("%s ", args[i]);
    }
    char tbuff[255];
    //serial_printf("at %s" KNORM "\n", get_current_time(tbuff));
    serial_printf("%s" KNORM "\n\n", HASH_BAR);
    LOGF("Starting %s\n", stage_name);

    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    double start_mem_free = 0;
    double start_mem_used = get_used_mem_gb();
    if (!MYTHREAD)
        start_mem_free = get_free_mem_gb();
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    double start_t = MPI_Wtime();
    if (!_dryrun) {
        int res = stage_main(nargs, args);
        if (res == EXIT_FAILURE)
            SDIE("%s stage failed with error code %d\n", stage_name, res);
    }
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    _stage_timings[_nstages] = MPI_Wtime() - start_t; //elapsed time
    _stage_names[_nstages] = strdup(stage_name);

    fflush(stdout);
    fflush(stderr);

//    double mem_leaked = get_used_mem_gb() - start_mem_used;
//    double tot_mem_leaked = reduce_double(mem_leaked, UPC_ADD, SINGLE_DEST);
    // now get max peak memory of any thread
    //double max_peak_mem = (double)reduce_int(get_max_mem_usage_mb(), UPC_MAX, SINGLE_DEST) / 1024;
    
    if (!MYTHREAD) {
        double end_mem_free = get_free_mem_gb();
        serial_printf("\n");
        serial_printf("%s" KNORM "\n", HASH_BAR);
        serial_printf("# Finished %s in %.2f s at %s\n",
                      stage_name, _stage_timings[_nstages], get_current_time(tbuff));
//        serial_printf("%s# Memory not freed: full system %.3f GB (on node 0), application %.3f GB "
//                      "(all nodes)" KNORM "\n", _use_col, start_mem_free - end_mem_free, tot_mem_leaked);
        serial_printf("# Memory not freed: full system %.3f GB (on node 0) " KNORM "\n",
                      start_mem_free - end_mem_free);
        serial_printf("# Memory remaining on node 0 after %s: %.3f GB" KNORM "\n", stage_name, end_mem_free);
        //serial_printf("%s# Peak memory per thread (incl. shared): %.2f GB" KNORM "\n", _use_col, max_peak_mem);
        serial_printf("%s" KNORM "\n", HASH_BAR);
        _stage_mem[_nstages] = start_mem_free - end_mem_free;
//        _stage_smaps_mem[_nstages] = tot_mem_leaked;
        _stage_smaps_mem[_nstages] = 0;
    //    _stage_peak_mem[_nstages] = max_peak_mem;
    }
    _nstages++;

#if defined(__APPLE__) && defined(__MACH__)
    // no /proc on Mac
#else
    FILE *fp_cr = fopen("/proc/self/clear_refs", "w");
    fprintf(fp_cr, "5\n");
    fclose(fp_cr);
#endif
    
    for (int i = 0; i < MAX_ARGS; i++)
        freeBuffer(argBuffers[i]);
    free_chk(argBuffers);
    free_chk(args);
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    return 1;
}

void print_timings(cfg_t *cfg)
{
	/*
    double tot_t = 0;
    for (int i = 0; i < _nstages; i++) 
        tot_t += _stage_timings[i];
    serial_printf("\n");
    serial_printf("%s%s" KNORM "\n", _use_col, HASH_BAR);
    serial_printf("%s# Completed %d stages in %.2f s on %d threads over %d nodes" KNORM "\n",
                  _use_col, _nstages, tot_t, THREADS, THREADS/cfg->cores_per_node);
    serial_printf("%s#    %-40s%8s    %8s %s  %s" KNORM "\n", 
                  _use_col, "STAGE", "time (s)", "Unfreed mem (GB)", "", "Peak mem (incl. shared)");
    serial_printf("%s#    %-40s%8s    %8s %8s      %8s" KNORM "\n", 
                  _use_col, "", "", "sys/node", "app/all", "GB/thread");
    double tot_time = 0;
    double tot_mem = 0;
    double tot_smaps_mem = 0;
    for (int i = 0; i < _nstages; i++) {
        serial_printf("%s#    %-40s%8.2f    %8.3f %8.3f      %8.3f" KNORM "\n", 
                      _use_col, _stage_names[i], _stage_timings[i], _stage_mem[i],
                      _stage_smaps_mem[i], _stage_peak_mem[i]);
        tot_time += _stage_timings[i];
        tot_mem += _stage_mem[i];
        tot_smaps_mem += _stage_smaps_mem[i];
    }
    serial_printf("%s#    %-40s%8.2f    %8.3f %8.3f" KNORM "\n", 
                  _use_col, "TOTAL", tot_time, tot_mem, tot_smaps_mem);
     */
}

char *get_str_from_file(const char *fname)
{
    if (_dryrun)
        return strdup("");
    FILE *f = fopen_track(fname, "r");
    if (!f)
        return NULL;
    char *val = (char*) calloc_chk(1, 1024);
    fscanf(f, "%s", val);
    fclose_track(f);
    return val;
}

/*
int64_t get_long_from_file_and_broadcast(const char *fname)
{
    int64_t val;
    if (!MYTHREAD) {
        char * str = get_str_from_file(fname);
        if (str) 
            val = atol(str);
        else
            val = 0;
        free_chk(str);
        broadcast_long(val, 0);
    } else {
        val = broadcast_long(-1, 0);
    }
    return val;
}
*/

void put_str_in_file(const char *fname, char *s)
{
    if (_dryrun)
        return;
    FILE *f = fopen_track(fname, "w");
    if (!f)
        SDIE("Could not write to file %s\n", fname);
    fprintf(f, "%s\n", s);
    fclose_track(f);
}

int64_t get_num_from_file(const char *prepend, const char *fname)
{
    if (_dryrun)
        return 0;

    char full_fname[1024];
    if (prepend[0]) {
        sprintf(full_fname, "%s-n%s.txt", prepend, fname);
        get_rank_path(full_fname, -1);
    } else {
        sprintf(full_fname, "n%s.txt", fname);
        get_rank_path(full_fname, -1);
    }
    char *val = get_str_from_file(full_fname);
    if (!val)
        SDIE("Error trying to read value from file %s\n", full_fname);
    int64_t num = atol(val);
    free_chk(val);
    return num;
}

/*
int64_t get_num_from_file_and_broadcast(const char *prepend, const char *fname)
{
    int64_t num;
    if (!MYTHREAD) {
        num = get_num_from_file(prepend, fname);
        broadcast_long(num, 0);
    } else {
        num = broadcast_long(-1, 0);
    }
    return num;
}
*/

int dummy_exec(char *stage_name)
{
	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    struct stat st;
    char buf[256];
    sprintf(buf, "%s-dummy.txt", stage_name);
    if (stat(buf, &st) == 0)
        return 0;
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    if (!MYTHREAD) {
        FILE *f = fopen_chk(buf, "w");
        fclose(f);
        serial_printf("\n");
        serial_printf("%s" KNORM "\n", HASH_BAR);
        char tbuff[255];
        serial_printf("# Starting %s (dummy) at %s" KNORM "\n",
                      stage_name);//, get_current_time(tbuff));
        serial_printf("%s" KNORM "\n\n", HASH_BAR);
    }
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    return 1;
}
