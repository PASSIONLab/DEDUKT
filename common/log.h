#ifndef __LOG_H
#define __LOG_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdint.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <time.h>

#include "../common/colors.h"
#include "../common/mpi_common.h"
// INCLUDE other common header files too

//#include "version.h"
#include "../common/StaticVars.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __FILENAME__
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#ifndef EXIT_FUNC
#define EXIT_FUNC(x) do { common_exit(x); } while (0)
#endif

#ifndef  __UPC_VERSION__
#define MYTHREAD get_rank()
#define THREADS get_num_ranks()
#endif



/* set up LOG, SLOG, DBG, DBG2, DIE, SDIE, WARN, etc macros */

/* verbosity levels: */
typedef enum _LOG_LEVELS {
                    LL_DIE = 0,  /*  FATAL/ERROR/DIE write to file and close file and print to stderr and flush and exit program */
                    LL_WARN, /*   WARN to file and print to stderr and flush */
                    LL_LOG,  /*   LOG/INFO to file and print to stderr */
                    LL_LOGF, /*   LOG/INFO just log to file */
                    LL_DBG,  /*   DEBUG1 to file only */
                    LL_DBG2, /*   DEBUG2 to file only */
                    LL_DBG3, /*   DEBUG3 to file only */
                    LL_MAX} _LOG_LEVELS_TYPE;
/*
   SLOG   only !MYTHREAD outputs ... all log if there is a file
   SDIE   only !MYTHREAD outputs ... all terminate
   SWARN  only !MYTHREAD outputs ... all log if there is a file

   everything >= HIPMER_VERBOSITY gets either printed to stderr or writen to _my_log or both
   everything <= 1 gets printed and flushed to stderr
*/

#ifdef HIPMER_VERBOSE
  #define HIPMER_VERBOSITY (HIPMER_VERBOSE + LL_WARN)
#else
  #ifndef DEBUG
    #define HIPMER_VERBOSITY LL_LOGF
  #else
    #define HIPMER_VERBOSITY LL_DBG
  #endif
#endif

inline time_t setLogPrefix(char *_LOG_prefix, int size, _LOG_LEVELS_TYPE level, int _printPrefix, const char *fmt, const char *filename, int line) {
   assert(size > 1);
   _LOG_prefix[0] = '\0';
   time_t _LOG_ltime = time(NULL);
   if (_printPrefix) {
     struct tm _LOG_result; char _LOG_stime[32];
     localtime_r(&_LOG_ltime, &_LOG_result); asctime_r(&_LOG_result, _LOG_stime);
     _LOG_stime[strlen(_LOG_stime) - 1] = '\0';
     char _LOG_logLabel[32];
     strcpy(_LOG_logLabel, "[INFO]");
     if (level >= LL_DBG) { snprintf(_LOG_logLabel, 32, "%s[DBG%d]%s", KLBLUE, (int)(level-LL_DBG+1), KNORM); }
     if (level < LL_LOG)  { snprintf(_LOG_logLabel, 32, "%s[%s]%s", (level <= LL_DIE ? KLRED : KRED), (level <= LL_DIE ? "DIE" : "WARN"), KNORM); }
     snprintf(_LOG_prefix, size, "Thread %d %s %s [%s:%d-%s]: ", (int) MYTHREAD, _LOG_logLabel, _LOG_stime, filename, (int) line, NULL);
    }
    return _LOG_ltime;
}
#define NO_LOG_PREFIX 256
#define _LOG(level, fmt, ...) do { \
   enum _LOG_LEVELS  _level = (enum _LOG_LEVELS) ((level) & ~(NO_LOG_PREFIX)); \
   int _printPrefix = ((level) & (NO_LOG_PREFIX)) == 0; \
   if (HIPMER_VERBOSITY >= _level) {                                        \
       char _LOG_prefix[192]; \
       int _haveLogFile = 0; int _flushLogFile = 0; \
       if (_sv) { if (MYSV._my_log) _haveLogFile = 1; if (MYSV.flushLog > 0 && (MYTHREAD % MYSV.flushLog) == 0) _flushLogFile=1; } \
       time_t _startLog = setLogPrefix(_LOG_prefix, 192, _level, _printPrefix, fmt, __FILENAME__, __LINE__); \
       if (_haveLogFile) { \
          fprintf(MYSV._my_log, "%ds %s" fmt, MYSV._logTime, _LOG_prefix, ##__VA_ARGS__); \
          if ( _level < LL_LOG || HIPMER_VERBOSITY >= LL_DBG || ((MYTHREAD) == 0) || _flushLogFile > 0 ) fflush(MYSV._my_log); \
       } \
       if (_level < LL_LOG) { \
          fprintf(stderr, "%s" fmt, _LOG_prefix, ##__VA_ARGS__); \
          fflush(stderr); \
       } else if (_level == LL_LOG || ((!_haveLogFile) && (!MYTHREAD))) { \
          fprintf(stdout, "%s" fmt, _LOG_prefix, ##__VA_ARGS__); \
          fflush(stdout); \
       } \
       if (_sv) { MYSV._logTime += time(NULL) - _startLog; } \
   } \
 } while(0)

#define LOG(fmt, ...)  _LOG(LL_LOG, fmt, ##__VA_ARGS__)
#define LOGN(fmt, ...)  _LOG((LL_LOG | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#define LOGF(fmt, ...) _LOG(LL_LOGF, fmt, ##__VA_ARGS__)
#define LOGFN(fmt, ...) _LOG((LL_LOGF | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#define TLOGF(test, fmt, ...) _LOG( ((test) ? LL_LOG : LL_LOGF), fmt, ##__VA_ARGS__)
#define TLOGFN(test, fmt, ...) _LOG( (((test) ? LL_LOG : LL_LOGF) | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#define SLOG(fmt, ...) _LOG((MYTHREAD == 0 ? LL_LOG : ( ((_sv != NULL) && (MYSV._my_log != NULL)) ? LL_LOGF : LL_DBG2 )), fmt, ##__VA_ARGS__)
#define SLOGN(fmt, ...) _LOG(((MYTHREAD == 0 ? LL_LOG : ( ((_sv != NULL) && (MYSV._my_log != NULL)) ? LL_LOGF : LL_DBG2 )) | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#define SDBG(fmt, ...) if (MYTHREAD == 0) { _LOG(LL_DBG, fmt, ##__VA_ARGS__); }
#define SDBGN(fmt, ...) if (MYTHREAD == 0) { _LOG((LL_DBG | NO_LOG_PREFIX), fmt, ##__VA_ARGS__); }

#define DBG(fmt, ...)  _LOG(LL_DBG, fmt, ##__VA_ARGS__)
#define DBGN(fmt, ...)  _LOG((LL_DBG | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#define DBG1(fmt, ...) _LOG(LL_DBG, fmt, ##__VA_ARGS__)
#define DBG1N(fmt, ...) _LOG((LL_DBG | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#define DBG2(fmt, ...) _LOG(LL_DBG2, fmt, ##__VA_ARGS__)
#define DBG2N(fmt, ...) _LOG((LL_DBG2 | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#define DBG3(fmt, ...) _LOG(LL_DBG3, fmt, ##__VA_ARGS__)
#define DBG3N(fmt, ...) _LOG((LL_DBG3 | NO_LOG_PREFIX), fmt, ##__VA_ARGS__)

#ifndef REDUCED_LOG_COUNT
#ifdef DEBUG
#define REDUCED_LOG_COUNT (THREADS)
#else
#define REDUCED_LOG_COUNT 128
#endif
#endif

#ifndef REDUCED_LOG_MODULUS
#define REDUCED_LOG_MODULUS ((THREADS+(REDUCED_LOG_COUNT-1))/(REDUCED_LOG_COUNT)) /* only up to REDUCED_LOG_COUNT (+3) logs are created */
#endif

#define OPEN_MY_LOG(fmt, ...) \
    do { \
        if (!_sv) INIT_STATIC_VARS; \
        assert(_sv); \
        if ((MYTHREAD % (REDUCED_LOG_MODULUS)) != 0 && \
            MYTHREAD != 0 && MYTHREAD != (THREADS-1) && MYTHREAD != THREADS/2 \
           ) break; /* only some threads will have a log: 0, mid, last & some capped modulus */ \
        char myLog_prefix[MAX_FILE_PATH], myLog_logfile_name[MAX_FILE_PATH];\
        snprintf(myLog_prefix, MAX_FILE_PATH, fmt, ##__VA_ARGS__); \
        snprintf(myLog_logfile_name, MAX_FILE_PATH, "%s-%d.log",  myLog_prefix, MYTHREAD); \
        char *rankpath = get_rank_path(myLog_logfile_name, MYTHREAD);\
        MYSV._my_log = fopen(rankpath, "a");                                \
        if (!MYSV._my_log) {\
            fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: Could not open file '%s': %s\n" KNORM, \
                    MYTHREAD, __FILENAME__, __LINE__, ""/*HIPMER_VERSION*/, rankpath, strerror(errno));\
            EXIT_FUNC(1);\
        }                \
        LOGF("Opened this file: %s.\n", rankpath); \
        fflush(MYSV._my_log); \
    } while (0)

#define CLOSE_MY_LOG  \
        do { \
            if (_sv) { if (MYSV._my_log) { \
                SLOG("Outstanding Memory: mallocs=%lld \tupc_all_alloc=%lld \tupc_alloc=%lld.  Untracked mallocs=%lld upc_allocs=%lld\n", (long long int) MYSV.outstandingMallocs, (long long int) MYSV.outstandingUPCAllAllocs, (long long int) MYSV.outstandingUPCAllocs, (long long int) MYSV.untracked_outstandingMallocs, (long long int) MYSV.untracked_outstandingUPCAllocs); \
                LOGF("Closing this file.\n"); \
                fclose_track(MYSV._my_log); \
                MYSV._my_log = NULL; \
            } } \
            fflush(stderr); fflush(stdout); \
        } while (0)

#define DIE(fmt,...) \
    do { \
        _LOG(LL_DIE, fmt, ##__VA_ARGS__); \
        CLOSE_MY_LOG; \
        EXIT_FUNC(1); \
    } while (0)

#define SDIE(fmt,...) \
    do { \
        if (!MYTHREAD) \
            _LOG(LL_DIE, fmt, ##__VA_ARGS__); \
        CLOSE_MY_LOG; \
        EXIT_FUNC(1); \
    } while (0)

#define WARN(fmt,...) \
    do { \
        if (_sv) MYSV._num_warnings++;\
        _LOG(LL_WARN, fmt, ##__VA_ARGS__); \
    } while (0)

#define SWARN(fmt,...) \
    do { \
        if (_sv && !MYTHREAD) MYSV._num_warnings++;\
        _LOG( (!MYTHREAD ? LL_WARN : (((_sv) && MYSV._my_log) ? LL_LOGF : LL_DBG) ), fmt, ##__VA_ARGS__); \
    } while (0)


#define serial_printf(fmt, ...)                     \
    do {                                            \
        if (MYTHREAD == 0) {                        \
            fprintf(stdout, fmt, ##__VA_ARGS__);    \
            fflush(stdout);                         \
            if (_sv != NULL) { if(MYSV._my_log != NULL) { if (strchr(fmt, '\n')) { LOGF(fmt, ##__VA_ARGS__); } else { LOGFN(fmt, ##__VA_ARGS__); } } }  \
        }                                           \
    } while (0)

#define sprintf_chk(s, max, fmt, ...) \
    do { \
        int linelen = snprintf(s, max, fmt, ##__VA_ARGS__); \
        if (linelen >= max) \
           WARN("Thread %d, buffer overflow printing '%s' at %s:%d-%s: (" fmt ")\n", MYTHREAD, fmt,  __FILENAME__, __LINE__, ""/*HIPMER_VERSION*/, ##__VA_ARGS__); \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif // __LOG_H

