#ifndef __STATIC_VARS_H
#define __STATIC_VARS_H

typedef struct {
    FILE *_my_log;
    FILE *_my_diags;
    FILE *_dbg_file;
    double fileIOTime, cardCalcTime, setupTime, storeTime;
    double barrierTime;
    int numBarriers;
    int _logTime;
    int _num_warnings;
    int cores_per_node;
    char *LOG_PREFIX;
    char flushLog;
    char logMemcheck;
    char _logMemcheckLast;
    int64_t outstandingMallocs;
    int64_t untracked_outstandingMallocs;
    int64_t outstandingUPCAllocs;
    int64_t untracked_outstandingUPCAllocs;
    int64_t outstandingUPCAllAllocs;
} _StaticVars;
extern _StaticVars *_sv; // This must be instantiated and initialized (NULL ok) by the primary executable
#define MYSV _sv[MYTHREAD]
#define LOG_SET_FLUSH(x) do { if (_sv) { MYSV.flushLog = x; } } while(0)
#define LOG_SET_MEMCHECK(x) do { \
    if (_sv) { \
        char tmp = MYSV._logMemcheckLast; \
        if (MYSV.logMemcheck) { MYSV._logMemcheckLast = MYSV.logMemcheck; } \
        MYSV.logMemcheck = ((x)<0) ? tmp : (x); \
    } \
} while (0)

#ifdef DEBUG
#define LOG_RESET_MEMCHECK LOG_SET_MEMCHECK(1)
#else
#define LOG_RESET_MEMCHECK LOG_SET_MEMCHECK(0)
#endif
#define MEMCHECK_COUNT0(x) do { if (_sv) { MYSV.untracked_outstandingMallocs += x; } } while(0)
#define MEMCHECK_COUNT(x) do { if (_sv) { MYSV.outstandingMallocs += x; } } while(0)
#define MEMCHECK_UPCCOUNT0(x) do { if (_sv) { MYSV.untracked_outstandingUPCAllocs += x; } } while(0)
#define MEMCHECK_UPCCOUNT(x) do { if (_sv) { MYSV.outstandingUPCAllocs += x; } } while(0)
#define MEMCHECK_UPCALLCOUNT(x) do { if (_sv) { MYSV.outstandingUPCAllAllocs += x; } } while(0)

/* set up StaticVars */

#if 0 // enable when MPI is gone...

  typedef shared[1] _StaticVars *StaticVars;

  #define INIT_STATIC_VARS \
    do { \
      assert(_sv == NULL) ; \
      UPC_ALL_ALLOC_CHK(_sv, THREADS, sizeof(_StaticVars)); \
      memset((char*) &_sv[MYTHREAD], 0, sizeof(_StaticVars)); \
      assert(MYSV._my_log == NULL); \
      MYSV.cores_per_node = get_cores_per_node(); \
      assert(MYSV.cores_per_node <= THREADS); \
      LOG_RESET_MEMCHECK; \
    } while (0)

  #define FREE_STATIC_VARS \
    do { \
      if (_sv != NULL) { \
          upc_barrier; \
          UPC_ALL_FREE(_sv); \
          _sv = NULL; \
      } \
    } while (0)

#else
  typedef _StaticVars *StaticVars;
  #if defined(__BERKELEY_UPC__) && (__BERKELEY_UPC_PTHREADS__ == 1)
    #ifndef HIPMER_PTHREADS 
    #error "__BERKELEY_UPC_PTHREADS__ is defined but HIPMER_PTHREADS is not!\n"
    #endif

    #define INIT_STATIC_VARS \
      do { \
        int cpn = get_cores_per_node(); \
        if (MYTHREAD % cpn == 0) { \
            assert(_sv == NULL); \
            _sv = (StaticVars) calloc_chk(cpn, sizeof(_StaticVars)); \
            for(int i = 0; i < cpn; i++) _sv[i].cores_per_node = cpn; \
            _sv -= MYTHREAD; \
        } \
        upc_barrier; \
        if (_sv == NULL) DIE("Could not initialize static variables!\n"); \
        assert(MYSV._my_log == NULL && MYSV.cores_per_node = cpn); \
        LOG_RESET_MEMCHECK; \
      } while (0)

    #define FREE_STATIC_VARS \
      do { \
        if (_sv == NULL) { \
            upc_barrier; \
        } else if (MYTHREAD % MYSV.cores_per_node == 0) { \
            upc_barrier; /* wait for all threads to get past this test */ \
            _sv += MYTHREAD; \
            free_chk(_sv); \
            _sv = NULL; \
        } else { \
            upc_barrier; /* past this test */ \
        } \
        assert(_sv == NULL); \
      } while(0)

  #else // NO PTHREADS
   
    #define INIT_STATIC_VARS \
    do { \
        assert(_sv == NULL); \
        _sv = (StaticVars) calloc_chk(1, sizeof(_StaticVars)); \
        assert(_sv != NULL); \
        _sv -= MYTHREAD; /* backtrack to make x[MYTHREAD] valid in this thread... */ \
        assert(MYSV._my_log == NULL); \
        MYSV.cores_per_node = get_cores_per_node(); \
        LOG_RESET_MEMCHECK; \
    } while (0)
 
    #define FREE_STATIC_VARS \
    do { \
       if (_sv != NULL) { \
           _sv += MYTHREAD; \
           free_chk(_sv); \
           _sv = NULL; \
       } \
    } while (0)

  #endif
#endif

#define SET_CORES_PER_NODE(cpn) \
 do { \
   if (!_sv) DIE("Can not SET_CORES_PER_NODE before _sv has been initialized!\n"); \
   int oldcpn = MYSV.cores_per_node; \
   if (oldcpn != 0 && oldcpn != cpn) SWARN("cores_per_node has a discrepancy!  Discovered %d but user overrides with %d\n", oldcpn, cpn); \
   MYSV.cores_per_node = cpn; \
   CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) ); \
 } while (0)

#define ADD_BARRIER_TIME(time) do { \
   if (_sv) { \
      MYSV.barrierTime += time; \
      MYSV.numBarriers++; \
   } \
 } while (0)

#define GET_BARRIER_TIME() ( _sv ? MYSV.barrierTime : 0.0 )

#endif // __STATIC_VARS_H
