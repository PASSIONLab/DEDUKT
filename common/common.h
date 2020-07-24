#ifndef __COMMON_H
#define __COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>
#ifdef __APPLE__
  #include <sys/param.h>
  #include <sys/mount.h>
#else
  #include <sys/vfs.h>
#endif
#include <sys/types.h>
#include <dirent.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdint.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <time.h>
#include <zlib.h>

typedef long long int lld;
typedef unsigned long long int llu;

//#include "version.h"
#include "defines.h"

#include "optlist.h"

#include "colors.h"
#include "StaticVars.h"

#include "common_rankpath.h"
#include "memory_chk.h"

#include "log.h"
#include "file.h"
#include "mach.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CHECK_ERR(cmd)                                                  \
    do {                                                                \
        int err;                                                        \
        if ((err = cmd) != 0)                                           \
            DIE("Thread %d, [%s:%d-%s]: " #cmd " failed, error %d: %s\n", MYTHREAD, __FILENAME__, __LINE__, ""/*HIPMER_VERSION*/, err, strerror(err)); \
    } while (0)


#define ASSERT_BOUNDS(x, bound) do { \
       assert(x >= 0); \
       assert(x < bound); \
   } while (0)

#ifdef DEBUG

  #define CHECK_BOUNDS(x, bound) do {                                     \
      if (((x) >= (bound)) || (((lld) (x)) < 0)) {                               \
         lld x1 = x, bound1 = bound;                                    \
         DIE("index out of range for %s >= %s: %lld > %lld\n", #x, #bound, x1, bound1); \
      } \
   } while (0)

  #define CHECK_UPPER_BOUND(x, bound) do {                                \
            if ((x) >= (bound)) {                                       \
                lld x1 = x, bound1 = bound;                             \
                DIE("index out of range for %s >= %s: %lld > %lld\n", #x, #bound, x1, bound1); \
            }                                                           \
        } while (0)
#else

  #define CHECK_BOUNDS(x, bound)
  #define CHECK_UPPER_BOUND(x, bound)

#endif

#define INT_CEIL(numerator, denominator) (((numerator) - 1) / (denominator) + 1);

static off_t get_file_size(const char *fname)
{
    struct stat s;
    if (stat(fname, &s) != 0) {
        WARN("could not stat %s: %s\n", fname, strerror(errno));
        return 0;
    }
    return s.st_size;
}

static int64_t __atoi64 (const char *nptr)
{
   int c;
   int64_t value;
   int sign;

   while (isspace((int)*nptr))
      ++nptr;

   c = (int)*nptr++;
   sign = c;
   if (c == '-' || c == '+')
      c = (int)*nptr++;

   value = 0;

   while (isdigit(c))
   {
      value = 10 * value + (c - '0');
      c = (int)*nptr++;
   }

   if (sign == '-')
      return -value;
   else
      return value;
}

/* Memory helpers (from memory_chk.h too */

#ifndef LOG_ALL_MEM
#ifdef DEBUG
#define LOG_ALL_MEM ( (MYTHREAD==0 || MYTHREAD==THREADS/2 || MYTHREAD==THREADS-1) && _sv != NULL && MYSV._my_log != NULL ) /* just first, middle and last threads and only if a log is defined */
#else
#define LOG_ALL_MEM 0
#endif
#endif

#define _free_chk(p) free_chk_at_line((void**) p, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 1)
#define free_chk(p) _free_chk(&(p))
static inline void free_chk_at_line(void **_ptr, const char *which_file, int which_line, int PAD_ALLOC, int log)
{
	if (!_ptr) DIE("[%s:%d-%s]: free_chk: called with invalid pointer to pointer!\n", which_file, which_line, "" /*HIPMER_VERSION*/);
	void *ptr = *_ptr;
	if (ptr) {
		char msg[256];
		size_t s = 0;
		if (PAD_ALLOC > 0) {
			ptr = pad_memory_check_at_line(ptr, msg, PAD_ALLOC, which_file, which_line);
		}
		if (ptr == NULL) {
		    DIE("[%s:%d-%s] Request to free a pointer, but padding is corrupt! %s\n", which_file, which_line, "" /*HIPMER_VERSION*/, msg);
        }
		s = pad_memory_real_size(*_ptr, PAD_ALLOC);
		if (log) { MEMCHECK_COUNT( - (PAD_ALLOC ? s : 1) ); } else { MEMCHECK_COUNT0( - (PAD_ALLOC ? s : 1) ); }
		if (log && _sv && LOG_ALL_MEM && MYSV.logMemcheck > 0) {
			DBGN("[%s:%d] %lld %lld free_chk(%p) realPtr=%p PAD=%llu size=%llu\n", which_file, which_line, (lld) MYSV.outstandingMallocs, (lld) MYSV.untracked_outstandingMallocs, *_ptr, ptr, (llu) PAD_ALLOC, (llu) s);
		}
		if (&MYSV != NULL && MYSV.outstandingMallocs < 0) WARN("[%s-%d] more freed than malloc: %lld\n", which_file, which_line, (lld) MYSV.outstandingMallocs);
		free(ptr);
	} else {
		WARN("[%s:%d-%s]: free_chk called on NULL pointer!\n", which_file, which_line, "" /*HIPMER_VERSION*/);
	}
	*_ptr = NULL;
}

#define malloc_chk(s) malloc_chk_at_line(s, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 1)
static inline void *malloc_chk_at_line(size_t _size, const char *which_file, int which_line, int PAD_ALLOC, int log)
{
	size_t size = _size;
	if (!size) {
		LOGF("[%s:%d-%s]: Request to malloc %llu bytes\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) size);
		return NULL;
	}
	size += PAD_ALLOC*2;
	void *_x = malloc(size);
	void *x = _x;
	if (!x) {
		DIE("[%s:%d-%s]: Could not malloc %llu bytes\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) size);
	}
	if (PAD_ALLOC > 0) {
		x = pad_memory_at_line(x, size, PAD_ALLOC, which_file, which_line);
	}
	if (log) { MEMCHECK_COUNT( PAD_ALLOC ? size : 1); } else { MEMCHECK_COUNT0( PAD_ALLOC ? size : 1); }
	if (log && (_sv) && LOG_ALL_MEM && MYSV.logMemcheck > 0) DBGN("[%s:%d] %lld malloc_chk(%llu) %p realPtr=%p realSize=%llu PAD=%llu\n", which_file, which_line, (lld) MYSV.outstandingMallocs, (llu) _size, x, _x, (llu) size, (llu) PAD_ALLOC);
	return x;
}

#define calloc_chk(n, s) calloc_chk_at_line(n, s, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 1)
static inline void *calloc_chk_at_line(size_t nmemb, size_t size, const char *which_file, int which_line, int PAD_ALLOC, int log)
{
	if (!size) {
		LOGF("[%s:%d-%s]: Request to calloc %llu bytes\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) size);
		return NULL;
	}
	if (!nmemb) {
		DIE("[%s:%d-%s]: Request to calloc %llu elements\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) nmemb);
	}
	size_t realSize = nmemb * size + 2*PAD_ALLOC;
	void *_x = malloc(realSize);
	void *x = _x;
	if (!x) {
		DIE("[%s:%d-%s]: Could not calloc %llu bytes\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) size * nmemb);
	}
	memset(x, 0, realSize);
	if (PAD_ALLOC > 0) {
		x = pad_memory_at_line(x, realSize, PAD_ALLOC, which_file, which_line);
	}
	if (log) { MEMCHECK_COUNT( PAD_ALLOC ? realSize : 1); } else { MEMCHECK_COUNT0( PAD_ALLOC ? realSize : 1); }
	if (log && _sv && LOG_ALL_MEM && MYSV.logMemcheck > 0) DBGN("[%s:%d] %lld calloc_chk(%llu, %llu) %p realPtr=%p realSize=%llu PAD=%llu\n", which_file, which_line, (lld) MYSV.outstandingMallocs, (llu) nmemb, (llu) size, x, _x, (llu) realSize, (llu) PAD_ALLOC);
	return x;
}

#define realloc_chk(p, s) realloc_chk_at_line(p, s, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 1)
static inline void *realloc_chk_at_line(void *_ptr, size_t _size, const char *which_file, int which_line, int PAD_ALLOC, int log)
{
	size_t size = _size;
        void *ptr = _ptr;
	if (ptr != NULL) {
		char msg[256];
		if (PAD_ALLOC > 0)
			ptr = pad_memory_check_at_line(ptr, msg, PAD_ALLOC, which_file, which_line);
		if (ptr == NULL) {
		    DIE("[%s:%d-%s]: Request to realloc %llu bytes, but padding is corrupt! %s\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) size, msg);
                }
		size_t oldSize = PAD_ALLOC ? pad_memory_real_size(_ptr, PAD_ALLOC) : 0;
		if (log) { MEMCHECK_COUNT( - (PAD_ALLOC ? oldSize : 1) ); } else { MEMCHECK_COUNT0( - (PAD_ALLOC ? oldSize : 1) ); }
	}
	if (!size) {
		LOGF("[%s:%d-%s]: Request to realloc %llu bytes\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) size);
		if (ptr) free_chk_at_line(&ptr, which_file, which_line, PAD_ALLOC, log);
		return NULL;
	}
	size += PAD_ALLOC * 2; 
	void *_x = realloc(ptr, size);
	void *x = _x;
	if (!x) {
		DIE("[%s:%d-%s]: Could not realloc %llu bytes\n", which_file, which_line, "" /*HIPMER_VERSION*/, (llu) size);
	}
	if (PAD_ALLOC > 0) {
		x = pad_memory_at_line(x, size, PAD_ALLOC, which_file, which_line);
	}
	if (log) { MEMCHECK_COUNT( PAD_ALLOC ? size : 1); } else { MEMCHECK_COUNT0( PAD_ALLOC ? size : 1); }
	if (log && _sv && LOG_ALL_MEM && MYSV.logMemcheck > 0) DBGN("[%s:%d] %lld realloc_chk(%llu) %p -> %p realPtrs: %p -> %p realSize=%llu PAD=%llu\n", which_file, which_line, (lld) MYSV.outstandingMallocs, (llu) _size, _ptr, x, ptr, _x, (llu) size, (llu) PAD_ALLOC);
	return x;
}


#define strndup_chk(s, n) strndup_chk_at_line(s, n, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 1)
#define strdup_chk(s) strndup_chk_at_line(s, strlen(s), __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 1)
static inline char *strndup_chk_at_line(const char * src, size_t n, const char *which_file, int which_line, int PAD_ALLOC, int log)
{
	if (log && _sv && LOG_ALL_MEM && MYSV.logMemcheck > 0) DBGN("[%s:%d] strdup(%lld)\n", which_file, which_line, (lld) n);
	char *ptr = (char*) malloc_chk_at_line( (n * sizeof(char))+1, which_file, which_line, PAD_ALLOC, log);
	if ( !strncpy(ptr, src, n*sizeof(char)) ) DIE("[%s:%d]: Could not strncpy(%lld)\n", which_file, which_line, (lld) n);
	ptr[n] = '\0';
	return ptr;
}

#ifndef FULL_MEMCHECK
  #define malloc_chk0(s) malloc_chk_at_line(s, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 0)
  #define calloc_chk0(n, s) calloc_chk_at_line(n, s, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 0)
  #define realloc_chk0(p, s) realloc_chk_at_line(p, s, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 0)
  #define strndup_chk0(s, n) strndup_chk_at_line(s, n, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 0)
  #define strdup_chk0(s) strndup_chk_at_line(s, strlen(s), __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 0)
  #define _free_chk0(p) free_chk_at_line((void**) p, __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 0)
  #define free_chk0(p) free_chk_at_line((void**) &(p), __FILENAME__, __LINE__, PAD_ALLOC_BYTES, 0)
#else
  #define malloc_chk0(s) malloc_chk(s)
  #define calloc_chk0(n, s) calloc_chk(n, s)
  #define realloc_chk0(p, s) realloc_chk(p, s)
  #define strndup_chk0(s, n) strndup_chk(s, n)
  #define strdup_chk0(s) strdup_chk(s)
  #define _free_chk0(p) _free_chk(p)
  #define free_chk0(p) free_chk(p)
#endif

/*
static size_t get_chunk_size(size_t element_size, int cores_per_node) {
    size_t chunk_size = (TARGET_CHUNK_PER_THREAD + element_size - 1) / element_size;
    while (chunk_size > MIN_CHUNK_SIZE && chunk_size * THREADS * element_size * cores_per_node > MAX_CHUNK_MEMORY_PER_NODE)
        chunk_size /= 2;
    assert(chunk_size >= MIN_CHUNK_SIZE /2 );
    assert(chunk_size >= 1);
    return chunk_size;
}
*/

static int purgeDir(const char *dir, lld *bytes, int *files, int *dirs) {
    uid_t uid = getuid(), euid = geteuid();
    DBG2("purgeDir(%s) as %d(%d)\n", dir, (int) uid, (int) euid);
    int count = 0, _files = 0, _dirs = 0, __files = 0, __dirs = 0;
    lld _bytes = 0, __bytes = 0;
    if (bytes == NULL) bytes = &__bytes;
    if (files == NULL) files = &__files;
    if (dirs == NULL) dirs = &__dirs;
    DIR *d = opendir(dir);
    // ignore permission denied errors - means directory owned by someone else
    if (!d && errno != 13) { 
        WARN("Could not opendir %s! %s\n", dir, strerror(errno)); 
        return count; 
    }

    //fprintf(stderr, "Opened dir %s\n", dir);
    struct dirent *result;
    int path_length;
    char path[MAX_FILE_PATH];
    struct stat s;
    while (( result = readdir(d)) != NULL ) {
        path_length = snprintf (path, MAX_FILE_PATH,
                                "%s/%s", dir, result->d_name);
        if (path_length >= MAX_FILE_PATH) {
            WARN("Path length has got too long.\n");
            continue;
        }
        stat(path, &s);
        if (s.st_uid != uid && s.st_uid != euid) {
            DBG2("Skipping %s/%s as I do not own it (%d)\n", dir, result->d_name, (int) s.st_uid);
            continue;
        }
        if(S_ISDIR(s.st_mode)) {

            _dirs++;
            if (result->d_name[0] != '.') { // skip '.', '..', and anything starting with '.'
                /* Recursively call purgeDir  */
                count += purgeDir(path, bytes, files, dirs);
            }

        } else if (S_ISREG(s.st_mode) || S_ISLNK(s.st_mode)) {
            DBG("Unlinking %s\n", path);
            _files++;
            if(unlink(path) == 0) {
              _bytes += s.st_size;
              count++;
            }
        }
    }
    if (closedir(d) == 0) {
        if (strcmp(dir, "/dev/shm") == 0 || strcmp(dir, "/dev/shm/") == 0 || strcmp(dir, "/dev/shm/.") == 0) {
            DBG("Finished purgeDir %s\n", dir);
        } else {
            DBG("rmdir %s\n", dir);
            if (rmdir(dir) != 0) WARN("Could not rmdir(%s)\n", dir);
        }
    }
    
    *bytes += _bytes;
    *files += _files;
    *dirs += _dirs;
    DBG("Exiting %s: purged %lld bytes, %d files, %d dirs (total %lld %d %d)\n", dir, (lld) _bytes, (int) _files, (int) _dirs, (lld) *bytes, (int) *files, (int) *dirs);
    return count;
}

static inline int isACGT(const char base) {
    if (base == 'A' || base == 'C' || base == 'G' || base == 'T') return 1;
    return 0;
}
static inline int isACGTN(const char base) {
    if (base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N') return 1;
    DBG("isACGTN(%c) is not ACGT or N\n", base);
    return 0;
}

static inline int check_seqs(char *seq, char *label, int64_t *num_ambig)
{
    int len = strlen(seq);
	for (int i = 0; i < len; i++) {
        char c = toupper(seq[i]);
        switch (c) {
        case 'A': case 'C': case 'G': case 'T': case 'N': case 'U':
            break;
        case 'R': //seq[i] = 'A'; break;
        case 'Y': //seq[i] = 'C'; break;
        case 'K': //seq[i] = 'T'; break;
        case 'M': //seq[i] = 'C'; break;
        case 'S': //seq[i] = 'G'; break;
        case 'W': //seq[i] = 'T'; break;
        case 'B': //seq[i] = 'C'; break;
        case 'D': //seq[i] = 'A'; break;
        case 'H': //seq[i] = 'T'; break;
        case 'V': //seq[i] = 'G'; break;
        case '-':
            // could be IUB/IUPAC coding with ambigous characters
            // if so, convert to an N. 
            (*num_ambig)++;
            seq[i] = 'N';
            break;
        default:
			DIE("Invalid base %c %d in sequence (len %d) %.80s\n%s\n", c, (int) c, len, seq, label);
        }
	}
    return len;
}

/*
static inline int check_seq_with_len(const char *seq, int len, char *label)
{
	for (int i = 0; i < len; i++) {
        char c = toupper(seq[i]);
		if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
			DIE("Invalid base %c %d in sequence at pos %d (len %d) %.80s\n%s\n", c, c, i, len, seq, label);
        }
	}
    return len;
}
*/

static inline char *createSHM(char *my_fname, int64_t fileSize) {
    // create shm file
    char *ret = NULL;
    int fd = shm_open(my_fname, O_CREAT|O_TRUNC|O_RDWR, S_IRUSR|S_IWUSR);
    if (fd != -1) {
      if (ftruncate(fd, fileSize) == -1) {
        WARN("Failed to truncate to %lld %s: %s\n", (lld) fileSize, my_fname, strerror(errno));
      } else if (fileSize > 0) {
        LOGF("Opened shm file %s for writing of size %lld\n", my_fname, (lld) fileSize);
        void *tmpaddr = mmap(NULL, fileSize, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
        if (tmpaddr == MAP_FAILED) {
            WARN("Could not mmap %s for writing: %s, fileSize %lld\n", my_fname, strerror(errno), (lld) fileSize);
        } else {
            ret = (char*) tmpaddr;
        }
      }
      close(fd);
    }
    return ret;
}


#define UPC_MEMGET_STR(dst, src, len) do {                              \
        upc_memget(dst, src, len);                                      \
        if (strlen(dst) != len - 1)                                     \
            DIE("Length mismatch when getting remote NULL terminated string: %lld != %lld\n", (lld) strlen(dst), (lld) len - 1); \
    } while (0)

static size_t printFoldedSequence(GZIP_FILE out, const char *finalSequence, size_t cur_length) {
    /* Print contig in FASTA FORMAT */
    assert(SEGMENT_LENGTH < 8192);
    size_t total_written = 0, towrite;
    const char *seqF = (char*) finalSequence;
    char fastaSegment[SEGMENT_LENGTH];
    while ( total_written < cur_length ) {
        if (total_written + SEGMENT_LENGTH-1 < cur_length) {
            towrite = SEGMENT_LENGTH-1;
        } else {
            towrite = cur_length - total_written;
        }
        memcpy(fastaSegment, seqF+total_written, towrite * sizeof(char));
        fastaSegment[towrite] = '\0';
        GZIP_PRINTF(out, "%s\n", fastaSegment);
        total_written += towrite;
    }
    assert(total_written == cur_length);
    return total_written;
}


#include "hash_funcs.h"
//#define hashkey(table_size, key, len) (MurmurHash3_x64_64(key, len) % (table_size))
//#define hashstr(table_size, key) (MurmurHash3_x64_64(key, strlen(key)) % (table_size))
#define HipmerHash(key, len) (MurmurHash3_x64_64(key, len))
#define hashkey(table_size, key, len) (HipmerHash(key, len) % (table_size))
#define hashstr(table_size, key) (HipmerHash(key, strlen(key)) % (table_size))

static inline void print_args(option_t *optList, const char *stage)
{
    if (!MYTHREAD) {
        option_t *olist = optList;
        option_t *opt = NULL;
        //serial_printf(KLWHITE "STAGE %s ", stage);
        serial_printf("STAGE %s ", stage);
        while (olist) {
            opt = olist;
            olist = olist->next;
            serial_printf("-%c ", opt->option);
            if (opt->argument)
                serial_printf("%s ", opt->argument);
        }
        //serial_printf(KNORM "\n");
        serial_printf("\n");
    }
#ifdef __UPC_VERSION__
    upc_barrier;
#endif
}

static inline int get_shared_heap_mb(void) 
{
    char *shared_heap_size = getenv("UPC_SHARED_HEAP_SIZE");
    if (!shared_heap_size)
        return 0;
    int len = strlen(shared_heap_size);
    if (!len)
        return 0;
    int64_t shared_heap_bytes = atoi(shared_heap_size);
    if (shared_heap_size[len - 1] == 'G') 
        shared_heap_bytes *= ONE_GB;
    else if (shared_heap_size[len - 1] == 'M') 
        shared_heap_bytes *= ONE_MB;
    else if (shared_heap_size[len - 1] == 'K') 
        shared_heap_bytes *= ONE_KB;
    else 
        shared_heap_bytes *= ONE_MB;
    return shared_heap_bytes / ONE_MB;
}

// for printing out diagnostics
static void init_diags(void)
{
    if (!_sv) return;
    if (!MYTHREAD) 
        MYSV._my_diags = fopen_chk("diags.log", "a");
}

static void fini_diags(void)
{
    if (!_sv) return;
    if (MYSV._my_diags) {
        fclose(MYSV._my_diags);
        MYSV._my_diags = NULL;
    }
}

#define ADD_DIAG(type, key, val)                                        \
    do {                                                                \
        if (_sv != NULL) {                                              \
            if (MYSV._my_diags != NULL) {                      \
                char *fname = strdup(__FILENAME__);                     \
                char *s = strchr(fname, '.');                           \
                if (s)                                                  \
                    s[0] = '\0';                                        \
                fprintf(MYSV._my_diags, "%s.%s\t" type "\n", fname, key, val); \
                fflush(MYSV._my_diags);                         \
                free(fname);                                            \
            }                                                           \
        }                                                               \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif // __COMMON_H

