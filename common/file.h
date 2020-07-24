#ifndef __FILE_H
#define __FILE_H

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
#include <zlib.h>

#include "../common/colors.h"
#include "../common/common_rankpath.h"
#include "../common/optlist.h"
#include "../common/StaticVars.h"
// INCLUDE other common header files too

//#include "version.h"
#include "log.h"

/****************************************************************************************/

#define fclose_track(f) fclose_track_at_line(f, __FILENAME__, __LINE__)
static inline int fclose_track_at_line(FILE *f, const char *which_file, int which_line)
{
#ifdef TRACK_FILE_OPENS
    fflush(stderr);
    LOG("\n%p/%d CLOSE. %s:%d\n", f, MYTHREAD, which_file, which_line);
    fflush(stderr);
#endif
    int64_t size = ftell(f);
    int ret = fclose(f);
    int isMyLog = _sv == NULL ? 0 : 1;
    if (isMyLog) isMyLog = f == MYSV._my_log;
    if (isMyLog) MYSV._my_log = NULL; // check to ensure _my_log is never used again even in this function
    if (ret != 0) WARN("Could not close file: %p! %s:%d\n", f, which_file, which_line);
    if (!isMyLog) LOGF("fclose_track(%p) finished. %lld bytes. %s:%d\n", f, (long long int) size, which_file, which_line);
    return ret;
}

#define fopen_track(p, m) fopen_track_at_line(p, m, __FILENAME__, __LINE__)
static FILE *fopen_track_at_line(const char *path, const char *mode, const char *which_file, int which_line)
{
    FILE *f = fopen(path, mode);
    if (!f) WARN("Could not open %s with '%s' mode: %s. (%s:%d)\n", path, mode, strerror(errno), which_file, which_line);
#ifdef TRACK_FILE_OPENS
    fflush(stderr);
    if (_sv != NULL && MYSV._my_log) LOGF("\n%p/%d OPEN %s (%s:%d)\n", f, MYTHREAD, path, which_file, which_line);
    fflush(stderr);
#else
    if (_sv != NULL && MYSV._my_log) LOGF("Opening %s with mode '%s': %s  (%s:%d)\n", path, mode, f == NULL ? "FAILED!" : "success.", which_file, which_line);
#endif
    return f;
}

#define fopen_chk(p, m) fopen_chk_at_line(p, m, __FILENAME__, __LINE__)
static inline FILE *fopen_chk_at_line(const char *path, const char *mode, const char *which_file, int which_line)
{
    FILE *f = fopen_track_at_line(path, mode, which_file, which_line);
    if (!f) {
        fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: Could not open file '%s': %s\n" KNORM, 
                MYTHREAD, which_file, which_line, NULL, path, strerror(errno));
        EXIT_FUNC(1);
    }
    return f;
}

#ifdef NO_GZIP

#define YES_GZIP 0
typedef FILE * GZIP_FILE;
#define GZIP_GETS(b,l,f) fgets(b,l,f)
#define GZIP_OPEN fopen_track
#define GZIP_OPEN_RANK_PATH fopen_rank_path
#define GZIP_CLOSE fclose_track
#define GZIP_REWIND rewind
#define GZIP_FWRITE fwrite
#define GZIP_FREAD fread
#define GZIP_FTELL ftell
#define GZIP_EOF feof
#define GZIP_PRINTF fprintf
#define GZIP_EXT ""

#else

#define YES_GZIP 1
typedef gzFile GZIP_FILE;
#define GZIP_GETS(b,l,f) gzgets(f,b,l)
#define GZIP_OPEN gzopen_track
#define GZIP_OPEN_RANK_PATH gzopen_rank_path
#define GZIP_CLOSE gzclose_track
#define GZIP_REWIND gzrewind
#define GZIP_FWRITE(ptr, size, n, f) gzwrite(f, ptr, ((long long int)size)*((long long int)n))
#define GZIP_FREAD(ptr, size, n, f) gzread(f, ptr, ((long long int)size)*((long long int)n))
#define GZIP_FTELL gztell
#define GZIP_EOF gzeof
#define GZIP_PRINTF /*** WARNING gzprintf is limited to 8191 bytes! ***/ gzprintf
#define GZIP_EXT ".gz"

#endif

#define gzopen_track(p,m) gzopen_track_at_line(p, m, __FILENAME__, __LINE__)
static inline gzFile gzopen_track_at_line(const char *path, const char *_mode, const char * which_file, int which_line)
{
    char mode[5];
    // choose fastest compression mode if not chosen already
    if (YES_GZIP && ( strcmp(_mode, "w") == 0 || strcmp(_mode, "wb") == 0 || strcmp(_mode, "a") == 0 || strcmp(_mode, "ab") == 0) ) {
        sprintf(mode, "%s1", _mode);
    } else {
        strcpy(mode, _mode);
    }
    LOGF("gzopening %s with %s mode (%s-%d)\n", path, mode, which_file, which_line); 
    gzFile f = gzopen(path, mode); 
#if ZLIB_VERNUM >= 0x1240
    gzbuffer(f, 128*1024);
#endif
    if (!f) WARN("Could not open %s with %s mode! %s\n", path, mode, strerror(errno)); 
    return f;
}
#define gzopen_chk(p, m) gzopen_chk_at_line(p,m,__FILENAME__, __LINE__)
static inline gzFile gzopen_chk_at_line(const char *path, const char *mode, const char *which_file, int which_line)
{
    gzFile f = gzopen_track(path, mode);
    if (!f) DIE("Could not open gzfile: %s\n", path);
    return f;
}
#define gzclose_track(gz) gzclose_track_at_line(gz, __FILENAME__, __LINE__)
static inline int gzclose_track_at_line(gzFile gz, const char * which_file, int which_line)
{
    int64_t size = gztell(gz);
    int status = gzclose(gz);
    if (Z_OK != status) WARN("Could not close a gzFile! %s:%d\n", which_file, which_line);
    LOGF("Closed gzip file. %lld bytes\n", (long long int) size);
    return status == Z_OK ? 0 : 1;
}
#define gzclose_chk(gz) gzclose_chk_at_line(gz, __FILENAME__, __LINE__)
static inline void gzclose_chk_at_line(gzFile gz, const char * which_file, int which_line)
{
    if (! gzclose_track_at_line( gz, which_file, which_line ) )
      DIE("Could not close gzfile! %s:%d\n", which_file, which_line);
}

#define fopen_rank_path(p, m, r) fopen_rank_path_at_line(p, m, r, __FILENAME__, __LINE__)
static inline FILE * fopen_rank_path_at_line(char *buf, const char *mode, int rank, const char *which_file, int which_line) {
    char *rankpath = get_rank_path(buf, rank);
    FILE *f = fopen_chk_at_line(rankpath, mode, which_file, which_line);
    return f;
}

#define gzopen_rank_path(p, m, r) gzopen_rank_path_at_line(p, m, r, __FILENAME__, __LINE__)
static inline gzFile gzopen_rank_path_at_line(char *buf, const char *mode, int rank, const char *which_file, int which_line) {
    char *rankpath = get_rank_path(buf, rank);
    gzFile f = gzopen_chk_at_line(rankpath, mode, which_file, which_line);
    return f;
}

#define link_chk(o,n) link_chk_at_line(o, n, __FILENAME__, __LINE__)
static inline void link_chk_at_line(const char *old_file, const char *new_file, const char *which_file, int which_line) {
    int ret = link(old_file, new_file);
    if (ret != 0) {
        if (errno == EEXIST) {
          ret = unlink(new_file);
          if (ret == 0) return link_chk(old_file, new_file); // try again!
        } else if (errno == EPERM) {
          ret = symlink(old_file, new_file); // try a symlink!
          if (ret != 0) WARN("Could not link %s to %s! %s\n", old_file, new_file, strerror(errno));
          return; // well at least we tried!
        }
    }
    if (ret != 0)
        DIE("Could not link %s to %s! [%s:%d] %s\n", old_file, new_file, which_file, which_line, strerror(errno));
}

#endif // __FILE_H

