#ifndef __EXEC_UTILS_H
#define __EXEC_UTILS_H

#include "config.h"

#ifndef DBG
#if 1
#define DBG(...)                                \
    do {                                        \
        fprintf(stderr, __VA_ARGS__);           \
        fflush(stderr);                         \
    } while (0)
#else
#define DBG(...)
#endif
#endif

#define HASH_BAR "########################################################################"

int exec_stage(char *stages, int is_scaff, int (*stage_main)(int, char**), const char *stage_name, ...);
void print_timings(cfg_t *cfg); // TODO: implementation currently commented out
int64_t get_num_from_file(const char *prepend, const char *fname);
//int64_t get_num_from_file_and_broadcast(const char *prepend, const char *fname);
char *get_str_from_file(const char *fname);
//int64_t get_long_from_file_and_broadcast(const char *fname);
void put_str_in_file(const char *fname, char *s);
int dummy_exec(char *stage_name);

#endif
