#ifndef __MACH_H
#define __MACH_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdint.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>

#include "log.h"
#include "file.h"
#include "defines.h"

typedef struct {
    size_t vm_size;
    size_t vm_rss;
    size_t share;
} memory_usage_t;

// prototypes
static int _get_cores_per_node_HW(void) ;
static void get_mem_usage(memory_usage_t *mem);
static int get_max_mem_usage_mb(void);
static double get_free_mem_gb(void);
static double get_used_mem_gb(void);

#if defined(__APPLE__) && defined(__MACH__)

/*
// use mach_info
#include <mach/vm_statistics.h>
#include <mach/task.h>
#include <mach/task_info.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
*/
#include <sys/types.h>
#include <sys/sysctl.h>

static int _get_cores_per_node_HW(void) {
    int count = 1;
    size_t count_len = sizeof(count);
    sysctlbyname("hw.physicalcpu", &count, &count_len, NULL, 0);
    return count;
}

static void get_mem_usage(memory_usage_t *mem) {
/*
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS != task_info(mach_task_self(),
                              TASK_BASIC_INFO, (task_info_t)&t_info, 
                              &t_info_count))
    {
        mem->vm_size = 0;
        mem->vm_rss = 0; 
    } else {
        // resident size is in t_info.resident_size;
        // virtual size is in t_info.virtual_size;
        mem->vm_size = t_info.virtual_size;
        mem->vm_rss *= t_info.resident_size;
    }
*/
    mem->vm_size = 0;
    mem->vm_rss = 0;
    mem->share = 0;

}
static int get_max_mem_usage_mb(void) {
    return 0;
}

static double get_free_mem_gb(void) {
    return 0.0;
}

static double get_used_mem_gb(void) {
    return 0.0;
}

#else // NOT __APPLE__

// use the /proc filesystem

static int _get_cores_per_node_HW(void)
{
    char buf[256];
    FILE *f = fopen_chk("/proc/cpuinfo", "r");
    int cpu_cores = -1;
    int physical_id = -1;
    while (fgets(buf, 256, f)) {
        if (strstr(buf, "cpu cores")) {
            cpu_cores = atoi(strstr(buf, ": ") + 1);
        } else if (strstr(buf, "physical id")) {
            physical_id = atoi(strstr(buf, ":") + 1);
        }
    }
    fclose(f);
    return cpu_cores * (physical_id + 1);
}

static void get_mem_usage(memory_usage_t *mem)
{
    FILE *f = fopen_chk("/proc/self/statm", "r");
    if (fscanf(f, "%lu %lu %lu", &(mem->vm_size), &(mem->vm_rss), &(mem->share)) != 3) {
        WARN("could not read memory usage\n");
        memset(mem, 0, sizeof(memory_usage_t));
        return;
    }
    fclose_track(f);
    int ps = getpagesize();
    mem->vm_size *= ps;
    mem->vm_rss *= ps;
    mem->share *= ps;
}

static int get_max_mem_usage_mb(void)
{
    int usage;
    FILE *f = fopen_chk("/proc/self/status", "r");
    char buf[1000];
    char units[5];
    while (fgets(buf, 999, f)) {
        if (strncmp(buf, "VmHWM:", 6) == 0) {
            sscanf(buf + 6, "%d %s\n", &usage, units);
            break;
        }
    }
    if (strcmp(units, "kB") == 0)
        usage /= 1024;
    fclose_track(f);
    return usage;
}


#define SCANFIELD(s, var)                                               \
    if (strncmp((buf), (s), strlen(s)) == 0) {                          \
        sscanf(buf + strlen(s), "%ld %s\n", &var, units);               \
    }
        

static double get_free_mem_gb(void)
{
    // flush file buffers to disk so memory usage is more accurate???
//    sync();
    char buf[1000];
    char units[5];
    long memfree = 0, buffers = 0, cached = 0, swapcached = 0, active = 0, inactive = 0,
        dirty = 0, writeback = 0, shmem = 0, slab = 0, kernelstack = 0, pagetables = 0, bounce = 0;
    FILE *f = fopen_chk("/proc/meminfo", "r");
    while (fgets(buf, 999, f)) {
        SCANFIELD("MemFree:", memfree);
        SCANFIELD("Buffers:", buffers);
        SCANFIELD("Cached:", cached);
        SCANFIELD("SwapCached:", swapcached);
        SCANFIELD("Active:", active);
        SCANFIELD("Inactive:", inactive);
        SCANFIELD("Dirty:", dirty);
        SCANFIELD("Writeback:", writeback);
        SCANFIELD("Shmem:", shmem);
        SCANFIELD("Slab:", slab);
        SCANFIELD("KernelStack:", kernelstack);
        SCANFIELD("PageTables:", pagetables);
        SCANFIELD("Bounce:", bounce);
    }
    fclose_track(f);
    /*
    memory_usage_t mem;
    get_mem_usage(&mem);
    serial_printf("%c: MemFree %ld Buffers %ld Cached %ld Swapcached %ld Active %ld Inactive %ld "
                  "Shmem %ld PageTables %ld VmRSS %ld\n", units[0],
                  memfree, buffers, cached, swapcached, active, inactive, shmem, pagetables,
                  mem.vm_rss);
    */
    memfree += buffers + cached + swapcached;
    if (units[0] == 'k') 
        return (double)memfree / ONE_MB;
    else
        return memfree;
}


static double get_used_mem_gb(void)
{
    /*
     ok, it seems the key here is to read /proc/<pid>/smaps and look at the Referenced field to
     find out the exact memory usage. Probably the only important entries are stack and heap.
     From man procfs, the section on /proc/[pid]/clear_refs:
 
        Clearing the PG_Referenced and ACCESSED/YOUNG bits provides a method to measure
        approximately how much memory a process is using.  One first inspects the values in the
        "Referenced" fields for the VMAs shown in /proc/[pid]/smaps to get an idea of the memory
        footprint of the process.  One then clears the PG_Referenced and ACCESSED/YOUNG bits and,
        after some measured time interval, once again inspects the values in the "Referenced" fields
        to get an idea of the change in memory footprint of the process during the measured
        interval.  If one is interested only in inspecting the selected mapping types, then the
        value 2 or 3 can be used instead of 1.

     */
    char buf[1000];
    char units[5];
    FILE *smaps_fd = fopen("/proc/self/smaps", "r");
    if (!smaps_fd) {
        WARN("Could not open /proc/self/smaps for checking memory usage\n");
        return 0;
    } else {
        double mem_used_gb = 0;
        int to_add = 0;
        unsigned int addr_start, addr_end;
        char perms[100], offset[100], dev[100], inode[100], pathname[1024];
        while (fgets(buf, 999, smaps_fd)) {
            int n = sscanf(buf, "%x-%x %s %s %s %s %s", &addr_start, &addr_end, perms, offset,
                           dev, inode, pathname);
            if (n == 7) {
                if (pathname[0] == '[') {
                    //serial_printf("found file mapping: %s", buf);
                    to_add = 1;
                }
            } else if (n == 6) {
                //serial_printf("found anon mapping: %s", buf);
                to_add = 1;
            }
//            if (strstr(buf, "Referenced") && to_add) {
            if (strstr(buf, "Pss:") && (buf[0] == 'P') && to_add) {
                //serial_printf("%s", buf);
                double mem_gb = 0;
                long mem_raw;
//                sscanf(buf + strlen("Referenced:"), "%ld %s\n", &mem_raw, units);
                sscanf(buf + strlen("Pss:"), "%ld %s\n", &mem_raw, units);
                if (units[0] == 'k')
                    mem_gb = (double)mem_raw / ONE_MB;
                else if (units[0] == 'm')
                    mem_gb = (double)mem_raw / ONE_KB;
                else
                    mem_gb = mem_raw;
                mem_used_gb += mem_gb;
                to_add = 0;
            }
        }
        fclose(smaps_fd);
        return mem_used_gb;
    }
}

#endif

static int get_cores_per_node(void) {
    int count = 0;
#if defined(__UPC__) && defined(_UPC_COMMON_H)
    count = divineThreadsPerNode();
    DBG("Found %d threads per node through UPC\n", count);
#ifdef DEBUG
    int upcCount = count;
    count = 0; // run the hardware routine and validate!
#endif
#endif
    if (count == 0) {
        count = _get_cores_per_node_HW();
    }
#if defined(__UPC__) && defined(DEBUG)
    if (upcCount > 0 && upcCount != count) {
        WARN("divineThreadsPerNode() == %d whereas _get_cores_per_node_HW looks discovered %d physical cores... What's up with that?\n", upcCount, count);
        count = upcCount; // because DEBUG ran extra code an Release would not have used _get_cores_per_node_HW()
    }
#endif
    LOGF("Found %d cores_per_node\n", count);
    return count;
}

#ifdef SHOW_MEMORY
static void print_memory_used(void) 
{
    if (!MYTHREAD) {
        memory_usage_t mem;
        get_mem_usage(&mem);
        serial_printf("Memory used (GB): VmSize %.2f, VmRSS %.2f, share %.2f\n",
                      (double)mem.vm_size / ONE_GB, (double)mem.vm_rss / ONE_GB, 
                      (double)mem.share / ONE_GB);
    }
}
#else
#define print_memory_used() /* noop */
#endif 

#endif // __MACH_H
