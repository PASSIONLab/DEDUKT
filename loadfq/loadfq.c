#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <string.h>
#include <libgen.h>
#include <sys/types.h>
#include <dirent.h>

#include "../version.h"

#include "loadfq.h"

#ifndef MAX_LIBRARIES
#define MAX_LIBRARIES 256 //128
#endif

static char library_names[MAX_LIBRARIES][LIB_NAME_LEN];
static int num_lib_names = 0;

static void get_library_names(char *s) 
{
    num_lib_names = 0;
    int len = strlen(s);
    if (len >= MAX_LIBRARIES*4) 
        DIE("Too long a library string\n"); 
    int start = 0;
    int name_i = 0;
    for (int i = 0; i < len; i++) {
        if (s[i] == ',') {
            int l = i - start;
            strncpy(library_names[name_i], s + start, l);
            library_names[name_i][l] = 0;
            num_lib_names++;
            serial_printf("Using lib %d: %s\n", name_i, library_names[name_i]);
            name_i++;
            if (name_i >= MAX_LIBRARIES) 
                DIE("Too many libraries!\n"); 
            start = i + 1;
        }
    }
    strcpy(library_names[name_i], s + start);
    serial_printf("Using lib %d: %s\n", name_i, library_names[name_i]);
    num_lib_names++;
}

int loadfq_main(int argc, char** argv)
{
    option_t *this_opt;
    option_t *opt_list = GetOptList(argc, argv, "D:f:N:B:");
    print_args(opt_list, __func__);
    char *libnames_ext = NULL;
    const char *base_dir = "/dev/shm";
    int coresPerNode = 1; //(MYSV.cores_per_node > 0)? MYSV.cores_per_node: 1; // restoring hipmer version
    int do_purge = 1;
    while (opt_list) {
        this_opt = opt_list;
        opt_list = opt_list->next;
        switch (this_opt->option) {
        case 'f':
            libnames_ext = this_opt->argument;
            serial_printf("string following -f is %s\n", this_opt->argument);
            serial_printf("libnames_ext set to %s\n", libnames_ext);
            break;
        case 'N':
            coresPerNode = atoi(this_opt->argument);
            SET_CORES_PER_NODE(coresPerNode);
            break;
        case 'B':
            base_dir = this_opt->argument;
            assert(base_dir != NULL);
            break;
        case 'D':
        	do_purge=0;
        	serial_printf("Unset do purge: will not purge /dev/shm before loading files.\n");
        	break;
        default:
            SDIE("Invalid option %c\n", this_opt->option);
      }
    }
    if (!libnames_ext) {
        SDIE("Usage: %s -f fofn_file_name\n", argv[0]);
    }
    serial_printf("libnames_ext=%s\n", libnames_ext);
    double t = MPI_Wtime();

    // delete any files from /dev/shm to ensure we have adequate space
    if (do_purge && MYTHREAD % coresPerNode == 0) {
        lld bytes = 0;
        int files = 0, dirs = 0;
        int purged = purgeDir("/dev/shm/per_thread", &bytes, &files, &dirs);
        if (purged) {
            LOG("purged %d files from /dev/shm/per_thread in %0.2f s (%0.3lf GB)\n",
                purged, (MPI_Wtime() - t), ((double) bytes) / ONE_GB);
            if (purged != files) WARN("Could not purge %d files in /dev/shm/per_thread!\n", files - purged);
        }
    }
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    // now create the /dev/shm directories
    char my_dir[MAX_FILE_PATH];
    strcpy(my_dir, "/dev/shm/.");
    if (MYTHREAD % coresPerNode == 0) 
        get_rank_path(my_dir, -1);
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    strcpy(my_dir, "/dev/shm/.");
    get_rank_path(my_dir, MYTHREAD);

    serial_printf("Loading FASTQ files from %s\n", libnames_ext);
    get_library_names(libnames_ext);

    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );

    for (int libi = 0; libi < num_lib_names; libi++) {
        double start = MPI_Wtime();
        char fofn[MAX_FILE_PATH], fname[MAX_FILE_PATH];
        snprintf(fofn, MAX_FILE_PATH, "%s", library_names[libi]);
        serial_printf("Processing fastq files from %s...\n", fofn);
        FILE *fofn_fd = fopen_chk(fofn, "r");
        while (fgets(fname, MAX_FILE_PATH, fofn_fd)) { 
            double start1 = MPI_Wtime();
            fname[strlen(fname)-1] = '\0';
            serial_printf("Loading %s into memory on %d threads\n", fname, THREADS);
            fq_reader_t fqr = create_fq_reader();
            int err = load_fq(fqr, fname, base_dir);
            if (err < 0)
                DIE("Could not load %s into memory: %s\n", fname, strerror(-err));
            double end1 = MPI_Wtime();
            double elapsed1 = end1-start1;
            double MB = 1.0/ONE_MB*(fqr->end_read - fqr->start_read);
            LOGF("Loading times for %s: read %0.3f MB in %0.3f s wrote in %0.3f s (%0.3f MB/s)\n", fname, MB, fqr->readTime, fqr->writeTime, (elapsed1 != 0 ? MB / elapsed1 : 0.0) );
            destroy_fq_reader(fqr);
        }
        if (fclose_track(fofn_fd) != 0) WARN("Could not close %s! %s\n", fofn, strerror(errno));
        LOGF("Loaded lib %d %s in %0.3f s\n", libi, fofn, MPI_Wtime() - start);
    }
    CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
    serial_printf("Overall time for %s is %.2f s\n", basename(argv[0]), (MPI_Wtime() - t) );
    return 0;
}


#ifndef SINGLE_EXEC

StaticVars _sv = NULL;
int main(int argc, char **argv) 
{
	CHECK_MPI( MPI_Init(&argc, &argv) );
    OPEN_MY_LOG("loadfq");
    serial_printf("Starting diBELLA:loadfq version %s on %d threads\n", DIBELLA_VERSION, THREADS);
    int ret = loadfq_main(argc, argv);
    MPI_Finalize();
    return ret;
}

#endif
