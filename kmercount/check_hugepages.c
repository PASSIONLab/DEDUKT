#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#ifndef USE_MPI
#define USE_MPI 1
#endif

#if USE_MPI
#include <mpi.h>
#endif

#if (  _SVID_SOURCE || _BSD_SOURCE || _POSIX_C_SOURCE >= 1 || _XOPEN_SOURCE || _POSIX_SOURCE )
#else
char *strtok_r(char *str, const char *delim, char **saveptr);
#endif

#if ( _SVID_SOURCE || _BSD_SOURCE || _XOPEN_SOURCE >= 500 || _XOPEN_SOURCE && _XOPEN_SOURCE_EXTENDED || _POSIX_C_SOURCE >= 200809L )
#else
char *strdup(const char *s);
#endif

#if ( _BSD_SOURCE || _XOPEN_SOURCE >= 500 || _POSIX_C_SOURCE >= 200112L )
#else
int gethostname(char *name, size_t len);
#endif

#define ONE_MB (1024L*1024L)

long parseMemFree(char * line, long pagesize, long minSize) {
	char *str1 = line, *savePtr = NULL;
	char *token = NULL;
	long curSize = pagesize;
	long memFree = 0;

	int tok_i = 0;
	while (token = strtok_r(str1, " \t", &savePtr)) {
		str1 = NULL;
		tok_i++;
		if (tok_i <= 4) continue;
		long num = atoi(token);
		if (curSize >= minSize) {
			memFree += curSize * num;
		}
		//printf("Found %d from %s at tok_i=%d curSize=%ld memFree=%ld MB\n", num, token, tok_i, curSize, memFree / ONE_MB);
		curSize *= 2;
	}
	
	return memFree;
}

const char *USAGE = "check_hugepages minPageSize(=2097152)  [ 'buddinfo prefix' [ threads testAlloc ... ] ]\n\tcat /proc/buddyinfo for the pattern you desire\n";
int main(int argc, char **argv) {

	// awk 'BEGIN{sum=0} /Node 0,/ {for(i=14; i<=NF; i++) sum+=$i * 2^(i-14) * 2} END{ print ENVIRON["SLURMD_NODENAME"]"\t"sum; } '
	/*
	Node 0, zone      DMA      0      0      0      0      2      1      1      1      0      1      3
	Node 0, zone    DMA32    185    151     81     47     74     52     48     35     13      8    232
	Node 0, zone   Normal 175089 996929 868062 581203 220417 158572  35187  28391  20019      0      0
	Node 1, zone   Normal 286738 614177 343565 144680    443      0      0      0      0      0      0
	*/

	int size = 1, rank = 0, i;
#if USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 2) {
		if (rank == 0) fprintf(stderr, "%s", USAGE);
		return 1;
	}

	char host[64];
	gethostname(host, 63);
	int maxOutSize = 2048;
	char *output = (char*) calloc(maxOutSize,1);;

	char *match = "";
	long minSize = 2*ONE_MB;
	if (argc >= 2)
		minSize = atol(argv[1]);
	if (argc >= 3)
		match = argv[2];
	int len = strlen(match);

	long pagesize = sysconf(_SC_PAGESIZE);

	FILE *f = fopen("/proc/buddyinfo", "r");
	char line[256];
	long freeHugePages = 0;
	while ( fgets(line, 255, f) ) {
		if (strstr(line, "Normal") == NULL) continue;
		if (len == 0 || strncmp(line, match, len) == 0) {
			char *parse = strdup(line);
			long memFree = parseMemFree(parse, pagesize, minSize);
			sprintf(output + strlen(output), "%ld MB on %dof%d %s %.*s\n", memFree / ONE_MB, rank, size, host, (len>0?len:22), line);
			freeHugePages += memFree;
			free(parse);
		}
	}
	fclose(f);
	sprintf(output + strlen(output), "%ld MB on %dof%d %s Total\n", freeHugePages / ONE_MB, rank, size, host);
	
#if USE_MPI
	char *recv = NULL;
	if (rank == 0) recv = malloc(size*maxOutSize);
	
	MPI_Gather(output, maxOutSize, MPI_CHAR, recv, maxOutSize, MPI_CHAR, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		for(i = 0; i < size; i++) {
			printf("%s", recv + i * maxOutSize);
		}
	}

	if (recv) free(recv);

#else
	printf("%s", output);
#endif

	if (argc >= 5) {
		int threads = atoi(argv[3]);
		double minTime = 100000.0, maxTime = 0.0, totalTime = 0.0;
		int64_t totMB = 0;
		# pragma omp parallel num_threads(threads) reduction(min:minTime) reduction(max:maxTime) reduction(+:totalTime) reduction(+:totMB)
		{
			int64_t **tests = (int64_t**) calloc(argc, sizeof(int64_t*)), j;
			int i;
			struct timespec start, end;
			clock_gettime(CLOCK_MONOTONIC, &start);
			for (i = 4; i < argc; i++) {
				int64_t mb = atoi(argv[i]);
				totMB += mb;
				int64_t *x = tests[i] = calloc(mb*1024L*1024L/sizeof(int64_t), sizeof(int64_t));
				for (j=0; j < mb*1024L*1024L/sizeof(int64_t); j++) x[j] = j+1;	
				if (x == NULL) fprintf(stderr, "%s thread %d of %d failed to calloc %ld mb\n", host, omp_get_thread_num(), omp_get_num_threads(), mb);
			}
			clock_gettime(CLOCK_MONOTONIC, &end);
			# pragma omp barrier
			for (i = 4; i < argc; i++) {
				if (tests[i]) free(tests[i]);
			}
			free(tests);
			minTime = maxTime = totalTime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1000000000.0;
		}
		
		output[0] = '\0';
		sprintf(output, "%0.6f %0.6f %0.6f on %dof%d %s Alloc %ld MB over %d threads min/max/avg\n", minTime, maxTime, totalTime / threads, rank, size, host, totMB, threads);
#if USE_MPI
		char *recv = NULL;
		int outSize = 128;
		if (rank == 0) recv = malloc(size * outSize);
		MPI_Gather(output,outSize, MPI_CHAR, recv, outSize, MPI_CHAR, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			for(i = 0; i < size; i++) {
				printf("%s", recv + i * outSize);
			}
		}
		if (recv) free(recv);
#else
		printf("%s", output); 
#endif
	}

	free(output);

		
#if USE_MPI
	MPI_Finalize();
#endif

	return 0;
}
