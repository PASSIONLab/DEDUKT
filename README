 ## Requirements 

1. Message Passing Interface (MPI) environment
        (a) Open MPI, https://www.open-mpi.org/
        (b) MPICH2, http://www.mcs.anl.gov/project/mpich-high-performance-portable-implementation-mpi

2. A working C/C++ compiler with OpenMP
       (a) e.g. GNU 6.2.0, 7.2.0, 7.4.0, 8.2.0

3. CMake >= 3.5 (C++ 11 support required)

## Building and Installing
  Overview
       The build process is automated with bash scripts and cmake. 
       Build environment scripts are provided (and automatically loaded) for several platforms.
       Customization is also available via the command line or by providing a custom environment script.
       
       For general-purpose systems, little if any environmental configuration is necessary.
       You may want to specify the C/C++ compiler(s) and the install directory.
       
       For more complex systems (e.g. HPC clusters), look at the provided 
       environment scripts, and copy-modify as needed.
       
  Default Build and Installation
       1. To build and install automatically, from the source directory run:
          $  sh install.sh
       
          Or to build and install the Debug version, from the source directory run:
          $  DEBUG=1 sh install.sh

       Some variables you may want to redefine are:
         BUILD_DIR     :    the path to your desired build directory
         INSTALL_DIR   :    the path to your desired build directory
         CC            :    the C compiler
         CXX           :    the C++ compiler 
       

       Additionally, to customize particular feature settings, add them to the cmake command as "-D<feature name1> -D<feature name 2> " and so on,
         or if using the install.sh script, as DEFS="-D<feature 1 name> -D<feature 2 name>" sh install.sh.
         
       Some (advanced) program-specific cmake variables of interest:
           TIGHT_READS ; Enable this optimization to reduce memory usage ONLY IF all reads in the input set are 
                         (strictly) shorter than 65,536 bps (OR) if k-mer/seed positions are not relevant for 
                         your use-case. The results will be incorrect if any of the input reads are longer than 
                         65,5365bps! This optimization is NOT enabled by default.
                         Example usage: DEFS="-D TIGHT_READS=1" sh install.sh
                                    OR:  cmake -D TIGHT_READS=1 /path/to/src

           MAX_NUM_READS_LIST ; Setting MAX_NUM_READS_LIST to a semi-colon separated list of integers will build program 
                                binaries with statically-sized arrays equal to the sizes specified in this list. Including
                                the reliable-kmer-upper-limit (see the -u runtime flag below) with which you intend to run
                                the program will save memory for that run. Including powers of 2 larger than your "-u" value
                                will pad the arrays, (potentially but not necessarily) improving memory access times.
                                At minimum, this flag should be set if the default values are less than your intended "-u" value.

                                Example usage:  DEFS="-DMAX_NUM_READS_LIST=4;8;16" sh install.sh
                                This will install dibella with 3 binaries per executable, corresponding to the 3 values specified.

           HEAVYHITTERS       ; Enable this optimization for data sets with an expected, relatively small, set of k-mers that heavily 
                                dominate the k-mer distribution. Also try if the reported load imbalance is significantly higher than 1.
                                Example usage:  DEFS="-DHEAVYHITTERS=1" sh install.sh
                                           OR:  cmake -D HEAVYHITTERS=1 /path/to/src

           BENCHMARKONLY      ; Set to 1 to measure computation and communication time with minimal intermediate output of any kind.
                                Example usage:  DEFS="-DBENCHMARKONLY=1" sh install.sh
                                           OR:  cmake -D BENCHMARKONLY=1 /path/to/src


  
  ## Binaries 
    Executable binaries - kmermatch-<u_max>-<k_max> and scripts can be found in $INSTALL_DIR.
    There are three primary executables: .
    Each of <u_max> and <k_max> are substituted for various integers in the actual installation.
    The former, <u_max>, corresponds to the maximumm reliable k-mer upper limit. 
    The pipeline (or pipeline stage) can be run with any -u value less than or equal to <u_max> with the respective binary.
    Similarly, <k_max>, corresponds to the maximum value of k (the next higher power of 2 e.g. 32 for k=17) with which the binary can be run.
    
    To compile the code base with values of <u_max> other than the default, specify a list of semi-colon separated values using the cmake flag "-DMAX_NUM_READS_LIST".
    You can do so on the command line (example below) or in your environment script (see the "Custom Build and Installation" section).
    
    $  DEFS="-DMAX_NUM_READS_LIST=8;16" sh install.sh 
    

Running diBELLA with SLURM or PBS
  See the provided run scripts under <path to source>/scripts.
  
  To run an individual job with SLURM, copy-modify and run scripts/sbatch_single.sh.
  
  To run a scaling job with SLURM, copy-modify and run scripts/sbatch_scaling.sh.
  
  For PBS, use qsub_single.sh and qsub_scaling.sh.

Running diBELLA Directly
  The pipeline can be run from the "main" executable with a MPI task executor (e.g. mpirun, aprun, srun,...).
  See the "Binaries" section above for an explanation of the exact names of executables (e.g. main-8-32).
  
  Example 1, running the full pipeline:  
  $  mpirun -n 4 main-8-32 -k 17 -u 8 -i list_of_fastqs.txt
  This will execute the entire pipeline using 4 threads (MPI tasks), a k-mer length of 17, and a reliable k-mer upper limit of 8.
  "list_of_fastqs.txt" is any text file containing a line-separated list of input fastq files.

  Example 2, running just the first stage:  
  $  mpirun -n 4 main-8-32 -k 17 -u 8 -i list_of_fastqs.txt -s kmermatch-17
  This will execute just the first stage of the pipeline, "kmermatch";
  "kmermatch" counts k-mers, prunes them (down to the reliable set), and uses them to compute read overlaps.
  The output files (by default for performance) will be written to separate per_thread directories.
  When running the entire pipeline from main, they will be located between stages automatically.
  To concatenate these files, for some other purpose, into one large file in the run directory (RUN_DIR) e.g.:
  $  cat $RUN_DIR/per_thread/*/*/overlaps-17_* > $RUN_DIR/overlaps-17
  
  Example 3, running just the alignment stage:
  $  mpirun -n 4 main-8-32 -k 17 -u 8 -i list_of_fastqs.txt -m overlaps-17 -o alignments-17 -s alignment-17
  This will execute just the alignment stage. If an input overlaps file (specified as "overlaps-17" in the example)
  is not provided (or does not contain overlap information generated by kmermatch), running this stage will fail. 
  The output files (by default for performance) will be written to separate per_thread directories.
  To concatenate these files into one large file in the run directory (RUN_DIR) e.g.:
  $  cat $RUN_DIR/per_thread/*/*/alignments-17_* > $RUN_DIR/alignments-17
  When using the provided run scripts, this is done automatically.

  Each stage can also be executed individually (without running "main"). 
  Additional options (not shown) are available for each stage run either way.
  Use the -h flag with a "main", "kmermatch", or "alignment" a binary to view the respective usage documentation.
    
  
  For performance analysis, you may be further interested in setting the CMake flag, 
  BENCHMARKONLY=1. Read more under "Custom Build and Installation".
