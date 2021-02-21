A distributed-memory parallel *K*-mer counter with NVIDIA GPU support. 

## Build requirements:
- GCC Compiler (tested with 6.4.0)
- CUDA SDK (tested with 10.2)
- Message Passing Interface (MPI) (tested with spectrum-mpi-10.3.1)
- CMake >= 3.5 (C++ 11 support required)

## Building and Installing

From the source directory run (the following example is for summit):
```
$  sh install.sh scripts/install-environments/summitgpu-env.sh
```
To build the Debug version run:
```
$  DEBUG=1 sh install.sh scripts/install-environments/summitgpu-env.sh
```
Some variables you may want to redefine are:
BUILD_DIR     :    the path to your desired build directory
INSTALL_DIR   :    the path to your desired build directory
CC            :    the C compiler
CXX           :    the C++ compiler 
More details can be founf at the bottom of this README.


## Binaries 
Executable binaries - kmermatch-<u_max>-<k_max> and scripts can be found in $INSTALL_DIR. 
Each of <u_max> and <k_max> are substituted for various integers in the actual installation.
The former, <u_max>, corresponds to the maximumm reliable k-mer upper limit. 
The pipeline (or pipeline stage) can be run with any -u value less than or equal to <u_max> with the respective binary.
Similarly, <k_max>, corresponds to the maximum value of k (the next higher power of 2 e.g. 32 for k=17) with which the binary can be run.

To compile the code base with values of <u_max> other than the default, specify a list of semi-colon separated values using the cmake flag "-DMAX_NUM_READS_LIST".
You can do so on the command line (example below) or in your environment script (see the "Custom Build and Installation" section).
```
$  DEFS="-DMAX_NUM_READS_LIST=8;16" sh install.sh 
```


## Run in interarctive session

Examples are based on summit cluster from ORNL. To run on 64 nodes (6 GPUs each):  

Example #1 (K-mer counter on GPU)  
```
$ jsrun --smpiargs=-gpu -n64 -c42 -g6 -a6 build/kmermatch-16-32 -k 17 -t 1 -u 8 -c 1 -i ${input} 
```
Example #2 (Supermer based K-mer counter on GPU)
```
$ jsrun --smpiargs=-gpu -n64 -c42 -g6 -a6 build/kmermatch-16-32 -k 17 -t 3 -n 9 -u 8 -c 1 -i ${input} 
```
```
Options:
-t: type of kmer counter:   
    0: kmer counter on CPU  
    1: kmer counter on NVIDIA GPU  
    2: supermer based kmer counter on CPU  
    3: supermer based kmer counter on NVIDIA GPU  
-k: K-mer length
-n: minimizer length
```
{input} is any text file containing a line-separated list of input fastq files. Use the -h flag with a "kmermatch" binary to view the respective usage documentation.
For performance analysis, you may be further interested in setting the CMake flag, 
BENCHMARKONLY=1. Read more under "Custom Build and Installation".
Additionally, to customize particular feature settings, add them to the cmake command as "-D<feature name1> -D<feature name 2> " and so on,
  or if using the install.sh script, as DEFS="-D<feature 1 name> -D<feature 2 name>" sh install.sh.


## Run as a batch job

Examples of running a batch job on summit using 64 nodes (384 GPUs, 6 GPUs each node) can be found in scripts/submit_64node.lsf

## More details on building and installing
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



