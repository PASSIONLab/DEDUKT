#ifndef COMMON_H
#define COMMON_H

// #include "Kmer.hpp"
// // #include "ParallelFASTQ.h"
// #include "Friends.h"
// #include "MPIType.h"
// #include "SimpleCount.h"
// #include "FriendsMPI.h"
// #include "Pack.h"
#include "../common/MPIUtils.h"
#include "../common/Buffer.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <libgen.h>
#include <algorithm>
#include <numeric>
#include <string>
#include <functional>

extern int nkmers_thisBatch;
extern int nkmers_processed;
extern double tot_alltoallv_GPU;
extern double tot_GPUsmer_alltoallv;

const int COMM_ITER=1;
const int BUFF_SCALE=4;

#endif