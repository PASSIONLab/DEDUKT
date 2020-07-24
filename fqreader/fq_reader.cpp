#include "../common/common.h"
#include "../common/memory_chk.h"
#include "../common/mpi_common.h"

extern StaticVars _sv;

//#include "fq_reader.h"
#include "fq_reader.c" // such a hackish way of getting cmake to build fq_reader as a c++ object (wouldn't cooperate with the not hackish methods)
