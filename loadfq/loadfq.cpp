#include "../common/optlist.h"
#include "../common/mpi_common.h"
#include "../common/MPIUtils.h"
#include "../fqreader/fq_reader.h"
#include "../common/defines.h"

#include "loadfq.c" // such a hackish way of getting cmake to build fq_reader as a c++ object (wouldn't cooperate with the not hackish methods)
