#----------------------------------------------------------------------------------------#
#
#        Handles the CMake options
#
#----------------------------------------------------------------------------------------#

include(MacroUtilities)
include(Compilers)

# options (always available)
add_option(${PROJECT_NAME}_USE_ARCH "Enable architecture specific flags" OFF)
add_option(${PROJECT_NAME}_USE_SANITIZER "Enable sanitizer" OFF)
add_option(${PROJECT_NAME}_USE_COVERAGE "Enable code coverage" OFF)
add_option(BUILD_STATIC_LIBS "Build static libraries" OFF)
add_option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# default build type for libraries
set(DEFAULT_LIBRARY_TYPE STATIC)
if(BUILD_SHARED_LIBS)
    set(BUILD_LIBRARY_TYPE SHARED)
endif()

# features
add_feature(CMAKE_BUILD_TYPE "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)")
add_feature(CMAKE_INSTALL_PREFIX "Installation prefix")
add_feature(${PROJECT_NAME}_C_FLAGS "C compiler flags")
add_feature(${PROJECT_NAME}_CXX_FLAGS "C++ compiler flags")
add_feature(CMAKE_C_STANDARD "C languae standard")
add_feature(CMAKE_CXX_STANDARD "C++11 STL standard")
add_feature(${PROJECT_NAME}_DEFINITIONS "Preprocessor definitions")
add_feature(GIT_PREFIX "Git prefix")
set(GIT_PREFIX DIBELLA)

# build options
add_option(BENCHMARKONLY "Enable benchmarking (disable superfluous output)" OFF)
add_option(TIGHT_READS "Enable tight reads" OFF)
add_option(KMERA_ONLY "Enable K-mera only" OFF)
add_option(NO_UPC_MEMORY_POOL "Disable UPC memory pool (deprecated)" OFF)
add_option(BROKEN_ALLOCATOR_REBIND "(deprecated)" OFF)

#
# to change any of below, just run cmake with "-D VAR=VALUE"
#
add_cache_feature(HIPMER_MAX_KMER_SIZE "32" STRING "Max kmer size") # was 256 for release, 160 for debug
#add_cache_feature(MAX_MAX_NUM_READS "8" STRING "Max number of reads")
add_cache_feature(MAX_NUM_READS_LIST "4;6;7;8;9;16" STRING "List of maximum number of reads per build")
add_cache_feature(HIPMER_MAX_FILE_PATH "256" STRING "")
add_cache_feature(HIPMER_MAX_READ_NAME "512" STRING "")
add_cache_feature(HIPMER_MAX_LIBRARIES "128" STRING "")
add_cache_feature(HIPMER_LIB_NAME_LEN "15" STRING "")
add_feature(HIPMER_KMER_LENGTHS "")
add_feature(HIPMER_BLOOM "")
add_feature(HIPMER_BLOOM64 "")

math(EXPR HIPMER_MAX_KMER_SIZE "((${HIPMER_MAX_KMER_SIZE} + 31) / 32) * 32")
if (NOT HIPMER_KMER_LENGTHS)
    set(tmplen 32)
    while( NOT ${HIPMER_MAX_KMER_SIZE} LESS ${tmplen} )
      set(HIPMER_KMER_LENGTHS ${HIPMER_KMER_LENGTHS} ${tmplen} )
      math(EXPR tmplen "(${tmplen} + 32)")
    endwhile()
endif()

#if(NOT MAX_NUM_READS_LIST)
#    set(tmplen 8)
#    while( NOT ${MAX_MAX_NUM_READS} LESS ${tmplen} )
#      set(MAX_NUM_READS_LIST ${MAX_NUM_READS_LIST} ${tmplen} )
#      math(EXPR tmplen "(${tmplen} * 2)")
#    endwhile()
#endif ()

if("${HIPMER_BLOOM}" STREQUAL "64")
    set(HIPMER_BLOOM64 1)
endif()

# SET MPI / UPC definitions and variables #TODO remove superfluous UPC definitions...
# default preprocessor definitions
set(${PROJECT_NAME}_DEFINITIONS
    PROFILE
    CONFIG_USE_COLORS
    CONFIG_SANITY_CHECK
    CONFIG_SHOW_PROGRESS
    CONFIG_CHECK_SEQS
    MAX_READ_NAME_LEN=${HIPMER_MAX_READ_NAME}
    MAX_FILE_PATH=${HIPMER_MAX_FILE_PATH}
    MAX_LIBRARIES=${HIPMER_MAX_LIBRARIES}
    LIB_NAME_LEN=${HIPMER_LIB_NAME_LEN}
    NO_PAD
    NO_GZIP=1)
#    __UPCXX_VERSION__)

# optional definitions
foreach(_OPT BENCHMARKONLY TIGHT_READS KMERA_ONLY NO_UPC_MEMORY_POOL)
    if(${_OPT})
        list(APPEND ${PROJECT_NAME}_DEFINITIONS ${_OPT})
    endif()
endforeach()

# inconsistent setting of preprocessor definition
if (BROKEN_ALLOCATOR_REBIND)
    add_definitions(-DBROKEN_ALLOCATOR_REBIND=${BROKEN_ALLOCATOR_REBIND})
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    list(APPEND ${PROJECT_NAME}_DEFINITIONS DEBUG=1 SEQAN_ENABLE_DEBUG=1 SEQAN_HAS_EXECINFO=0)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	list(APPEND ${PROJECT_NAME}_DEFINITIONS SEQAN_ENABLE_TESTING=0 SEQAN_ENABLE_DEBUG=0 SEQAN_HAS_EXECINFO=0)
endif()

#----------------------------------------------------------------------------------------#
#
#        UPCXX installation options
#
#----------------------------------------------------------------------------------------#
#
#set(UPCXX_CODE_MODES O3 debug)
#set(UPCXX_GASNET_MODES par seq)
#set(UPCXX_BACKEND_MODES mpi udp smp)

#if(NOT UPCXX_CODEMODE)
#    set(_CODEMODE "$ENV{UPCXX_CODEMODE}" )
#    if(NOT _CODEMODE AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
#        set(_CODEMODE "debug")
#    elseif(NOT _CODEMODE)
#        set(_CODEMODE "O3")
#    endif()
#endif()

#add_cache_feature(UPCXX_GASNET "par" STRING "UPCXX gasnet modes (choices: [${UPCXX_GASNET_MODES}])")
#add_cache_feature(UPCXX_BACKEND "mpi" STRING "UPCXX backend modes (choices: [${UPCXX_BACKEND_MODES}])")
#add_cache_feature(UPCXX_CODEMODE "${_CODEMODE}" STRING "UPCXX codemode (choices: [${UPCXX_CODE_MODES}])")

# if all of these were specified, place this at the beginning of search options
#set(UPCXX_SEARCH_PATHS upcxx.${UPCXX_CODEMODE}.gasnet_${UPCXX_GASNET}.${UPCXX_BACKEND})
#set(UPCXX_PREFIX_PATHS )
#foreach(_PREFIX ${CMAKE_PREFIX_PATH})
#    list(APPEND UPCXX_PREFIX_PATHS ${_PREFIX}/upcxx.${UPCXX_CODEMODE}.gasnet_${UPCXX_GASNET}.${UPCXX_BACKEND})
#endforeach()

# loop over all options and build a search path if above is not set
# but prefix the list with the cache values to give preferential searching
# for those settings (end up with some duplicates but not a big deal)
#
# the "default" will be O3 + par + mpi
#
#foreach(_CODE ${UPCXX_CODEMODE} ${UPCXX_CODE_MODES})
#    foreach(_GASNET ${UPCXX_GASNET} ${UPCXX_GASNET_MODES})
#        foreach(_BACKEND ${UPCXX_BACKEND} ${UPCXX_BACKEND_MODES})
#            list(APPEND UPCXX_SEARCH_PATHS upcxx.${_CODE}.gasnet_${_GASNET}.${_BACKEND})
#            foreach(_PREFIX ${CMAKE_PREFIX_PATH})
#                list(APPEND UPCXX_PREFIX_PATHS ${_PREFIX}/upcxx.${_CODE}.gasnet_${_GASNET}.${_BACKEND})
#            endforeach()
#        endforeach()
#    endforeach()
#endforeach()

#set(UPCXX_CURRENT_CONFIG "${UPCXX_GASNET}.${UPCXX_BACKEND}.${UPCXX_CODEMODE}")
#set(UPCXX_LAST_CONFIG "${UPCXX_GASNET}.${UPCXX_BACKEND}.${UPCXX_CODEMODE}" CACHE STRING "Last configuration")


