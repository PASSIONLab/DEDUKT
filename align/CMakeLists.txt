#
# pairwise alignment libraries 
#
if (NOT SEQAN_INCLUDE_PATH)
    set(SEQAN_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../libs/seqan-old)
endif()
list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${SEQAN_INCLUDE_PATH} ${ARGP_INCLUDE_DIRS})

if (NOT EDLIB_INCLUDE_PATH)
    set(EDLIB_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../libs/edlib)
endif()
list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${EDLIB_INCLUDE_PATH} ${ARGP_INCLUDE_DIRS})
###########################################

FOREACH(MAX_NUM_READS ${MAX_NUM_READS_LIST})
	FOREACH(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})

		set(NAME alignment-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		get_kmer_defs(KMER_DEFS)
        create_executable(
            INSTALL
            TARGET_NAME         ${NAME}
            OUTPUT_NAME         ${NAME}
            SOURCES             ${CMAKE_SOURCE_DIR}/align/main.cpp ${CMAKE_SOURCE_DIR}/align/Aligner.cpp 
            HEADERS				align.h Aligner.h
            DEFINITIONS         MAX_NUM_READS=${MAX_NUM_READS} ${KMER_DEFS}
            CXXFLAGS			${MPI_CXX_COMPILE_FLAGS}
            INCLUDE_DIRECTORIES ${ARGP_INCLUDE_DIRS} ${SEQAN_INCLUDE_PATH} ${EDLIB_INCLUDE_PATH}
            LINK_LIBRARIES      ${RT_LIBRARIES} ${libc} ${MPI_LIBRARIES} ${ZLIB_LIBRARIES}  
            OBJECT_LIBRARIES    common fqreader #KmerObjects-${HIPMER_KMER_LEN} MURMURHASH 
            INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/bin) # not sure if MPI_CXX_COMPILE_FLAGS and MPI_CXX_LINK_FLAGS need to be explicitly added here or in parent
		add_readid_defs(alignment-${MAX_NUM_READS}-${HIPMER_KMER_LEN})

		# library object to be linked with main/main without -DSINGLE_EXEC, maybe there's a better way of doing this
	    #add_library(align_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN} OBJECT ${EDLIB_INCLUDE_PATH}/src/edlib.cpp ${CMAKE_SOURCE_DIR}/align/Aligner.cpp ${CMAKE_SOURCE_DIR}/align/main.cpp)
        set(NAME align_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
        create_library(
            TYPE				OBJECT
            TARGET_NAME         ${NAME}
            OUTPUT_NAME         ${NAME}
            SOURCES             ${CMAKE_SOURCE_DIR}/align/main.cpp ${CMAKE_SOURCE_DIR}/align/Aligner.cpp 
            DEFINITIONS         MAX_NUM_READS=${MAX_NUM_READS} ${KMER_DEFS} SINGLE_EXEC=1
            INCLUDE_DIRECTORIES ${ARGP_INCLUDE_DIRS} ${SEQAN_INCLUDE_PATH} ${EDLIB_INCLUDE_PATH}
            LINK_LIBRARIES      ${RT_LIBRARIES} ${libc} ${MPI_LIBRARIES} ${ZLIB_LIBRARIES}
         )
		add_readid_defs(${NAME})

	ENDFOREACH()
ENDFOREACH()
