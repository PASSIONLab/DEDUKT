include_directories(${MPI_C_INCLUDE_PATH})
SET_SOURCE_FILES_PROPERTIES(main.cpp exec_utils.c config.c PROPERTIES LANGUAGE "CXX")

if (NOT SEQAN_INCLUDE_PATH)
    set(SEQAN_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../libs/seqan-old)
endif()

if (NOT EDLIB_INCLUDE_PATH)
    set(EDLIB_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../libs/edlib)
endif()

# if ( CMAKE_COMPILER_IS_GNUCC )
# set(MPI_CXX_COMPILE_FLAGS  "${MPI_CXX_COMPILE_FLAGS} -Wall -Wextra")
# endif()

FOREACH(MAX_NUM_READS ${MAX_NUM_READS_LIST})
  foreach(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})
  	set(NAME main-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
  	get_kmer_defs(KMER_DEFS)
  	create_executable(
  		INSTALL
  		TARGET_NAME				${NAME}
  		OUTPUT_NAME				${NAME}
  		SOURCES					main.cpp exec_utils.c config.c
  		HEADERS					config.h exec_utils.h
  		DEFINITIONS				${KMER_DEFS}
  		PROPERTIES				LANGUAGE CXX LINKER_LANGUAGE CXX
  		INCLUDE_DIRECTORIES		${ARGP_INCLUDE_DIRS} ${SEQAN_INCLUDE_PATH} # do I really need seqan here again?
  		OBJECT_LIBRARIES		common fqreader MPIType	LIBBLOOM MURMURHASH loadfq_objs KmerObjects-${HIPMER_KMER_LEN} kmermatch_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN} readufx-${MAX_NUM_READS}-${HIPMER_KMER_LEN} ReadOverlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN} align_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN}
  		LINK_LIBRARIES			cuda_objs  ${ARGP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} 
  		CXXFLAGS				${MPI_CXX_COMPILE_FLAGS} ${MPI_CXX_LINK_FLAGS} 
  		INSTALL_DESTINATION		${CMAKE_INSTALL_PREFIX}/bin
  	)
    add_readid_defs(${NAME}) #TODO not sure this actually works

#    add_executable(main-${MAX_NUM_READS}-${HIPMER_KMER_LEN}
#		$<TARGET_OBJECTS:main_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN}>
#        $<TARGET_OBJECTS:loadfq_objs>
#		 $<TARGET_OBJECTS:OPTLIST>
#        $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}> 
#        $<TARGET_OBJECTS:KmerHandling-${HIPMER_KMER_LEN}> 
#        $<TARGET_OBJECTS:readufx-${MAX_NUM_READS}-${HIPMER_KMER_LEN}>  
#        $<TARGET_OBJECTS:kmermatch_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN}>
#        $<TARGET_OBJECTS:MURMURHASH>
#        $<TARGET_OBJECTS:LIBBLOOM>
#        $<TARGET_OBJECTS:MPIType>
#        $<TARGET_OBJECTS:_MPI_FQ_OBJS> #TODO... cmake doesn't recognize MPI_FQ_OBS
#        $<TARGET_OBJECTS:HASH_FUNCS>
#        $<TARGET_OBJECTS:CPU_AFFINITY>
#        $<TARGET_OBJECTS:FASTA>
#		 $<TARGET_OBJECTS:overlap_obs-${MAX_NUM_READS}-${HIPMER_KMER_LEN}>
#		 $<TARGET_OBJECTS:align_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN}>
#    )
	
#	set_target_properties(main-${MAX_NUM_READS}-${HIPMER_KMER_LEN} 
#        PROPERTIES LINKER_LANGUAGE "CXX" 
#        COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
#        LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")

  endforeach() # HIPMER_KMER_LENGTHS
endforeach() # MAX_NUM_READS