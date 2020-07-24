add_definitions(-DHISTOGRAM ) #TODO seems this should be moved into Options.cmake or somewhere else

FOREACH(MAX_NUM_READS ${MAX_NUM_READS_LIST})
	foreach(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})
		set(NAME overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		get_kmer_defs(KMER_DEFS)
	
		add_executable(overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN} ${CMAKE_SOURCE_DIR}/readoverlap/main.cpp
			${CMAKE_SOURCE_DIR}/readoverlap/ReadOverlapper.cpp
			$<TARGET_OBJECTS:OPTLIST> 
			$<TARGET_OBJECTS:MPIType>
			$<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}>
			$<TARGET_OBJECTS:MURMURHASH> # required for KmerObjects
			$<TARGET_OBJECTS:HASH_FUNCS> # required for KmerObjects
			$<TARGET_OBJECTS:readufx-${MAX_NUM_READS}-${HIPMER_KMER_LEN}>
		)
		add_dependencies(overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
		set_target_properties(overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN}
	        PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};${HIPMER_MPI_DEFINITIONS}"
	        COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
	        LINK_FLAGS "${MPI_CXX_LINK_FLAGS}"
	    )
	    target_link_libraries(overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN} ${MPI_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} MPI_FQ_OBJS Buffer MemoryChk ReadIdToProcId)
	    add_kmer_defs(overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
	    add_readid_defs(overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
	    	set(myTargets ${myTargets} overlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
	
	endforeach() # HIPMER_KMER_LENGTHS
ENDFOREACH() # MAX_NUM_READS

FOREACH(targetName ${myTargets})
    install(TARGETS ${targetName} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/ )
ENDFOREACH(targetName)
