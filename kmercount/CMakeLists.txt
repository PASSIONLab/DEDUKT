#TODO shouldn't need to add the MPI compile or link flags or include directory explicitly to each target - should be handled by Packages.cmake

#TODO move
find_package(OpenMP)
if (OPENMP_FOUND)
#TODO shouldn't need this...
  set(OpenMP_EXTRA_LINKER_FLAGS "${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
message("Set OpenMP_EXTRA_LINKER_FLAGS=${OpenMP_EXTRA_LINKER_FLAGS}")


# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif()

# TODO shouldn't need this here...
if(DISCOVER_LIBC)
  find_library(libc c)
else()
  set(libc "")
endif()

# TODO move to options...?
if (HIPMER_KHASH)
  add_definitions(-DKHASH)
endif()

add_subdirectory(DataStructures/libbloom)



set(NAME MPIType)
create_library(
	TYPE			OBJECT
	TARGET_NAME		${NAME}
	OUTPUT_NAME		${NAME}
	SOURCES			MPIType.cpp
	HEADERS			MPIType.h
	CXXFLAGS		${MPI_CXX_COMPILE_FLAGS}
)

# TODO not sure where to move these
add_definitions(-DHEAVYHITTERS -DHISTOGRAM )
if (MAXHITTERS)
  add_definitions(-DMAXHITTERS=${MAXHITTERS})
endif()


# KmerObjects
FOREACH(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})
	set(NAME KmerObjects-${HIPMER_KMER_LEN})
	get_kmer_defs(KMER_DEFS)
	create_library(
		TYPE			OBJECT
		TARGET_NAME		${NAME}
		OUTPUT_NAME		${NAME}
		HEADERS			Kmer.hpp KmerIterator.hpp Kmer.hpp # includes common/common.h and common/hash_funcs.h
		SOURCES			Kmer.cpp KmerIterator.cpp
		DEFINITIONS		${KMER_DEFS}
	)
ENDFOREACH()

add_library(gpu_kmer STATIC
	simple.cu 
	linearprobing.cu
 	simple.cuh 
 	linearprobing.h 
  )

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(gpu_kmer PUBLIC )

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties( gpu_kmer
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON
                      )


FOREACH(MAX_NUM_READS ${MAX_NUM_READS_LIST})
	FOREACH(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})
		get_kmer_defs(KMER_DEFS)

#readufx
		set(NAME readufx-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		create_library(
			TYPE				OBJECT
			TARGET_NAME			${NAME}
			OUTPUT_NAME			${NAME}
object			SOURCES				readufx.cpp
			HEADERS				readufx.h
			DEFINITIONS			${KMER_DEFS}
			LINK_LIBRARIES		${ZLIB_LIBRARIES}
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS})
		add_readid_defs(${NAME})

#kmermatch
		set(NAME kmermatch-${MAX_NUM_READS}-${HIPMER_KMER_LEN}-NOUPC)
		create_library(
			TYPE				OBJECT
			TARGET_NAME			${NAME}
			OUTPUT_NAME			${NAME}
			SOURCES				kmermatch.cpp
			DEFINITIONS			${KMER_DEFS}
			LINK_LIBRARIES		${ZLIB_LIBRARIES}
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS})
		add_readid_defs(${NAME})
		
# b2tufx
		set(NAME b2tufx-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		create_executable(
			INSTALL
			TARGET_NAME				${NAME}
			OUTPUT_NAME				${NAME}
			SOURCES					bin2textUFX.cpp
			DEFINITIONS				${KMER_DEFS}
			CXXFLAGS				${MPI_CXX_COMPILE_FLAGS}
			OBJECT_LIBRARIES		common fqreader KmerObjects-${HIPMER_KMER_LEN}
			LINK_LIBRARIES			${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} 
			INSTALL_DESTINATION 	${CMAKE_INSTALL_PREFIX}/bin
		)
		add_readid_defs(${NAME})
	
# b2tcnt
 		set(NAME b2tcnt-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		create_executable(
			INSTALL
			TARGET_NAME			${NAME}
			OUTPUT_NAME			${NAME}
			SOURCES				bin2textCount.cpp
			DEFINITIONS			${KMER_DEFS}
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS}
			OBJECT_LIBRARIES	common KmerObjects-${HIPMER_KMER_LEN}
			LINK_LIBRARIES		${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} 
			INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
		)
		add_readid_defs(${NAME})
		
	#    IF (HIPMER_FULL_BUILD)
	#      # SNPidentify
	#      add_executable(SNPidentify-${HIPMER_KMER_LEN} SNPidentify.cpp $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}> $<TARGET_OBJECTS:HASH_FUNCS>)
	##      add_dependencies(SNPidentify-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
	#      set_target_properties(SNPidentify-${HIPMER_KMER_LEN} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
	#      target_link_libraries(SNPidentify-${HIPMER_KMER_LEN} ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc})
	#      add_kmer_defs(SNPidentify-${HIPMER_KMER_LEN})
	#      set(myTargets ${myTargets} SNPidentify-${HIPMER_KMER_LEN})
	#    ENDIF()
	
#overlaps_objs / _obs
		set(NAME ReadOverlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		create_library(
			TYPE 				OBJECT
			TARGET_NAME			${NAME}
			OUTPUT_NAME			${NAME}
			SOURCES				${PROJECT_SOURCE_DIR}/readoverlap/ReadOverlapper.cpp ${PROJECT_SOURCE_DIR}/readoverlap/main.cpp
			HEADERS				${PROJECT_SOURCE_DIR}/readoverlap/ReadOverlapper.h ${PROJECT_SOURCE_DIR}/readoverlap/overlap.h
			DEFINITIONS			${KMER_DEFS} SINGLE_EXEC=1
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS})
		add_readid_defs(${NAME})
		
# kmermatch
		set(NAME kmermatch-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		create_executable(
			INSTALL
			TARGET_NAME			${NAME}
			OUTPUT_NAME			${NAME}
			DEFINITIONS			${KMER_DEFS}
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS} ${CMAKE_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}
			OBJECT_LIBRARIES	common fqreader MURMURHASH readufx-${MAX_NUM_READS}-${HIPMER_KMER_LEN} ReadOverlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN} LIBBLOOM-NOUPC MPIType KmerObjects-${HIPMER_KMER_LEN} kmermatch-${MAX_NUM_READS}-${HIPMER_KMER_LEN}-NOUPC
			LINK_LIBRARIES		${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} gpu_kmer  
			INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
		)
		add_readid_defs(${NAME})

#kmermatch_objs for unified exec        
		set(NAME kmermatch_objs-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		create_library(
			TYPE 				OBJECT
			TARGET_NAME			${NAME}
			OUTPUT_NAME			${NAME}
			SOURCES				kmermatch.cpp
			DEFINITIONS			${KMER_DEFS} SINGLE_EXEC=1
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS} ${CMAKE_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}
#			OBJECT_LIBRARIES	common fqreader ReadOverlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN} LIBBLOOM-NOUPC MPIType KmerObjects-${HIPMER_KMER_LEN} kmermatch-${MAX_NUM_READS}-${HIPMER_KMER_LEN}-NOUPC
			LINK_LIBRARIES		${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} gpu_kmer)
		add_readid_defs(${NAME})
	
	ENDFOREACH() # HIPMER_KMER_LENGTHS
ENDFOREACH() # MAX_NUM_READS

#set_source_files_properties(check_hugepages.c PROPERTIES LANGUAGE "C" )
#ADD_EXECUTABLE(check_hugepages check_hugepages.c)
#SET_TARGET_PROPERTIES(check_hugepages PROPERTIES
#                          COMPILE_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS} ${OpenMP_C_FLAGS}"
#                          LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OpenMP_EXTRA_LINKER_FLAGS}")
#TARGET_LINK_LIBRARIES(check_hugepages  ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc})
#INSTALL(TARGETS check_hugepages DESTINATION  ${CMAKE_INSTALL_PREFIX}/bin/ )

set(NAME check_hugepages)
create_executable(
	INSTALL
	TARGET_NAME         ${NAME}
	OUTPUT_NAME         ${NAME}
	SOURCES				check_hugepages.c
	DEFINITIONS			USE_MPI=1
	PROPERTIES			LANGUAGE C LINKER_LANGUAGE C LINK_FLAGS ${OpenMP_EXTRA_LINKER_FLAGS}
	CFLAGS				${MPI_C_COMPILE_FLAGS} ${OpenMP_C_FLAGS} 
	INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}
	LINK_LIBRARIES		${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} ${MPI_C_LIBRARIES} #gpu_kmer
	INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
