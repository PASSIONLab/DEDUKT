#TODO shouldn't need to add the MPI compile or link flags or include directory explicitly to each target - should be handled by Packages.cmake

#TODO move
find_package(OpenMP)
if (OPENMP_FOUND)
#TODO shouldn't need this...
  set(OpenMP_EXTRA_LINKER_FLAGS "${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
message("Set OpenMP_EXTRA_LINKER_FLAGS=${OpenMP_EXTRA_LINKER_FLAGS}")

set(CUDA_SEPARABLE_COMPILATION ON)
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
# 	supermer.cu 
# 	linearprobing.cu
#  	supermer.h 
# 	common.h
# 	common_gpu.h
#  	linearprobing.h 
  	simple.cpp
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
			SOURCES				kmermatch.cu supermer.cpp supermer.cu kmercounter_GPU.cu linearprobing.cu 
			HEADERS				common.h supermer.h common_gpu.h linearprobing.h
			DEFINITIONS			${KMER_DEFS}
			LINK_LIBRARIES		${ZLIB_LIBRARIES}
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS}
			CUDAFLAGS	        -ccbin=${CMAKE_CXX_COMPILER} -Xcompiler -fopenmp -arch=sm_70 -maxrregcount=32 #TODO some of these flags shouldn't be hard coded here
 			PROPERTIES          COMPILE_LANGUAGE CUDA ${CUDA_PROPERTIES} ) 
		add_readid_defs(${NAME})
			
		
# kmermatch
		set(NAME kmermatch-${MAX_NUM_READS}-${HIPMER_KMER_LEN})
		create_executable(
			INSTALL
			TARGET_NAME			${NAME}
			OUTPUT_NAME			${NAME}
			DEFINITIONS			${KMER_DEFS}
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS} ${CMAKE_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}
			OBJECT_LIBRARIES	common fqreader MURMURHASH readufx-${MAX_NUM_READS}-${HIPMER_KMER_LEN} LIBBLOOM-NOUPC MPIType KmerObjects-${HIPMER_KMER_LEN} kmermatch-${MAX_NUM_READS}-${HIPMER_KMER_LEN}-NOUPC
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
			SOURCES				kmermatch.cu supermer.cpp supermer.cu kmercounter_GPU.cu 
			HEADERS				common.h supermer.h common_gpu.h
			DEFINITIONS			${KMER_DEFS} SINGLE_EXEC=1
			CXXFLAGS			${MPI_CXX_COMPILE_FLAGS} ${CMAKE_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}
#			OBJECT_LIBRARIES	common fqreader ReadOverlap-${MAX_NUM_READS}-${HIPMER_KMER_LEN} LIBBLOOM-NOUPC MPIType KmerObjects-${HIPMER_KMER_LEN} kmermatch-${MAX_NUM_READS}-${HIPMER_KMER_LEN}-NOUPC
			LINK_LIBRARIES		${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} gpu_kmer
			CUDAFLAGS	        -ccbin=${CMAKE_CXX_COMPILER} -Xcompiler -fopenmp -arch=sm_70 -maxrregcount=32 #TODO some of these flags shouldn't be hard coded here
 			PROPERTIES          COMPILE_LANGUAGE CUDA ${CUDA_PROPERTIES} ) 
		add_readid_defs(${NAME})
	
	ENDFOREACH() # HIPMER_KMER_LENGTHS
ENDFOREACH() # MAX_NUM_READS


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
	LINK_LIBRARIES		${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${libc} ${MPI_C_LIBRARIES} gpu_kmer
	INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
