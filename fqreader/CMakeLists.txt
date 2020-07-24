#add_definitions(-DFQ_READER_BUFFER_SIZE=${HIPMER_READ_BUFFER})
#include_directories(${MPI_CXX_INCLUDE_PATH})
#set(MPI_FQ_FILES fq_reader_mpi.cpp)
#set_source_files_properties(fq_reader_mpi.cpp PROPERTIES OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/common/mpi_common.h")
#add_library(_MPI_FQ_OBJS OBJECT ${MPI_FQ_FILES})
#add_library(MPI_FQ_OBJS $<TARGET_OBJECTS:_MPI_FQ_OBJS> $<TARGET_OBJECTS:_MemoryChk>)
#add_dependencies(MPI_FQ_OBJS REPLACE_VERSION_H)
#set_target_properties(MPI_FQ_OBJS PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}" )

#set(UPC_FQ_FILES fq_reader_upc.c)
#set_source_files_properties(${UPC_FQ_FILES} PROPERTIES LANGUAGE "UPC" OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/hipmer/common/upc_common.h")
#add_library(UPC_FQ_OBJS OBJECT ${UPC_FQ_FILES})
#add_dependencies(UPC_FQ_OBJS REPLACE_VERSION_H)
#SET_SOURCE_FILES_PROPERTIES(${UPC_FQ_FILES} PROPERTIES OBJECT_DEPENDS
#"${CMAKE_SOURCE_DIR}/src/hipmer/fqreader/fq_reader.c" OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/hipmer/common/common.h")

FILE(GLOB headers ${PROJECT_SOURCE_DIR}/common/*.h)

create_library(
	TARGET_NAME fqreader
    TYPE OBJECT     # object library
    OUTPUT_NAME fqreader
    SOURCES fq_reader.cpp # fq_reader.c
    HEADERS fq_reader.h ${headers}
    INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/common
#    OBJECT_LIBRARIES common
)