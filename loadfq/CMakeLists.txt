# fix linking when objects are from multiple languages
#set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "")
#set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "")
#set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++")
#set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "")

#SET(MPI_CODE loadfq.c unloadfq.c)
#SET_SOURCE_FILES_PROPERTIES(${MPI_CODE} PROPERTIES LANGUAGE "CXX" OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/common/mpi_common.h")
#SET_SOURCE_FILES_PROPERTIES(${MPI_CODE} PROPERTIES OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/common/mpi_common.h")
#SET_SOURCE_FILES_PROPERTIES(${MPI_CODE} PROPERTIES OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/fqreader/fq_reader.h")
#SET_SOURCE_FILES_PROPERTIES(${MPI_CODE} PROPERTIES OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/fqreader/fq_reader.c")

FILE(GLOB headers ${PROJECT_SOURCE_DIR}/common/*.h ${PROJECT_SOURCE_DIR}/common/*.hpp)

# The standalone loadfq executable isn't totally necessary
# and it currently has a nondeterministically occurring linking error
# where CMake generates "-o . " instead of "-o loadfq" despite 
# having the target and output names set.
#
#set(NAME loadfq)
#create_executable(
#    INSTALL
#	TARGET_NAME         ${NAME}
#	OUTPUT_NAME         ${NAME}
#	SOURCES             loadfq.cpp
#	HEADERS             loadfq.h #${headers}
#	PROPERTIES LANGUAGE CXX LINKER_LANGUAGE CXX
#	LINK_LIBRARIES      ${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${MPI_CXX_LIBRARIES}
#	OBJECT_LIBRARIES    common fqreader
#	INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/fqreader
#	INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
#)

set(NAME loadfq_objs)
create_library(
	TYPE OBJECT     		# object library
	TARGET_NAME				${NAME}
	OUTPUT_NAME				${NAME}
	HEADERS					loadfq.h ${headers}
	SOURCES					loadfq.cpp
	DEFINITIONS				SINGLE_EXEC=1
	LINK_LIBRARIES			${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${MPI_CXX_LIBRARIES}
	INCLUDE_DIRECTORIES 	${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/common ${PROJECT_SOURCE_DIR}/fqreader
)

# not currently used / supported
#set(NAME unloadfq)
#create_executable(
#	INSTALL
#	TARGET_NAME				${NAME}
#	OUTPUT_NAME				${NAME}
#	HEADERS					${headers}
#	SOURCES					unloadfq.c
#	LINK_LIBRARIES			${ZLIB_LIBRARIES} ${RT_LIBRARIES} ${MPI_CXX_LIBRARIES}
#	OBJECT_LIBRARIES		common fqreader
#	INCLUDE_DIRECTORIES 	${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/fqreader
#	INSTALL_DESTINATION 	${CMAKE_INSTALL_PREFIX}/bin)

