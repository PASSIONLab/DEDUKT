################################################################################
#
#        Handles the external packages
#
################################################################################

include(FindPackageHandleStandardArgs)


################################################################################
#
#        MPI
#
################################################################################

    # MPI C compiler from environment
    set(_ENV MPICC)
    if(NOT DEFINED MPI_C_COMPILER AND NOT "$ENV{${_ENV}}" STREQUAL "")
        message(STATUS "Setting MPI C compiler to: $ENV{${_ENV}}")
        set(MPI_C_COMPILER $ENV{${_ENV}} CACHE FILEPATH "MPI C compiler")
    endif()

    # MPI C++ compiler from environment
    set(_ENV MPICC)
    if(NOT DEFINED MPI_CXX_COMPILER AND NOT "$ENV{${_ENV}}" STREQUAL "")
        message(STATUS "Setting MPI C++ compiler to: $ENV{${_ENV}}")
        set(MPI_CXX_COMPILER $ENV{${_ENV}} CACHE FILEPATH "MPI C++ compiler")
    endif()
    unset(_ENV)

    find_package(MPI REQUIRED)
    set(MPI_LIBRARIES )

        # Add the MPI-specific compiler and linker flags
        to_list(_FLAGS "${MPI_C_COMPILE_FLAGS}")
        foreach(_FLAG ${_FLAGS})
            add_c_flag_if_avail("${_FLAG}")
        endforeach()
        unset(_FLAGS)
        to_list(_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
        foreach(_FLAG ${_FLAGS})
            add_cxx_flag_if_avail("${_FLAG}")
            message(STATUS "checking ${_FLAG}")
        endforeach()
        unset(_FLAGS)
        add(CMAKE_EXE_LINKER_FLAGS "${MPI_CXX_LINK_FLAGS}")
        list(APPEND EXTERNAL_INCLUDE_DIRS
            ${MPI_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})

        foreach(_TYPE C_LIBRARIES CXX_LIBRARIES EXTRA_LIBRARY)
            set(_TYPE MPI_${_TYPE})
            if(${_TYPE})
                list(APPEND EXTERNAL_LIBRARIES ${${_TYPE}})
            endif()
        endforeach()

        if(NOT MPIEXEC_EXECUTABLE AND MPIEXEC)
            set(MPIEXEC_EXECUTABLE ${MPIEXEC} CACHE FILEPATH "MPI executable")
        endif()

        if(NOT MPIEXEC_EXECUTABLE AND MPI_EXECUTABLE)
            set(MPIEXEC_EXECUTABLE ${MPI_EXECUTABLE} CACHE FILEPATH "MPI executable")
        endif()

################################################################################
#
#        Git
#
################################################################################
set(DIBELLA_GIT_VERSION_FILE ${CMAKE_SOURCE_DIR}/DIBELLA_VERSION)

find_package(Git QUIET)
if (GIT_FOUND AND IS_DIRECTORY ${CMAKE_SOURCE_DIR}/.git)
	add_custom_target(GET_GIT_VERSION ALL
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD > ${DIBELLA_GIT_VERSION_FILE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Retrieving git version")
endif()

set(DIBELLA_VERSION_FILE_TEMPLATE ${CMAKE_SOURCE_DIR}/version.h.in)
set(DIBELLA_VERSION_FILE ${CMAKE_SOURCE_DIR}/version.h)
set(${GIT_PREFIX}_GIT_VERSION_FILE DIBELLA_VERSION)
configure_file(${CMAKE_SOURCE_DIR}/makeVersionFile.cmake.in ${CMAKE_BINARY_DIR}/makeVersionFile.cmake @ONLY)
add_custom_target(REPLACE_VERSION_H ALL
		COMMAND ${CMAKE_COMMAND}
			-DDIBELLA_GIT_VERSION_FILE=${DIBELLA_GIT_VERSION_FILE}
			-DDIBELLA_VERSION_FILE=${DIBELLA_VERSION_FILE}
			-DDIBELLA_VERSION_FILE_TEMPLATE=${DIBELLA_VERSION_FILE_TEMPLATE}
			-P ${CMAKE_BINARY_DIR}/makeVersionFile.cmake
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		DEPENDS GET_GIT_VERSION
		COMMENT "Building ${DIBELLA_VERSION_FILE}")
add_dependencies(REPLACE_VERSION_H GET_GIT_VERSION)

################################################################################
#
#        zlib
#
################################################################################
find_package(ZLIB 1.2.3 REQUIRED)
list(APPEND EXTERNAL_INCLUDE_DIRS ${ZLIB_INCLUDE_DIRS})
list(APPEND EXTERNAL_LIBRARIES} ${ZLIB_LIBRARIES})

find_package(MPI REQUIRED)

################################################################################
#
#        External variables
#
################################################################################

# including the directories
safe_remove_duplicates(EXTERNAL_INCLUDE_DIRS ${EXTERNAL_INCLUDE_DIRS})
safe_remove_duplicates(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES})
foreach(_DIR ${EXTERNAL_INCLUDE_DIRS})
    include_directories(SYSTEM ${_DIR})
endforeach()

# include dirs
set(${PROJECT_NAME}_INCLUDE_DIRECTORIES ${EXTERNAL_INCLUDE_DIRS})

# link libs
set(${PROJECT_NAME}_LINK_LIBRARIES ${EXTERNAL_LIBRARIES})

set(${PROJECT_NAME}_PROPERTIES
    C_STANDARD                  11
    C_STANDARD_REQUIRED         ON
    CXX_STANDARD                14
    CXX_STANDARD_REQUIRED       ON
)

set(CUDA_PROPERTIES
    CUDA_STANDARD               14
    CUDA_STANDARD_REQUIRED      ON
    CUDA_SEPARABLE_COMPILATION  ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    CUDA_USE_STATIC_CUDA_RUNTIME OFF)
