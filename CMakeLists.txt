cmake_minimum_required (VERSION 3.5 FATAL_ERROR) # was 2.8

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    set(MSG "")
    message(STATUS "Warning! Building from the source directory is not recommended")
    message(STATUS "If unintented, please remove 'CMakeCache.txt' and 'CMakeFiles'")
    message(STATUS "and build from a separate directory")
    message(WARNING "In-source build")
endif()

project(dibella LANGUAGES C CXX CUDA)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

include(MacroUtilities)
include(Compilers)
include(Options)
include(GNUInstallDirs)
include(ProjectSettings)
include(BuildSettings)
include(Packages)
include(Macros)

add_subdirectory(common)
add_subdirectory(fqreader)
add_subdirectory(loadfq)
add_subdirectory(kmercount)
#add_subdirectory(readoverlap) # not used after combining into kmermatch
# add_subdirectory(align)
#add_subdirectory(main)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/scripts
	DESTINATION bin
	PATTERN "scripts/*"
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
		GROUP_EXECUTE GROUP_READ
		WORLD_EXECUTE WORLD_READ)

print_features()
