FILE(GLOB headers ${CMAKE_CURRENT_LIST_DIR}/*.h ${CMAKE_CURRENT_LIST_DIR}/*.hpp)

create_library(TARGET_NAME common
    TYPE OBJECT     # object library
    OUTPUT_NAME common
    HEADERS ${headers}
    SOURCES memory_chk.c Buffer.c hash_funcs.c optlist.c ReadIdToProcId.cpp FileOpenerMPI.cpp
    INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/common
    LINK_LIBRARIES ${ZLIB_LIBRARIES}
)