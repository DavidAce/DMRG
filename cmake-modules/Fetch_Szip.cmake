
message(STATUS "SZIP will be installed into ${INSTALL_DIRECTORY}/szip on first build.")
include(ExternalProject)
ExternalProject_Add(library_SZIP
        URL      https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
        PREFIX              "${INSTALL_DIRECTORY}/szip"
        UPDATE_DISCONNECTED 1
        TEST_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DSZIP_ENABLE_ENCODING:BOOL=ON
        )

ExternalProject_Get_Property(library_SZIP INSTALL_DIR)
add_library(SZIP INTERFACE)
add_dependencies(SZIP library_SZIP)
set(SZIP_LIBRARY        ${INSTALL_DIR}/lib/libszip-static.a)
set(SZIP_INCLUDE_DIRS   ${INSTALL_DIR}/include)

set_target_properties(SZIP PROPERTIES
        TARGET_LINK_LIBRARIES         "${SZIP_LIBRARY}"
        TARGE_INCLUDE_DIRECTORIES     "${SZIP_INCLUDE_DIRS}"
        )

target_link_libraries(${PROJECT_NAME} PRIVATE SZIP)
