
message(STATUS "SZIP will be installed into ${INSTALL_DIRECTORY}/szip on first build.")
include(ExternalProject)
ExternalProject_Add(external_SZIP
        URL      https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
        PREFIX      ${BUILD_DIRECTORY}/szip
        INSTALL_DIR ${INSTALL_DIRECTORY}/szip
        UPDATE_DISCONNECTED 1
        TEST_COMMAND ""
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        -DSZIP_ENABLE_ENCODING:BOOL=ON
        )

ExternalProject_Get_Property(external_SZIP INSTALL_DIR)
add_library(SZIP INTERFACE)
add_dependencies(SZIP external_SZIP)
set(SZIP_LIBRARY        ${INSTALL_DIR}/lib/libszip-static.a)
set(SZIP_INCLUDE_DIRS   ${INSTALL_DIR}/include)
target_link_libraries(SZIP INTERFACE ${INSTALL_DIR}/lib/libszip-static.a)
target_include_directories(SZIP SYSTEM INTERFACE ${SZIP_INCLUDE_DIRS})

