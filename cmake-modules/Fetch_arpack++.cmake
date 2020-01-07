


include(cmake-modules/FindArpack++.cmake)
find_Arpackpp()

if (TARGET arpack::arpack++)
    message(STATUS "Arpack++ found")
elseif(NOT ${DOWNLOAD_METHOD} MATCHES "none")
    message(STATUS "Arpack++ will be installed into ${CMAKE_BINARY_DIR}/dmrg-deps-install/arpack++ on first build.")
    include(ExternalProject)
    ExternalProject_Add(external_ARPACK++
            GIT_REPOSITORY      https://github.com/m-reuter/arpackpp.git
            GIT_TAG             2.3.0
            GIT_PROGRESS false
            GIT_SHALLOW true
            PREFIX      ${CMAKE_BINARY_DIR}/dmrg-deps-build/arpack++
            INSTALL_DIR ${CMAKE_BINARY_DIR}/dmrg-deps-install/arpack++
#            CMAKE_ARGS
#            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
#            -DARPACK_LIB=${ARPACK_LIBRARIES_GENERATOR}
#            -LAPACK_LIBRARIES=${LAPACK_LIBRARIES_GENERATOR}
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            INSTALL_COMMAND ""
            CONFIGURE_COMMAND ""
#            BUILD_COMMAND
#            ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include && find <INSTALL_DIR>/include -maxdepth 1 -type l -delete &&
#            ${CMAKE_COMMAND} -E create_symlink <SOURCE_DIR>/include <INSTALL_DIR>/include/arpack++
            BUILD_COMMAND
            ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/include <INSTALL_DIR>/include/arpack++
            DEPENDS lapacke::lapacke arpack::arpack gfortran::gfortran
            )

    ExternalProject_Get_Property(external_ARPACK++ INSTALL_DIR)
    add_library(arpack::arpack++ INTERFACE IMPORTED)
    target_link_libraries(arpack::arpack++ INTERFACE lapacke::lapacke arpack::arpack)
    target_include_directories(arpack::arpack++ SYSTEM INTERFACE ${INSTALL_DIR}/include)
    add_dependencies(arpack::arpack++ external_ARPACK++)
else()
    message(FATAL_ERROR "Dependency Arpack++ not found and DOWNLOAD_METHOD = ${DOWNLOAD_METHOD}")
endif()



