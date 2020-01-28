


include(cmake-modules/FindArpack++.cmake)
find_Arpackpp()

if (TARGET arpack::arpack++)
    message(STATUS "Arpack++ found")
elseif(NOT ${DOWNLOAD_METHOD} MATCHES "none")
    message(STATUS "Arpack++ will be installed into ${CMAKE_BINARY_DIR}/dmrg-deps-install/arpack++ on first build.")
    include(ExternalProject)
    ExternalProject_Add(external_ARPACK++
            URL         https://github.com/m-reuter/arpackpp/archive/2.3.0.tar.gz
            URL_MD5     1b09e35b6c44e118003922643b99978a
            PREFIX      ${CMAKE_BINARY_DIR}/dmrg-deps-build/arpack++
            INSTALL_DIR ${CMAKE_BINARY_DIR}/dmrg-deps-install/arpack++
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            INSTALL_COMMAND ""
            CONFIGURE_COMMAND ""
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



