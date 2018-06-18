

enable_language(Fortran)
include(cmake-modules/FindGFortran.cmake)

find_package(GSL)
if(GSL_FOUND)
    message(STATUS "GSL FOUND IN SYSTEM: ${GSL_LIBRARIES}")
    add_library(GSL UNKNOWN IMPORTED)
else()
    message(STATUS "GSL will be installed into ${INSTALL_DIRECTORY}/gsl on first build.")
    include(ExternalProject)
    ExternalProject_Add(library_GSL
            URL      http://ftp.acc.umu.se/mirror/gnu.org/gnu/gsl/gsl-2.4.tar.gz
            PREFIX              "${INSTALL_DIRECTORY}/gsl"
            CONFIGURE_COMMAND
                cd <SOURCE_DIR> &&
                pwd &&
                ./configure --enable-silent-rules CFLAGS= --enable-shared=no --prefix=<INSTALL_DIR>
            BUILD_COMMAND
                cd <SOURCE_DIR> &&
                pwd &&
                ${CMAKE_MAKE_PROGRAM} --quiet --silent
            INSTALL_COMMAND
                cd <SOURCE_DIR> &&
                pwd &&
                ${CMAKE_MAKE_PROGRAM} --quiet install
            )

    ExternalProject_Get_Property(library_GSL INSTALL_DIR)
    add_library(GSL UNKNOWN IMPORTED)
#    add_library(GSLcblas UNKNOWN IMPORTED)
    add_dependencies(GSL library_GSL)
#    add_dependencies(GSLcblas library_GSL)
    set(GSL_LIBRARY        ${INSTALL_DIR}/lib/libgsl${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(GSL_CBLAS_LIBRARY  ${INSTALL_DIR}/lib/libgslcblas${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
    set(GSL_INCLUDE_DIRS ${INSTALL_DIR}/include)
endif()

set_target_properties(GSL PROPERTIES
        IMPORTED_LOCATION             "${GSL_LIBRARY}"
        INTERFACE_LINK_LIBRARIES      "${GSL_CBLAS_LIBRARY}"
        INCLUDE_DIRECTORIES           "${GSL_INCLUDE_DIRS}"
        )

target_link_libraries(${PROJECT_NAME} PRIVATE GSL)
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
