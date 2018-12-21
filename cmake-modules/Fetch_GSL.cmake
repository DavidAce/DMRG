find_package(GSL)
if (NOT GSL_LIBRARY OR NOT GSL_CBLAS_LIBRARY OR NOT GSL_INCLUDE_DIRS)
    # Try finding arpack as module library
    message(STATUS "SEARCHING FOR ARPACK IN LOADED MODULES")

    find_library(GSL_LIBRARY
            NAMES libgsl${CUSTOM_SUFFIX}
            PATHS "$ENV{GSL_DIR}/lib"
            )
    find_library(GSL_CBLAS_LIBRARY
            NAMES libgslcblas${CUSTOM_SUFFIX}
            PATHS "$ENV{GSL_DIR}/lib"
            )
    find_path(GSL_INCLUDE_DIRS
            NAMES gsl_blas.h
            PATHS "$ENV{GSL_DIR}/include/gsl"
            )
endif()



if(GSL_LIBRARY AND GSL_CBLAS_LIBRARY AND GSL_INCLUDE_DIRS)
    set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
    message(STATUS "GSL FOUND IN SYSTEM: ${GSL_LIBRARIES}")
    add_library(GSL INTERFACE IMPORTED)
else()
    message(STATUS "GSL will be installed into ${INSTALL_DIRECTORY}/gsl on first build.")
    include(ExternalProject)
    ExternalProject_Add(library_GSL
            URL      http://ftp.acc.umu.se/mirror/gnu.org/gnu/gsl/gsl-2.4.tar.gz
            PREFIX              "${INSTALL_DIRECTORY}/gsl"
            CONFIGURE_COMMAND
                cd <SOURCE_DIR> &&
                pwd &&
                ./configure --enable-silent-rules CFLAGS= --enable-shared=yes --prefix=<INSTALL_DIR>
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
    add_library(GSL INTERFACE)
    add_dependencies(GSL library_GSL)
    set(GSL_LIBRARY        ${INSTALL_DIR}/lib/libgsl${CUSTOM_SUFFIX})
    set(GSL_CBLAS_LIBRARY  ${INSTALL_DIR}/lib/libgslcblas${CUSTOM_SUFFIX})
    set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
    set(GSL_INCLUDE_DIRS ${INSTALL_DIR}/include)
endif()

set_target_properties(GSL PROPERTIES
        INTERFACE_LINK_LIBRARIES      "${GSL_LIBRARY};${GSL_CBLAS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIRS}"
        )

#target_link_libraries(${PROJECT_NAME} PRIVATE GSL)
