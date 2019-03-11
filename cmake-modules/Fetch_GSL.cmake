find_package(GSL)

if(GSL_FOUND)
    if("${GSL_LIBRARY}" MATCHES "${CUSTOM_SUFFIX}")
        message("Found GSL with correct suffix")
    else()
        set(GSL_FOUND OFF)
        set(GSL_LIBRARY OFF)
        set(GSL_CBLAS_LIBRARY OFF)
        message("Found GSL with incorrect suffix")
        find_library(GSL_LIBRARY
                NAMES libgsl${CUSTOM_SUFFIX}
                PATHS ${GSL_LIBDIR}
                )
        find_library(GSL_CBLAS_LIBRARY
                NAMES libgslcblas${CUSTOM_SUFFIX}
                PATHS ${GSL_LIBDIR}
                )
        endif()
endif()

if (NOT GSL_LIBRARY OR NOT GSL_CBLAS_LIBRARY OR NOT GSL_INCLUDE_DIRS)
    # Try finding arpack as module library
    message(STATUS "SEARCHING FOR GSL IN LOADED MODULES")

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
    add_library(GSL INTERFACE)
else()
    message(STATUS "GSL will be installed into ${INSTALL_DIRECTORY}/gsl on first build.")
    include(ExternalProject)
    ExternalProject_Add(external_GSL
            URL      http://ftp.acc.umu.se/mirror/gnu.org/gnu/gsl/gsl-2.4.tar.gz
            PREFIX      ${BUILD_DIRECTORY}/gsl
            INSTALL_DIR ${INSTALL_DIRECTORY}/gsl
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

    ExternalProject_Get_Property(external_GSL INSTALL_DIR)
    add_library(GSL INTERFACE)
    add_dependencies(GSL external_GSL)
    set(GSL_LIBRARY        ${INSTALL_DIR}/lib/libgsl${CUSTOM_SUFFIX})
    set(GSL_CBLAS_LIBRARY  ${INSTALL_DIR}/lib/libgslcblas${CUSTOM_SUFFIX})
    set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
    set(GSL_INCLUDE_DIRS ${INSTALL_DIR}/include)
endif()

target_link_libraries(GSL INTERFACE ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
target_include_directories(GSL INTERFACE ${GSL_INCLUDE_DIRS})


#target_link_libraries(${PROJECT_NAME} PRIVATE GSL)
